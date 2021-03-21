/*
 * Vector Network Analyzer Library
 * Copyright © 2020, 2021 D Scott Guthridge <scott_guthridge@rompromity.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

// #define DEBUG

/*
 * DIVIDE_AND_ROUND_UP: divide a by b and round up
 */
#define DIVIDE_AND_ROUND_UP(a, b)	(((a) + (b) - 1) / (b))

/*
 * leakage_term_t: used to collect leakage terms
 */
typedef struct leakage_term {
    double complex lt_sum;	/* sum of the samples */
    double lt_sumsq;	        /* sum of the squares of the magnitudes */
    int lt_count;		/* count of samples */
} leakage_term_t;

/*
 * solve_state: iterate over vnacal_new_coefficient_t's
 */
typedef struct solve_state {
    /*
     * Common Data
     */
    vnacal_new_t *ss_vnp;		/* new calibration structure */
    leakage_term_t *ss_leakage_vector;	/* leakage term vector */
    int *ss_leakage_map;		/* leakage term map */
    double complex **ss_p_vector;	/* unknown parameter by index, findex */

    /*
     * Equation Iterator
     */
    enum {
	EI_START_SYSTEM,		/* start new system */
	EI_START_EQUATION,		/* start new equation */
	EI_TERM,			/* in term list */
	EI_END_EQUATION,		/* no remaining terms */
	EI_END_SYSTEM			/* no remaining equations */
    } ss_state;				/* iterator state */
    int ss_findex;			/* current frequency index */
    vnacal_new_system_t *ss_vnsp;	/* currrent system */
    vnacal_new_equation_t *ss_vnep;	/* current equation */
    vnacal_new_measurement_t *ss_vnmp;	/* current measurement */
    vnacal_new_coefficient_t *ss_vncp;	/* current term */
    double complex *ss_m_matrix;	/* matrix of cached m values */
    double complex *ss_s_matrix;	/* matrix of cached s values */
    uint32_t *ss_m_bitmap;		/* valid cells in ss_m_matrix */
    uint32_t *ss_s_bitmap;		/* valid cells in ss_s_matrix */
} solve_state_t;

/*
 * start_new_system: prepare to iterate over system
 *   @ssp:    pointer state structure
 *   @vnsp:   system of equations to iterate over
 */
static void start_new_system(solve_state_t *ssp, vnacal_new_system_t *vnsp)
{
    ssp->ss_state  = EI_START_SYSTEM;
    ssp->ss_vnsp   = vnsp;
    ssp->ss_vnep   = NULL;
    ssp->ss_vnmp   = NULL;
    ssp->ss_vncp   = NULL;
}

/*
 * get_equation: position to the next equation in the system
 *   @ssp: pointer state structure
 */
static bool get_equation(solve_state_t *ssp)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    vnacal_type_t type = VL_TYPE(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);

    for (;;) {
	vnacal_new_equation_t *vnep;
	vnacal_new_measurement_t *vnmp;
	int eq_row, eq_column, eq_cell;

	/*
	 * If we're starting a new system, set ss_vnep to the first
	 * equation.  If we're already started, advance to the next
	 * equation.  If we're at the end, keep returning false.
	 */
	switch (ssp->ss_state) {
	case EI_START_SYSTEM:
	    ssp->ss_vnep = ssp->ss_vnsp->vns_equation_list;
	    ssp->ss_state = EI_START_EQUATION;
	    break;

	case EI_START_EQUATION:
	case EI_TERM:
	case EI_END_EQUATION:
	    ssp->ss_vnep = ssp->ss_vnep->vne_next;
	    ssp->ss_state = EI_START_EQUATION;
	    break;

	case EI_END_SYSTEM:
	    return false;
	}
	if ((vnep = ssp->ss_vnep) == NULL) {
	    ssp->ss_state = EI_END_SYSTEM;
	    return false;
	}
	vnmp      = vnep->vne_vnmp;
	eq_row    = vnep->vne_row;
	eq_column = vnep->vne_column;
	eq_cell   = eq_row * s_columns + eq_column;

	/*
	 * If we're in a calibration type that doesn't handle off-diagonal
	 * leakage terms within the linear system and there's no signal
	 * path through the calibration standard for the current equation,
	 * then skip the equation.
	 */
	if (type != VNACAL_T16 && type != VNACAL_U16 && eq_row != eq_column &&
		!vnmp->vnm_reachability_matrix[eq_cell]) {
	    continue;
	}

	/*
	 * Success
	 */
	ssp->ss_vncp = NULL;
	break;
    }
    return true;
}

/*
 * get_coefficient: get the next coefficient in the system
 *   @ssp: pointer state structure
 */
static bool get_coefficient(solve_state_t *ssp)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    vnacal_new_equation_t *vnep = ssp->ss_vnep;
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_coefficient_t *vncp;

    /*
     * If we're starting a new equation, set ss_vncp to the first
     * term.  If we're already in the term list, advance to the next
     * term.  If we're at the end of the equation, keep returning false;
     */
    switch (ssp->ss_state) {
    case EI_START_SYSTEM:
	abort();			/* must call get_equation first */
	break;

    case EI_START_EQUATION:
	ssp->ss_vncp = ssp->ss_vnep->vne_coefficient_list;
	ssp->ss_state = EI_TERM;
	break;

    case EI_TERM:
	ssp->ss_vncp = ssp->ss_vncp->vnc_next;
	break;

    case EI_END_EQUATION:
    case EI_END_SYSTEM:
	return false;
    }
    if ((vncp = ssp->ss_vncp) == NULL) {
	ssp->ss_state = EI_END_EQUATION;
	return false;
    }

    /*
     * If we're starting a new measurement, clear the cached m and
     * s matrices.
     */
    if (vnep->vne_vnmp != ssp->ss_vnmp) {
	(void)memset((void *)ssp->ss_m_matrix, 0,
		m_rows * m_columns * sizeof(double complex));
	(void)memset((void *)ssp->ss_s_matrix, 0,
		s_rows * s_columns * sizeof(double complex));
	(void)memset((void *)ssp->ss_m_bitmap, 0,
		DIVIDE_AND_ROUND_UP(m_rows * m_columns, 32) *
		sizeof(uint32_t));
	(void)memset((void *)ssp->ss_s_bitmap, 0,
		DIVIDE_AND_ROUND_UP(s_rows * s_columns, 32) *
		sizeof(uint32_t));
	ssp->ss_vnmp = vnmp;
    }

    return true;
}


/*
 * get_m: get the current measurement
 *   @ssp: pointer state structure
 */
static double complex get_m(solve_state_t *ssp)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const int findex = ssp->ss_findex;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_equation_t *vnep = ssp->ss_vnep;
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
    const int m_cell = vncp->vnc_m_cell;

    assert(m_cell >= 0);
    if (!(ssp->ss_m_bitmap[m_cell / 32] & (1U << (m_cell % 32)))) {
	const int m_row    = m_cell / m_columns;
	const int m_column = m_cell % m_columns;
	leakage_term_t *leakage_vector = ssp->ss_leakage_vector;
	double complex leakage = 0.0;
	double complex value;

	if (m_row != m_column && leakage_vector != NULL) {
	    const int term = ssp->ss_leakage_map[m_cell];
	    const leakage_term_t *ltp = &leakage_vector[term];

	    if (ltp->lt_count > 0) {
		leakage = ltp->lt_sum / ltp->lt_count;
	    }
	}
	value = vnmp->vnm_m_matrix[m_cell][findex] - leakage;
	ssp->ss_m_matrix[m_cell] = value;
	ssp->ss_m_bitmap[m_cell / 32] |= (1U << (m_cell % 32));
    }
    return ssp->ss_m_matrix[m_cell];
}

/*
 * get_s: get the current s-parameter value
 *   @ssp: pointer state structure
 */
static double complex get_s(solve_state_t *ssp)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const int findex = ssp->ss_findex;
    vnacal_new_equation_t *vnep = ssp->ss_vnep;
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
    const int s_cell = vncp->vnc_s_cell;
    const double frequency = vnp->vn_frequency_vector[findex];

    assert(s_cell >= 0);
    if (!(ssp->ss_s_bitmap[s_cell / 32] & (1U << (s_cell % 32)))) {
	vnacal_new_parameter_t *vnprp;
	double complex value;

	vnprp = vnmp->vnm_s_matrix[vncp->vnc_s_cell];
	assert(vnprp != NULL);
	if (vnprp->vnpr_unknown) {
	    value = ssp->ss_p_vector[vnprp->vnpr_unknown_index][findex];
	} else {
	    value = _vnacal_get_parameter_value_i(vnprp->vnpr_parameter,
		    frequency);
	}
	ssp->ss_s_matrix[s_cell] = value;
	ssp->ss_s_bitmap[s_cell / 32] |= (1U << (s_cell % 32));
    }
    return ssp->ss_s_matrix[s_cell];
}

/*
 * analytic_solve: solve for the error terms where all s-parameters are known
 *   @ssp: pointer to state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 */
static int analytic_solve(solve_state_t *ssp, double complex *x_vector,
	int x_length)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int unknowns = vlp->vl_t_terms - 1;

    /*
     * For each system of equations...
     */
    assert(x_length == vnp->vn_systems * unknowns);
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	const int offset = sindex * unknowns;
	double complex a_matrix[vnp->vn_max_equations][unknowns];
	double complex b_vector[vnp->vn_max_equations];
	int eq_count = 0;

	/*
	 * Build the coefficient matrix (a) and right-hand side vector (b).
	 */
	for (int i = 0; i < vnp->vn_max_equations; ++i) {
	    for (int j = 0; j < unknowns; ++j) {
		a_matrix[i][j] = 0.0;
	    }
	    b_vector[i] = 0.0;
	}
	start_new_system(ssp, vnsp);
	while (get_equation(ssp)) {
	    while (get_coefficient(ssp)) {
		const vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
		double complex value = vncp->vnc_negative ? -1.0 : 1.0;
		int m_cell = vncp->vnc_m_cell;
		int s_cell = vncp->vnc_s_cell;

		if (m_cell >= 0) {
		    value *= get_m(ssp);
		}
		if (s_cell >= 0) {
		    value *= get_s(ssp);
		}
		if (vncp->vnc_coefficient == -1) {
		    b_vector[eq_count] = value;
		} else {
		    a_matrix[eq_count][vncp->vnc_coefficient] = value;
		}
	    }
	    ++eq_count;
	}

	/*
	 * Solve for the unknowns, using LU decomposition if a_matrix
	 * is square, or QR decomposition if the system is overdetermined.
	 */
	if (eq_count < unknowns) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "insufficient number of standards to solve error terms");
	    return -1;
	}
	if (eq_count == unknowns) {
	    double complex determinant;

	    determinant = _vnacommon_mldivide(&x_vector[offset],
		    &a_matrix[0][0], b_vector, unknowns, 1);
	    if (determinant == 0.0 || !isnormal(cabs(determinant))) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	} else {
	    int rank;

	    rank = _vnacommon_qrsolve(&x_vector[offset], &a_matrix[0][0],
		    b_vector, eq_count, unknowns, 1);
	    if (rank < unknowns) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	}
    }
    return 0;
}

/*
 * iterative_solve: solve for both error terms and unknown s-parameters
 *   @ssp: pointer to state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 *
 * This is an implementation of the algorithm described in H.Van Hamme
 * and M. Vanden Bossche, "Flexible vector network analyzer calibration
 * with accuracy bounds using an 8-term or a 16-term error correction
 * model," in IEEE Transactions on Microwave Theory and Techniques,
 * vol. 42, no. 6, pp. 976-987, June 1994, doi: 10.1109/22.293566.
 *
 * Instead of calculating the error bounds on the error parameters,
 * we use calc_rms_error to calculate error in terms of measurement error.
 */
static int iterative_solve(solve_state_t *ssp, double complex *x_vector,
	int x_length)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const int findex = ssp->ss_findex;
    vnacal_t *vcp = vnp->vn_vcp;
    const int p_length = vnp->vn_unknown_parameters;
    const int correlated = vnp->vn_correlated_parameters;
    vnacal_new_m_error_t *m_error_vector = vnp->vn_m_error_vector;
    const double frequency = vnp->vn_frequency_vector[findex];
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const bool is_t = VL_IS_T(vlp);
    double complex p0_vector[p_length];
    bool is_correlated_vector[p_length];
    int w_terms = 0;
    double noise = 0.0;
    double tracking = 0.0;
    double complex *w_term_vector = NULL;
    double *w_vector = NULL;
    int equations = 0;
    int p_equations;
    int j_rows;
    double previous_error = 0.0;
    int rv = -1;

    /*
     * Count equations not skipped due to unreachability and check
     * if enough equations were given.
     */
    assert(x_length == vnp->vn_systems * (vlp->vl_t_terms - 1));
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];

	start_new_system(ssp, vnsp);
	while (get_equation(ssp)) {
	    ++equations;
	}
    }
    if (equations < x_length + p_length) {
	_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: not enough "
		"standards given to solve the system");
	return -1;
    }
    p_equations = equations - x_length;
    j_rows = p_equations + correlated;

    /*
     * Copy the initial p_vector into p0_vector.
     */
    for (int i = 0; i < p_length; ++i) {
	p0_vector[i] = ssp->ss_p_vector[i][findex];
    }

    /*
     * Determine which unknown parameters are of correlated type.
     */
    (void)memset((void *)is_correlated_vector, 0, sizeof(is_correlated_vector));
    for (vnacal_new_parameter_t *vnprp = vnp->vn_unknown_parameter_list;
	    vnprp != NULL; vnprp = vnprp->vnpr_next_unknown) {
	if (vnprp->vnpr_parameter->vpmr_type == VNACAL_CORRELATED) {
	    assert(vnprp->vnpr_unknown_index >= 0);
	    assert(vnprp->vnpr_unknown_index < p_length);
	    is_correlated_vector[vnprp->vnpr_unknown_index] = true;
	}
    }

    /*
     * Iterate using Gauss-Newton to find the unknown parameters, p_vector.
     */
    for (int iteration = 0; /*EMPTY*/; ++iteration) {
	/* coefficient matrix for solving the error terms (x_vector) */
	double complex a_matrix[equations][x_length];

	/* right-hand side vector for solving the error terms */
	double complex b_vector[equations];

	/* orthogonal matrix from QR decomposition of a_matrix */
	double complex q_matrix[equations][equations];

	/* upper-triangular matrix from QR decompositin of a_matrix */
	double complex r_matrix[equations][x_length];

	/* jacobian matrix */
	double complex j_matrix[j_rows][p_length];

	/* right-hand side residual vector for Gauss-Newton */
	double complex k_vector[j_rows];

	/* difference vector from Gauss-Newton */
	double complex d_vector[p_length];

	/* current equation */
	int equation;

	/* rank of a_matrix */
	int rank;

	/* used to accumulate the RMS error */
	double current_error;

	/* used to scale d_vector for stability */
	double scale_factor;

	/*
	 * Build a_matrix and b_vector.
	 */
	for (int i = 0; i < equations; ++i) {
	    for (int j = 0; j < x_length; ++j) {
		a_matrix[i][j] = 0.0;
	    }
	    b_vector[i] = 0.0;
	}
	equation = 0;
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    start_new_system(ssp, vnsp);
	    while (get_equation(ssp)) {
		while (get_coefficient(ssp)) {
		    const vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
		    double complex v = vncp->vnc_negative ? -1.0 : 1.0;
		    int m_cell = vncp->vnc_m_cell;
		    int s_cell = vncp->vnc_s_cell;

		    if (m_cell >= 0) {
			v *= get_m(ssp);
		    }
		    if (w_vector != NULL) {
			v *= w_vector[equation];
		    }
		    if (s_cell >= 0) {
			v *= get_s(ssp);
		    }
		    if (vncp->vnc_coefficient == -1) {
			b_vector[equation] = v;
		    } else {
			a_matrix[equation][offset + vncp->vnc_coefficient] = v;
		    }
		}
		++equation;
	    }
	}
	assert(equation == equations);

	/*
	 * If we're using weighted equations and haven't yet done so,
	 * allocate the weight term vector and weight vector, and set
	 * up associated variables.
	 *
	 * Currently, we use weights only when we have to: when correlated
	 * unknown parameters are present and we need the weights to
	 * normalize the residuals between the correlated parameter
	 * equations and the other equations.  This may need revisiting.
	 */
	if (correlated != 0 && w_vector == NULL) {
	    w_terms = is_t ? VL_M_COLUMNS(vlp) : VL_M_ROWS(vlp);
	    noise = m_error_vector[findex].vnme_noise;
	    tracking = m_error_vector[findex].vnme_tracking;
	    w_term_vector = malloc(w_terms * sizeof(double complex));
	    if (w_term_vector == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"malloc: %s", strerror(errno));
		goto out;
	    }
	    if ((w_vector = malloc(equations * sizeof(double))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"malloc: %s", strerror(errno));
		goto out;
	    }
	    for (int i = 0; i < equations; ++i) {
		w_vector[i] = 1.0;
	    }
	}

	/*
	 * Decompose a_matrix into q_matrix and r_matrix, destroying
	 * a_matrix.
	 */
	rank = _vnacommon_qr(*a_matrix, *q_matrix, *r_matrix,
		equations, x_length);
	if (rank < x_length) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "singular linear system");
	    goto out;
	}

	/*
	 * Solve for x_vector.
	 */
	_vnacommon_qrsolve2(x_vector, *q_matrix, *r_matrix, b_vector,
		equations, x_length, 1);

	/*
	 * If there are no unknown parameters, we're done.
	 */
	if (p_length == 0)
	    break;

	/*
	 * Create the Jacobian matrix (j_matrix) and right-hand side
	 * (k_vector) for Gauss-Newton.
	 */
	for (int i = 0; i < j_rows; ++i) {
	    for (int j = 0; j < p_length; ++j) {
		j_matrix[i][j] = 0.0;
	    }
	    k_vector[i] = 0.0;
	}
	equation = 0;
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    start_new_system(ssp, vnsp);
	    while (get_equation(ssp)) {
		vnacal_new_equation_t *vnep = ssp->ss_vnep;
		vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;

		if (w_term_vector != NULL) {
		    for (int i = 0; i < w_terms; ++i) {
			w_term_vector[i] = 0.0;
		    }
		}
		while (get_coefficient(ssp)) {
		    const vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
		    int coefficient = vncp->vnc_coefficient;
		    int m_cell = vncp->vnc_m_cell;
		    int s_cell = vncp->vnc_s_cell;
		    vnacal_new_parameter_t *vnprp = NULL;

		    /*
		     * Apply the coefficient's contribution to the current
		     * row of the Jacobian matrix.
		     */
		    if (s_cell >= 0 && (vnprp =
				vnmp->vnm_s_matrix[s_cell])->vnpr_unknown) {
			double complex v = vncp->vnc_negative ? -1.0 : 1.0;
			int m_cell = vncp->vnc_m_cell;
			int unknown = vnprp->vnpr_unknown_index;

			if (m_cell >= 0) {
			    v *= get_m(ssp);
			}
			if (w_vector != NULL) {
			    v *= w_vector[equation];
			}
			assert(coefficient >= 0);
			v *= x_vector[offset + coefficient];
			for (int k = 0; k < p_equations; ++k) {
			    j_matrix[k][unknown] -=
				conj(q_matrix[equation][x_length + k]) * v;
			}
		    }

		    /*
		     * Accumulate the contribution to the weight vector.
		     */
		    if (w_term_vector != NULL && m_cell >= 0) {
			double t;
			double complex v;
			int i = is_t ? m_cell % VL_M_COLUMNS(vlp) :
			               m_cell / VL_M_COLUMNS(vlp);

			t = tracking * cabs(get_m(ssp));
			v = sqrt(noise * noise + t * t);
			if (vncp->vnc_negative) {
			    v *= -1.0;
			}
			if (coefficient >= 0) {
			    v *= x_vector[offset + coefficient];
			} else {
			    v *= -1.0;
			}
			if (s_cell >= 0) {
			    v *= cabs(get_s(ssp));
			}
			assert(i < w_terms);
			w_term_vector[i] += v;
		    }
		}

		/*
		 * Build the right-hand-side, k_vector.
		 */
		for (int k = 0; k < p_equations; ++k) {
		    k_vector[k] -= conj(q_matrix[equation][x_length + k]) *
			b_vector[equation];
		}

		/*
		 * Calculate the new weight.
		 */
		if (w_vector != NULL) {
		    double u = 0.0;

		    for (int i = 0; i < w_terms; ++i) {
			double complex v = w_term_vector[i];

			u += creal(v * conj(v));
		    }
		    u /= w_terms;
		    u = sqrt(u);
		    if (u < noise) {	/* avoid divide by zero */
			u = noise;
		    }
		    w_vector[equation] = 1.0 / u;
		}
		++equation;
	    }
	}
	assert(equation == equations);

	/*
	 * Add a row to j_matrix and k_vector for each correlated parameter.
	 */
	if (correlated != 0) {
	    int j_row = p_equations;
	    vnacal_new_parameter_t *vnprp1;

	    for (vnprp1 = vnp->vn_unknown_parameter_list; vnprp1 != NULL;
		    vnprp1 = vnprp1->vnpr_next_unknown) {
		vnacal_parameter_t *vpmrp1 = vnprp1->vnpr_parameter;
		vnacal_new_parameter_t *vnprp2;
		double coefficient;

		/*
		 * Skip if not a correlated parameter.
		 */
		if (vpmrp1->vpmr_type != VNACAL_CORRELATED)
		    continue;

		/*
		 * If the correlate is an unknown parameter, then place
		 * them in j_matrix with opposite signs to set them equal
		 * with one over sigma weight.	If the correlate is known,
		 * then put the known value into k_vector with one over
		 * sigma weight.
		 */
		coefficient = 1.0 / vpmrp1->vpmr_sigma;
		vnprp2 = vnprp1->vnpr_correlate;
		j_matrix[j_row][vnprp1->vnpr_unknown_index] = coefficient;
		if (vnprp2->vnpr_unknown) {
		    j_matrix[j_row][vnprp2->vnpr_unknown_index] = -coefficient;
		} else {
		    k_vector[j_row] = coefficient * vnprp2->vnpr_known_value;
		}
		++j_row;
	    }
	    assert(j_row == j_rows);
	}

	/*
	 * Solve the j_matrix, k_vector system to create d_vector, the
	 * correction to the unknown parameters, ss_p_vector.
	 */
	if (j_rows == p_length) {
	    double complex determinant;

	    determinant = _vnacommon_mldivide(d_vector,
		    &j_matrix[0][0], k_vector, p_length, 1);
	    if (determinant == 0.0 || !isnormal(cabs(determinant))) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	} else {
	    int rank;

	    rank = _vnacommon_qrsolve(d_vector, &j_matrix[0][0],
		    k_vector, j_rows, p_length, 1);
	    if (rank < p_length) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular linear system");
		return -1;
	    }
	}
#ifdef DEBUG
	(void)printf("# findex %d iteration %d\n", findex, iteration);
	(void)printf("# d =");
	for (int i = 0; i < p_length; ++i) {
	    (void)printf(" %f%+f",
		    creal(d_vector[i]),
		    cimag(d_vector[i]));
	}
	(void)printf("\n");
#endif /* DEBUG */

	/*
	 * If the norm of d_vector is more than half that of p_vector,
	 * scale the correction so that p_vector changes by only half
	 * its magnitude at a time, or by at most 0.5 if the magnitude of
	 * p_vector is less than one.  This greatly improves convergence
	 * by avoiding the situation where we continually overcorrect
	 * the elements of p_vector such that they form an unstable
	 * oscillation about zero.
	 */
	{
	    double sum_p_sq = 0.0;
	    double sum_d_sq = 0.0;

	    for (int i = 0; i < p_length; ++i) {
		double complex p, d;

		p = ssp->ss_p_vector[i][findex];
		sum_p_sq += creal(p * conj(p));
		d = d_vector[i];
		sum_d_sq += creal(d * conj(d));
	    }
	    if (sum_p_sq < 1.0) {
		sum_p_sq = 1.0;
	    }
	    if (sum_d_sq >= sum_p_sq * 0.25) {
		scale_factor = sqrt(sum_p_sq / sum_d_sq) * 0.5;
#ifdef DEBUG
		printf("# limiter: %f\n", scale_factor);
#endif /* DEBUG */
	    } else {
		scale_factor = 1.0;
	    }
	}
#ifdef DEBUG
	(void)printf("# p =");
#endif /* DEBUG */
	for (int i = 0; i < p_length; ++i) {
	    ssp->ss_p_vector[i][findex] += scale_factor * d_vector[i];
	    if (!is_correlated_vector[i] && p0_vector[i] != 0.0 &&
		    creal(ssp->ss_p_vector[i][findex] / p0_vector[i]) < 0.0) {
#ifdef DEBUG
		printf(" [applied flip]");
#endif /* DEBUG */
		ssp->ss_p_vector[i][findex] *= -1.0;
	    }
#ifdef DEBUG
	    (void)printf(" %f%+f",
		    creal(ssp->ss_p_vector[i][findex]),
		    cimag(ssp->ss_p_vector[i][findex]));
#endif /* DEBUG */
	}
#ifdef DEBUG
	(void)printf("\n");
#endif /* DEBUG */

	/*
	 * Calculate the RMS error of d_vector and stop if the
	 * error is no longer decreasing.
	 */
	current_error = 0.0;
	for (int i = 0; i < p_length; ++i) {
	    double e = cabs(d_vector[i]);
	    current_error += e * e;
	}
	current_error = sqrt(current_error / p_length);
#ifdef DEBUG
	(void)printf("# previous_error %e current_error %e\n",
		previous_error, current_error);
#endif /* DEBUG */
	if (current_error < 1.0e-5 && iteration > 0 &&
		current_error >= previous_error) {
#ifdef DEBUG
	    (void)printf("# converged %e\n", previous_error);
	    for (int i = 0; i < x_length; ++i) {
		printf("# x[%d] = %9.6f%+9.6fj\n",
			i, creal(x_vector[i]), cimag(x_vector[i]));
	    }
	    if (w_vector != NULL) {
		for (int i = 0; i < equations; ++i) {
		    printf("# w[%d] = %10f\n", i, w_vector[i]);
		}
	    }
#endif /* DEBUG */
	    break;
	}
	previous_error = current_error;

	/*
	 * Limit the number of iterations.
	 */
	if (iteration >= 50) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "system failed to converge at %e Hz", frequency);
	    goto out;
	}
#ifdef DEBUG
	if (previous_error > 1000.0) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "system diverged at %e Hz", frequency);
	    goto out;
	}
#endif /* DEBUG */
    }
    rv = 0;

out:
    free((void *)w_vector);
    free((void *)w_term_vector);
    return rv;
}

/*
 * calc_rms_error: calculate the RMS error of the solution, normalized to 1
 *   @ssp: pointer to state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 */
static double calc_rms_error(solve_state_t *ssp,
	double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const int findex = ssp->ss_findex;
    leakage_term_t *leakage_vector = ssp->ss_leakage_vector;
    const int correlated = vnp->vn_correlated_parameters;
    const vnacal_new_m_error_t *m_error_vector = vnp->vn_m_error_vector;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const bool is_t = VL_IS_T(vlp);
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    int w_terms = is_t ? VL_M_COLUMNS(vlp) : VL_M_ROWS(vlp);
    double noise = m_error_vector[findex].vnme_noise;
    double tracking = m_error_vector[findex].vnme_tracking;
    double complex w_term_vector[w_terms];
    double squared_error = 0.0;
    int count = 0;

    /*
     * Accumulate squared weighted residuals from the linear system.
     */
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	int offset = sindex * (vlp->vl_t_terms - 1);

	start_new_system(ssp, vnsp);
	while (get_equation(ssp)) {
	    double complex residual = 0.0;
	    double u = 0.0;

	    for (int i = 0; i < w_terms; ++i) {
		w_term_vector[i] = 0.0;
	    }
	    while (get_coefficient(ssp)) {
		const vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
		int coefficient = vncp->vnc_coefficient;
		int m_cell = vncp->vnc_m_cell;
		int s_cell = vncp->vnc_s_cell;
		double complex v = vncp->vnc_negative ? -1.0 : 1.0;

		if (coefficient >= 0) {
		    assert(offset + coefficient < x_length);
		    v *= x_vector[offset + coefficient];
		} else {
		    v *= -1.0;
		}
		if (s_cell >= 0) {
		    v *= get_s(ssp);
		}
		if (m_cell >= 0) {
		    double complex m = get_m(ssp);
		    double t;

		    int i = is_t ? m_cell % VL_M_COLUMNS(vlp) :
				   m_cell / VL_M_COLUMNS(vlp);

		    t = tracking * cabs(m);
		    assert(i < w_terms);
		    w_term_vector[i] += v * sqrt(noise * noise + t * t);
		    v *= m;
		}
		residual += v;
	    }
	    for (int i = 0; i < w_terms; ++i) {
		double complex v = w_term_vector[i];

		u += creal(v * conj(v));
	    }
	    u /= w_terms;
	    if (u < noise * noise) {	/* avoid divide by zero */
		u = noise * noise;
	    }
	    squared_error += creal(residual * conj(residual)) / u;
	    ++count;
	}
    }

    /*
     * Accumulate the error from correlated parameters.
     */
    if (correlated != 0) {
	vnacal_new_parameter_t *vnprp1;

	for (vnprp1 = vnp->vn_unknown_parameter_list; vnprp1 != NULL;
		vnprp1 = vnprp1->vnpr_next_unknown) {
	    vnacal_parameter_t *vpmrp1 = vnprp1->vnpr_parameter;

	    if (vpmrp1->vpmr_type == VNACAL_CORRELATED) {
		vnacal_new_parameter_t *vnprp2 = vnprp1->vnpr_correlate;
		double complex v;

		v = ssp->ss_p_vector[vnprp1->vnpr_unknown_index][findex];
		if (vnprp2->vnpr_unknown) {
		    v -= ssp->ss_p_vector[vnprp2->vnpr_unknown_index][findex];
		} else {
		    v -= vnprp2->vnpr_known_value;
		}
		squared_error += creal(v * conj(v)) /
		    (vpmrp1->vpmr_sigma * vpmrp1->vpmr_sigma);
		++count;
	    }
	}
    }

    /*
     * Accumulate variance from leakage parameter measurements.
     */
    if (leakage_vector != NULL) {
	for (int row = 0; row < m_rows; ++row) {
	    for (int column = 0; column < m_columns; ++column) {
		if (row != column) {
		    const int m_cell = row * m_columns + column;
		    const int term = ssp->ss_leakage_map[m_cell];
		    const leakage_term_t *ltp = &leakage_vector[term];
		    double value;

		    if (ltp->lt_count > 1) {
			double complex sum_x = ltp->lt_sum;
			const int n = ltp->lt_count;
			double n_mean_squared;

			n_mean_squared = creal(sum_x * conj(sum_x)) / n;
			value = (ltp->lt_sumsq - n_mean_squared);
			if (m_error_vector != NULL) {
			    double weight = 1.0 / (noise * noise +
				    n_mean_squared / n * tracking * tracking);

			    value *= weight;
			}
			if (value > 0.0) {
			    squared_error += value;
			}
			count += n - 1;
		    }
		}
	    }
	}
    }
    assert(!isnan(squared_error));
    assert(squared_error >= 0.0);
    return sqrt(squared_error / count);
}

/*
 * convert_ue14_to_e12: convert UE14 error terms to E12 error terms
 *   @e: input and output vector
 *   @vlp_in: input layout
 *   @vlp_out: output layout
 *   @ss.ss_leakage_map: map m_cell to external leakage term
 *
 * Do the matrix conversion:
 *   El = -Um^-1 Ui + El_in
 *   Er =  Um^-1
 *   Et =  Us - Ux Um^-1 Ui
 *   Em =  Ux Um^-1
 *
 *   but as a m_columns long sequence of independent m_rows x 1 systems
 *   with a row rotation.  At the same time, normalize Er and Et such that
 *   Et is the identity matrix.  Because Um, Ui, Ux, Us, Er and Et are all
 *   diagonal matrices, the conversion doesn't require a lot of computation.
 *   But we have to be careful with the indices given that each m_column
 *   represents an indepdent system.
 */
static int convert_ue14_to_e12(double complex *e,
	const vnacal_layout_t *vlp_in, const vnacal_layout_t *vlp_out,
	const int *_leakage_map)
{
    const int m_rows     = VL_M_ROWS(vlp_in);
    const int m_columns  = VL_M_COLUMNS(vlp_in);
    const int e12_terms  = VL_ERROR_TERMS(vlp_out);
    const double complex *el_in = &e[VL_EL_OFFSET(vlp_in)];
    double complex e_out[e12_terms];

    assert(VL_TYPE(vlp_in) == _VNACAL_E12_UE14);
    assert(VL_TYPE(vlp_out) == VNACAL_E12);
    assert(VL_M_ROWS(vlp_out) == m_rows);
    assert(VL_M_COLUMNS(vlp_out) == m_columns);
    for (int m_column = 0; m_column < m_columns; ++m_column) {
	const double complex *um = &e[VL_UM14_OFFSET(vlp_in, m_column)];
	const double complex *ui = &e[VL_UI14_OFFSET(vlp_in, m_column)];
	const double complex *ux = &e[VL_UX14_OFFSET(vlp_in, m_column)];
	const double complex *us = &e[VL_US14_OFFSET(vlp_in, m_column)];
	const double complex n = us[0] - ui[0] * ux[m_column] / um[m_column];
	double complex *el = &e_out[VL_EL12_OFFSET(vlp_out, m_column)];
	double complex *er = &e_out[VL_ER12_OFFSET(vlp_out, m_column)];
	double complex *em = &e_out[VL_EM12_OFFSET(vlp_out, m_column)];

	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    const int m_cell = m_row * m_columns + m_column;

	    /*
	     * Test for singular system.
	     */
	    if (um[m_row] == 0.0) {
		errno = EDOM;
		return -1;
	    }

	    /*
	     * Convert leakage term.
	     */
	    if (m_row == m_column) {
		el[m_row] = -ui[0] / um[m_column];
	    } else {
		el[m_row] = el_in[_leakage_map[m_cell]];
	    }

	    /*
	     * Convert reflection tracking term, normalizing
	     * to make the transmission tracking term 1.
	     */
	    er[m_row] = n / um[m_row];

	    /*
	     * Convert port match term.
	     */
	    em[m_row] = ux[m_row] / um[m_row];
	}
    }

    /*
     * Copy result.
     */
    (void)memcpy((void *)e, (void *)e_out,
	    e12_terms * sizeof(double complex));

    return 0;
}

/*
 * _vnacal_new_solve_internal: solve for the error parameters
 *   @vnp: pointer to vnacal_new_t structure
 */
int _vnacal_new_solve_internal(vnacal_new_t *vnp)
{
    vnacal_t *const vcp = vnp->vn_vcp;
    const int frequencies = vnp->vn_frequencies;
    const int unknown_parameters = vnp->vn_unknown_parameters;
    const vnacal_layout_t *const vlp_in = &vnp->vn_layout;
    const vnacal_layout_t *vlp_out = vlp_in;
    vnacal_layout_t vl_e12;
    const vnacal_type_t type_in = VL_TYPE(vlp_in);
    vnacal_type_t type_out = type_in;
    const int m_rows = VL_M_ROWS(vlp_in);
    const int m_columns = VL_M_COLUMNS(vlp_in);
    const int s_rows = VL_S_ROWS(vlp_in);
    const int s_columns = VL_S_COLUMNS(vlp_in);
    const int leakage_terms = VL_EL_TERMS(vlp_in);
    const int error_terms_in = VL_ERROR_TERMS(vlp_in);
    int error_terms_out = error_terms_in;
    solve_state_t ss;
    vnacal_calibration_t *calp = NULL;
    int rc = -1;

    /*
     * Init the state structure.
     */
    (void)memset((void *)&ss, 0, sizeof(solve_state_t));
    ss.ss_vnp = vnp;
    ss.ss_state = EI_END_SYSTEM;

    /*
     * Make sure the frequency vector was given.
     */
    if (!vnp->vn_frequencies_valid) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_solve: "
		"calibration frequency " "vector must be given");
	return -1;
    }

    /*
     * If there are leakage terms outside of the linear system, allocate
     * a vector of  structures to accumulate them, and a map from m_cell
     * to leakage term.
     */
    if (leakage_terms != 0) {
	int term = 0;

	if ((ss.ss_leakage_vector = calloc(leakage_terms,
			sizeof(leakage_term_t))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    goto out;
	}
	if ((ss.ss_leakage_map = calloc(m_rows * m_columns,
			sizeof(int))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    goto out;
	}
	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		const int m_cell = m_row * m_columns + m_column;

		if (m_row == m_column) {
		    ss.ss_leakage_map[m_cell] = -1;
		} else {
		    ss.ss_leakage_map[m_cell] = term++;
		}
	    }
	}
    }

    /*
     * If there are unknown parameters, allocate a vector, one entry per
     * frequency, of vectors of unknown parameter values.
     */
    if (unknown_parameters != 0) {
	if ((ss.ss_p_vector = calloc(unknown_parameters,
			sizeof(double complex *))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    goto out;
	}
	for (int i = 0; i < unknown_parameters; ++i) {
	    if ((ss.ss_p_vector[i] = calloc(frequencies,
			    sizeof(double complex))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s",
			strerror(errno));
		goto out;
	    }
	}
    }

    /*
     * Allocate matrices to cache calculated m and s values to save
     * from having to continually recalculate
     */
    if ((ss.ss_m_matrix = calloc(m_rows * m_columns,
		    sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto out;
    }
    if ((ss.ss_s_matrix = calloc(s_rows * s_columns,
		    sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto out;
    }
    if ((ss.ss_m_bitmap = calloc(DIVIDE_AND_ROUND_UP(m_rows * m_columns,
		    32), sizeof(uint32_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto out;
    }
    if ((ss.ss_s_bitmap = calloc(DIVIDE_AND_ROUND_UP(s_rows * s_columns,
		    32), sizeof(uint32_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto out;
    }

    /*
     * If the type is _VNACAL_E12_UE14, set up a different output layout
     * for automatically converting to VNACAL_E12.
     */
    if (type_in == _VNACAL_E12_UE14) {
	_vnacal_layout(&vl_e12, VNACAL_E12, m_rows, m_columns);
	vlp_out = &vl_e12;
	type_out = VNACAL_E12;
	error_terms_out = VL_ERROR_TERMS(vlp_out);
    }

    /*
     * Create the vnacal_calibration_t structure.
     */
    if ((calp = _vnacal_calibration_alloc(vcp, type_out, m_rows, m_columns,
		    frequencies, error_terms_out)) == NULL) {
	goto out;
    }
    (void)memcpy((void *)calp->cal_frequency_vector,
	    (void *)vnp->vn_frequency_vector,
	    frequencies * sizeof(double));
    calp->cal_z0 = vnp->vn_z0;

    /*
     * For each frequency, solve for the error parameters.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	const double frequency = vnp->vn_frequency_vector[findex];
	const int x_length = vnp->vn_systems * (vlp_in->vl_t_terms - 1);
	double complex x_vector[x_length];
	double complex e_vector[error_terms_out];
	int eterm_index = 0;

	/*
	 * Save the current frequency index in the state structure.
	 */
	ss.ss_findex = findex;

	/*
	 * Initialize the unknown parameter values.
	 */
	for (vnacal_new_parameter_t *vnprp = vnp->vn_unknown_parameter_list;
		vnprp != NULL; vnprp = vnprp->vnpr_next_unknown) {
	    ss.ss_p_vector[vnprp->vnpr_unknown_index][findex] =
		_vnacal_get_parameter_value_i(vnprp->vnpr_parameter, frequency);
	}

	/*
	 * If TE10, UE10 or UE14, compute leakage terms outside of the
	 * linear system.
	 */
	if (leakage_terms != 0) {
	    leakage_term_t *ltp;
	    vnacal_new_measurement_t *vnmp;

	    for (int term = 0; term < leakage_terms; ++term) {
		ltp = &ss.ss_leakage_vector[term];
		ltp->lt_sum   = 0.0;
		ltp->lt_sumsq = 0.0;
		ltp->lt_count = 0;
	    }
	    for (vnmp = vnp->vn_measurement_list; vnmp != NULL;
		    vnmp = vnmp->vnm_next) {
		for (int row = 0; row < m_rows; ++row) {
		    for (int column = 0; column < m_columns; ++column) {
			const int m_cell = row * m_columns + column;
			const int s_cell = row * s_columns + column;
			const int term = ss.ss_leakage_map[m_cell];
			double complex m;

			if (row == column) {
			    continue;
			}
			if (vnmp->vnm_m_matrix[m_cell] == NULL) {
			    continue;
			}
			if (vnmp->vnm_reachability_matrix[s_cell]) {
			    continue;
			}
			m = vnmp->vnm_m_matrix[m_cell][findex];
			ltp = &ss.ss_leakage_vector[term];
			ltp->lt_sum   += m;
			ltp->lt_sumsq += creal(m * conj(m));
			++ltp->lt_count;
		    }
		}
	    }
	    for (int term = 0; term < leakage_terms; ++term) {
		ltp = &ss.ss_leakage_vector[term];

		if (ltp->lt_count == 0) {
		    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			    "leakage term system is singular");
		    goto out;
		}
	    }
	}

	/*
	 * Solve the system.  If there are no unknown parameters, we can
	 * solve analytically using LU or QR decomposition.  If there
	 * are unknown parameters, then the system is non-linear and we
	 * have to use an iterative approach.
	 */
	if (unknown_parameters == 0) {
	    if (analytic_solve(&ss, x_vector, x_length) == -1) {
		goto out;
	    }
	} else {
	    if (iterative_solve(&ss, x_vector, x_length) == -1) {
		goto out;
	    }
	}

	/*
	 * If the mesurement error was given, calculate the RMS
	 * error of the solution and fail if it's too high.
	 */
	if (vnp->vn_m_error_vector != NULL) {
	    double error;

	    error = calc_rms_error(&ss, x_vector, x_length);
	    if (error > 6.0) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"too much error");
		goto out;
	    }
#ifdef DEBUG
	    (void)printf("# findex %d error %e\n", findex, error);
#endif /* DEBUG */
	}

	/*
	 * Copy from x_vector to e_vector, inserting the unity term.
	 */
        for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
            int unity_index = _vl_unity_offset(vlp_in, sindex);
	    const int offset = sindex * (vlp_in->vl_t_terms - 1);

            for (int k = 0; k < unity_index; ++k) {
                e_vector[eterm_index++] = x_vector[offset + k];
            }
            e_vector[eterm_index++] = 1.0;
            for (int k = unity_index; k < vlp_in->vl_t_terms - 1; ++k) {
                e_vector[eterm_index++] = x_vector[offset + k];
            }
	}

	/*
	 * Add the leakage terms.
	 */
	for (int term = 0; term < leakage_terms; ++term) {
	    leakage_term_t *ltp = &ss.ss_leakage_vector[term];
	    double complex el;

	    if (ltp->lt_count != 0) {
		el = ltp->lt_sum / ltp->lt_count;
	    } else {
		el = 0.0;
	    }
	    e_vector[eterm_index++] = el;
	}
	assert(eterm_index == error_terms_in);

	/*
	 * If _VNACAL_E12_UE14, convert to VNACAL_E12.
	 */
	if (type_in == _VNACAL_E12_UE14) {
	    rc = convert_ue14_to_e12(e_vector, vlp_in, vlp_out,
		    ss.ss_leakage_map);
	    if (rc == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular system");
		goto out;
	    }
	}

	/*
	 * Copy the error terms to the calibration structure.
	 */
	for (int term = 0; term < error_terms_out; ++term) {
	    calp->cal_error_term_vector[term][findex] = e_vector[term];
	}
    }

    /*
     * If we solved for unknown parameters, store them into the
     * corresponding parameter structures.
     */
    for (vnacal_new_parameter_t *vnprp = vnp->vn_unknown_parameter_list;
	    vnprp != NULL; vnprp = vnprp->vnpr_next_unknown) {
	vnacal_parameter_t *vpmrp = vnprp->vnpr_parameter;
	int index = vnprp->vnpr_unknown_index;

	assert(vpmrp->vpmr_type == VNACAL_UNKNOWN ||
	       vpmrp->vpmr_type == VNACAL_CORRELATED);
	free((void *)vpmrp->vpmr_gamma_vector);
	vpmrp->vpmr_gamma_vector = NULL;
	if (vpmrp->vpmr_frequencies != frequencies) {
	    free((void *)vpmrp->vpmr_frequency_vector);
	    vpmrp->vpmr_frequency_vector = NULL;
	    vpmrp->vpmr_frequencies = 0;
	    if ((vpmrp->vpmr_frequency_vector = calloc(frequencies,
			    sizeof(double))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"calloc: %s", strerror(errno));
		goto out;
	    }
	    vpmrp->vpmr_frequencies = frequencies;
	}
	(void)memcpy((void *)vpmrp->vpmr_frequency_vector,
		(void *)vnp->vn_frequency_vector,
		frequencies * sizeof(double));
	assert(ss.ss_p_vector[index] != NULL);
	vpmrp->vpmr_gamma_vector = ss.ss_p_vector[index];
	ss.ss_p_vector[index] = NULL;
    }

    /*
     * Save the solved error terms into the vnacal_new_t structure.
     */
    _vnacal_calibration_free(vnp->vn_calibration);
    vnp->vn_calibration = calp;
    calp = NULL;		/* ownership transferred to vnacal_new_t */
    rc = 0;

out:
    _vnacal_calibration_free(calp);
    free((void *)ss.ss_s_bitmap);
    free((void *)ss.ss_m_bitmap);
    free((void *)ss.ss_s_matrix);
    free((void *)ss.ss_m_matrix);
    if (ss.ss_p_vector != NULL) {
	for (int i = unknown_parameters - 1; i >= 0; --i) {
	    free((void *)ss.ss_p_vector[i]);
	}
	free((void *)ss.ss_p_vector);
    }
    free((void *)ss.ss_leakage_map);
    free((void *)ss.ss_leakage_vector);
    return rc;
}

/*
 * vnacal_new_solve: solve for the error parameters
 *   @vnp: pointer to vnacal_new_t structure
 */
int vnacal_new_solve(vnacal_new_t *vnp)
{
    if (vnp == NULL) {
	errno = EINVAL;
	return -1;
    }
    return _vnacal_new_solve_internal(vnp);
}
