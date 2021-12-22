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

#ifndef DBL_EPSILON
#define DBL_EPSILON	1.0e-9
#endif /* DBL_EPSILON */

//#define DEBUG

/*
 * PHI_INV: inverse of the golden ratio
 * PHI_INV2: inverse of the golden ratio squared
 */
#define PHI_INV		0.61803398874989484820
#define PHI_INV2	0.38196601125010515180

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
 * get_m_by_cell: get measurement by cell for the current equation
 */
static double complex get_m_by_cell(solve_state_t *ssp, int m_cell)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_columns = VL_M_COLUMNS(vlp);
    const int findex = ssp->ss_findex;
    vnacal_new_equation_t *vnep = ssp->ss_vnep;
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;

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
 * get_m: get measurement associated with the current coefficient (or -1)
 *   @ssp: pointer state structure
 */
static double complex get_m(solve_state_t *ssp)
{
    vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
    const int m_cell = vncp->vnc_m_cell;

    return get_m_by_cell(ssp, m_cell);
}

/*
 * get_s: get s-parameter value associated with the current coefficient (or -1)
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
 * calc_weights: compute w_vector from x_vector, p_vector
 *   @ssp: pointer to state structure
 *   @x_vector: current error parameter solution
 *   @w_vector: new weight vector
 *
 * TODO: we're not calculating weights according to the more elegant
 * solution given in the Van hamme paper.  We need to do some rework of
 * the way we represent equations in order to do it right. This is kind
 * of close, though.
 */
static void calc_weights(solve_state_t *ssp,
	const double complex *x_vector, double *w_vector)
{
    vnacal_new_t *vnp = ssp->ss_vnp;
    const int findex = ssp->ss_findex;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_cells = VL_M_COLUMNS(vlp) * VL_M_ROWS(vlp);
    vnacal_new_m_error_t *m_error_vector = vnp->vn_m_error_vector;
    const double noise = m_error_vector[findex].vnme_noise;
    const double tracking = m_error_vector[findex].vnme_tracking;
    double complex m_weight_vector[m_cells];
    int equation = 0;

    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	int offset = sindex * (vlp->vl_t_terms - 1);

	start_new_system(ssp, vnsp);
	while (get_equation(ssp)) {
	    double u = 0.0;

	    for (int m_cell = 0; m_cell < m_cells; ++m_cell) {
		m_weight_vector[m_cell] = 0.0;
	    }
	    while (get_coefficient(ssp)) {
		const vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
		int coefficient = vncp->vnc_coefficient;
		int m_cell = vncp->vnc_m_cell;
		int s_cell = vncp->vnc_s_cell;

		if (m_cell >= 0) {
		    double complex v = 1.0;

		    if (vncp->vnc_negative) {
			v *= -1.0;
		    }
		    if (s_cell >= 0) {
			v *= get_s(ssp);
		    }
		    if (coefficient >= 0) {
			v *= x_vector[offset + coefficient];
		    } else {
			v *= -1.0;
		    }
		    assert(m_cell < m_cells);
		    m_weight_vector[m_cell] += v;
		}
	    }

	    /*
	     * Calculate the new weight.
	     */
	    for (int m_cell = 0; m_cell < m_cells; ++m_cell) {
		double complex v = m_weight_vector[m_cell];

		if (v != 0.0) {
		    double complex m;
		    double temp;

		    /*
		     * Add the squared error contributed by m_cell.
		     */
		    m = get_m_by_cell(ssp, m_cell);
		    temp = noise * noise +
			tracking * tracking * creal(m * conj(m));
		    u += creal(v * conj(v)) * temp;
		}
	    }
	    u = sqrt(u);
	    if (u < noise) {	/* avoid divide by zero */
		u = noise;
	    }
	    w_vector[equation] = 1.0 / u;
	    ++equation;
	}
    }
}

/*
 * iterative_solve: solve for both error terms and unknown s-parameters
 *   @ssp: pointer to state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 *
 * This implementation is based on the algorithm described in H.Van Hamme
 * and M. Vanden Bossche, "Flexible vector network analyzer calibration
 * with accuracy bounds using an 8-term or a 16-term error correction
 * model," in IEEE Transactions on Microwave Theory and Techniques,
 * vol. 42, no. 6, pp. 976-987, June 1994, doi: 10.1109/22.293566.
 * There are a few differences, however.  For example, instead of
 * calculating the error bounds on the error parameters, we use
 * calc_rms_error to calculate error in terms of measurement error.
 */
static int iterative_solve(solve_state_t *ssp, double complex *x_vector,
	int x_length)
{
    /* pointer to vnacal_new_t structore */
    vnacal_new_t *vnp = ssp->ss_vnp;

    /* current frequency index */
    const int findex = ssp->ss_findex;

    /* current frequency */
    const double frequency = vnp->vn_frequency_vector[findex];

    /* pointer to vnacal_t structure */
    vnacal_t *vcp = vnp->vn_vcp;

    /* number of unknown (including correlated) parameters */
    const int p_length = vnp->vn_unknown_parameters;

    /* number of correlated parameters */
    const int correlated = vnp->vn_correlated_parameters;

    /* pointer to vnacal_layout_t structure */
    const vnacal_layout_t *vlp = &vnp->vn_layout;

    /* map indicating which of the unknown parameters are correlated type */
    bool is_correlated_vector[p_length];

    /* weight vector */
    double *w_vector = NULL;

    /* weight vector that was used to create best solution */
    double *best_w_vector = NULL;

    /* best error parameters */
    double complex best_x_vector[x_length];

    /* best unknown parameters */
    double complex best_p_vector[p_length];

    /* Gauss-Newton correction vector generated from the best solution */
    double complex best_d_vector[p_length];

    /* sum of squares of best_d_vector, initially infinite */
    double best_sum_d_squared = INFINITY;

    /* count of equations in the linear error term system */
    int equations = 0;

    /* count of "excess" equations used to solve for the unknown standards */
    int p_equations;

    /* count of rows in the Jacobian matrix */
    int j_rows;

    /* current number of iterations in the backtracking line search */
    int backtrack_count = 0;

    /* return status of this function */
    int rv = -1;

    /*
     * Count the equations not skipped due to unreachability and check
     * if enough were given to solve all the unknowns.
     */
    assert(x_length == vnp->vn_systems * (vlp->vl_t_terms - 1));
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];

	start_new_system(ssp, vnsp);
	while (get_equation(ssp)) {
	    ++equations;
	}
    }
    if (equations + correlated < x_length + p_length) {
	_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: not enough "
		"standards given to solve the system");
	return -1;
    }
    p_equations = equations - x_length;
    j_rows = p_equations + correlated;

    /*
     * Determine which unknown parameters are of correlated type.
     */
    (void)memset((void *)is_correlated_vector, 0, sizeof(is_correlated_vector));
    for (vnacal_new_parameter_t *vnprp = vnp->vn_unknown_parameter_list;
	    vnprp != NULL; vnprp = vnprp->vnpr_next_unknown) {
	assert(vnprp->vnpr_unknown_index >= 0);
	assert(vnprp->vnpr_unknown_index < p_length);
	if (vnprp->vnpr_parameter->vpmr_type == VNACAL_CORRELATED) {
	    is_correlated_vector[vnprp->vnpr_unknown_index] = true;
	}
    }

    /*
     * Iterate using Gauss-Newton to find the unknown parameters, ss_p_vector.
     */
    for (int iteration = 0; /*EMPTY*/; ++iteration) {
	/* coefficient matrix of the linear error term system */
	double complex a_matrix[equations][x_length];

	/* right-hand side of the linear error term system */
	double complex b_vector[equations];

	/* orthogonal matrix from QR decomposition of a_matrix */
	double complex q_matrix[equations][equations];

	/* upper-triangular matrix from QR decompositin of a_matrix */
	double complex r_matrix[equations][x_length];

	/* Jacobian matrix for Gauss-Newton */
	double complex j_matrix[j_rows][p_length];

	/* right-hand side residual vector for Gauss-Newton */
	double complex k_vector[j_rows];

	/* difference vector from Gauss-Newton */
	double complex d_vector[p_length];

	/* current equation index */
	int equation;

	/* rank of a_matrix */
	int rank;

	/* sum of squares of d_vector */
	double sum_d_squared;

	/*
	 * Build a_matrix and right-hand-side b_vector.  This linear
	 * system is built from the measurements of the calibration
	 * standards added to the vnacal_new_t structure via the
	 * vnacal_new_add_* functions.	It's used to solve for the error
	 * parameters, x_vector, given estimates of any unknown standards.
	 *
	 * Note that in calibration types other than T16 and U16, the
	 * leakage equations are excluded from the system.  For example,
	 * a double reflect standard in 2x2 T8 contributes only two
	 * equations intead of four.  In TE10 and UE10, the other two
	 * are used to compute leakage terms -- that's done outside of
	 * this function.
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

	    /*
	     * The start_new_system, get_equation and get_coefficient
	     * functions are abstract iterators that systematically
	     * go through the equations added via vnacal_new_add_*.
	     *
	     * In the case of UE14 (used to solve classic E12 SOLT),
	     * each column of the measurement matrix forms an independent
	     * linear system with its own separate error terms.  These
	     * independent systems, however, share the same unknown
	     * calibration parameters (p_vector), and for simplicity
	     * of solving them, we create one big (possibly sparse)
	     * matrix equation representing them all.
	     */
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
	 * Find the QR decomposition of a_matrix, creating q_matrix and
	 * r_matrix, destroying a_matrix.
	 *
	 * Conceptually, Q and R are partitioned as follows:
	 *
	 *   [ Q1 Q2 ] [ R
	 *               0 ]
	 *
	 * with dimensions:
	 *   Q1: equations x x_length
	 *   Q2: equations x (equations - x_length)
	 *   R:  x_length  x x_length
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
	 *   R x = Q^H b, where Q^H is the conjugate transpose of Q
	 */
	_vnacommon_qrsolve2(x_vector, *q_matrix, *r_matrix, b_vector,
		equations, x_length, 1);

	/*
	 * If measurement error was given (via vnacal_new_set_m_error),
	 * then we weight the equations in the system to compensate for
	 * uncertainty in the measurements.
	 *
	 * If we haven't already done so, allocate w_vector and
	 * best_w_vector, compute weights based on the initial guesses
	 * of ss_p_vector, and restart the loop from the top.
	 */
	if (vnp->vn_m_error_vector != NULL && w_vector == NULL) {
	    if ((w_vector = malloc(equations * sizeof(double))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"malloc: %s", strerror(errno));
		goto out;
	    }
	    if (p_length != 0) {
		if ((best_w_vector = calloc(equations,
				sizeof(double))) == NULL) {
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "calloc: %s", strerror(errno));
		    goto out;
		}
	    }
	    for (int i = 0; i < equations; ++i) {
		w_vector[i] = 1.0;
	    }
	    calc_weights(ssp, x_vector, w_vector);
	    --iteration;
	    continue;
	}

	/*
	 * If there are no unknown parameters, we're done.
	 */
	if (p_length == 0)
	    goto success;

	/*
	 * At this point, we know that the system is nonlinear.
	 * We have two sets of variables to solve, the error terms,
	 * x_vector, and the unknown calibration parameters, p_vector.
	 * The a_matrix depends on p_vector; consequently, the
	 * system A x = b contains products of p and x variables
	 * and is nonlinear.  It is, however, a separable nonlinear
	 * least squares problem that can be solved using the variable
	 * projection method as described by Golub and LeVeque, 1979
	 * http://faculty.washington.edu/rjl/pubs/GolubLeVeque1979/
	 * GolubLeVeque1979.pdf
	 *
	 * Using this method, we make an initial guess for p, solve x as
	 * a linear system, project the remaining equations into a new
	 * space that lets us construct the Jacobian in terms of p only,
	 * use Gauss-Newton to improve our estimate of p and repeat from
	 * the solve for x step until we have suitable convergence.
	 *
	 * Following is a brief derivation of the variable projection method.
	 * 
	 * Our goal is to minimize the system A(p) x = b in a least-squares
	 * sense, where A(p) is matrix valued function of vector p, b is a
	 * known vector, and x and p are the unknown vectors we need to find
	 * to minimize:
	 *
	 *     || b - A(p) x ||^2
	 * 
	 * There must be an orthogonal matrix Q that diagonalizes A to R.
	 * Both of the new resulting matrices still depend on p.
	 *
	 *   A(p) = Q(p) R(p)
	 *
	 * Partition Q(p) and R(p) as follows:
	 *
	 *   A(p) = [ Q1(p) Q2(p) ] [ R1(p) ]
	 *                          [   0   ]
	 *        = Q1(p) R1(p)
	 *
	 * It also follows that:
	 *
	 *   Q2(p)^H A(p) = 0
	 *
	 * where ^H is the conjugate transpose.
	 *
	 * Solve A(p) x = b for x in a least-squares sense:
	 *
	 *   A(p)        x = b
	 *   Q1(p) R1(p) x = b
	 *   R1(p)       x = Q1(p)^H b
	 *               x = R1(p)^-1 Q1(p)^H b
	 *
	 * From the invariance of the 2-norm under orthagonal
	 * transformations, we can multiply the inside by Q^H
	 * without changing the norm:
	 *
	 *     || b - A(p) x ||^2
	 *
	 *   = || Q(p)^H (b - A(p) x) ||^2
	 *
	 *   = || Q1(p)^H b - Q1(p)^H A(p) x ||^2
	 *     || Q2(p)^H b - Q2(p)^H A(p) x ||^2
	 *
	 * But Q1(p)^H A(p) = R(p), and Q2(p)^H A(p) = 0, so
	 *
	 *   = || Q1(p)^H b - R(p) x ||^2
	 *     || Q2(p)^H b - 0      ||
	 *
	 * and because R(p) x = Q1(p)^H b from above, Q1(p)^H b - R(p) x = 0
	 *
	 *   = || 0         ||
	 *     || Q2(p)^H b ||
	 *
	 * so the residual we need to minimize is simply:
	 *
	 *   Q2(p)^H b
	 *
	 * For Gauss-Newton, we need the Jacobian of the residual above
	 * with respect to p.  Note that our choice of using tm11 or
	 * um11 for the unity term in the T or U error parameters,
	 * respectively, ensures that b never depends on p.
	 *
	 * Recall from above that Q2(p)^H A(p) = 0.  If we take the
	 * derivative, then from the product rule, we get:
	 *
	 *   Q2(p)^H' A(p) +  Q2(p)^H A(p)' = 0
	 *
	 * Re-arranging:
	 *
	 *   Q2(p)^H' A(p) = -Q2(p)^H A(p)'
	 *
	 * Using A(p) = Q1(p) R(p):
	 *
	 *   Q2(p)^H' Q1(p) R(p) = -Q2(p)^H A(p)'
	 *
	 * Multiply on the right by R(p)^-1 Q1(p)^H b:
	 *
	 *   Q2(p)^H' Q1(p) Q1(p)^H b = -Q2(p)^H d A(p)' R(p)^-1 Q1(p)^H b
	 *
	 * We'd really like Q1(p) Q1(p)^H to cancel, but they don't in
	 * this direction.  But Kaufman 1975 suggests the approximation:
	 *
	 *   Q2(p)^H' ≈ -Q2(p)^H A(p)' A(p)^+
	 *   where A(p)^+ is the pseudoinverse of A(p), or R(p)^-1 Q(p)^H
	 *
	 * Using the approximation, we can treat Q1(p) Q1(p)^H as if they
	 * do cancel:
	 *
	 *   Q2(p)^H' b ≈ -Q2(p)^H A(p)' R(p)^-1 Q1(p)^H b
	 *
	 * From above, R(p)^-1 Q1(p)^H b = x:
	 *
	 *   Q2(p)^H' b ≈ -Q2(p)^H A(p)' x
	 *
	 * We can easily find A(p)' since it's just the coefficients
	 * of A that contain the given p.  The result is our Jacobian
	 * (j_matrix):
	 *
	 *   J(p) ≈ -Q2(p)^H A(p)' x
	 *
	 * And the right hand side residual for Gauss-Newton is:
	 *
	 *   k(p) =  Q2(p)^H b
	 *
	 * To find the correction in p, we use QR decomposition to solve:
	 *
	 *   J(p) d = k(p)
	 *
	 * and apply the correction:
	 *
	 *   p -= d
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

		while (get_coefficient(ssp)) {
		    const vnacal_new_coefficient_t *vncp = ssp->ss_vncp;
		    int coefficient = vncp->vnc_coefficient;
		    int s_cell = vncp->vnc_s_cell;
		    vnacal_new_parameter_t *vnprp = NULL;

		    /*
		     * Apply this coefficient's contribution to the
		     * current row of the Jacobian matrix.  What we're
		     * doing here is computing -Q2(p)^H A'(p) x, but
		     * doing the the first matrix multiplication with
		     * loop nesting inverted from the usual order so
		     * that we can go row by row through A.
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
		}

		/*
		 * Build the right-hand-side vector of residuals, k_vector.
		 */
		for (int k = 0; k < p_equations; ++k) {
		    k_vector[k] -= conj(q_matrix[equation][x_length + k]) *
			b_vector[equation];
		}
		++equation;
	    }
	}
	assert(equation == equations);

	/*
	 * Add an additional row to j_matrix and k_vector for each
	 * correlated parameter.
	 *
	 * When the parameter is correlated with a constant parameter,
	 * we have an equation of the form:
	 *
	 *   1/sigma p_i = 1/sigma constant
	 *
	 * when a correlated parameter is correlated with another unknown
	 * parameter, we have an equation of the form:
	 *
	 *   1/sigma p_i - 1/sigma p_j = 0
	 *
	 * We store the Jacobian of the coefficient matrix (just the
	 * 1/sigma terms) into j_matrix and the residuals into k_matrix.
	 * In the terminology of the Van Hamme paper, the elements in
	 * j_matrix are E matrix, and the elements of
	 * k_vector are the residuals E*p - f.
	 */
	if (correlated != 0) {
	    int j_row = p_equations;
	    vnacal_new_parameter_t *vnprp1;

	    for (vnprp1 = vnp->vn_unknown_parameter_list; vnprp1 != NULL;
		    vnprp1 = vnprp1->vnpr_next_unknown) {
		vnacal_parameter_t *vpmrp1 = vnprp1->vnpr_parameter;
		vnacal_new_parameter_t *vnprp2;
		int pindex1;
		double coefficient;

		/*
		 * Skip if not a correlated parameter.
		 */
		if (vpmrp1->vpmr_type != VNACAL_CORRELATED)
		    continue;

		/*
		 * Place the partial derivative of the correlated
		 * parameter into j_matrix and it's contribution to the
		 * residual into k_vector, both weighted by sigma^-1.
		 * If the correlate is an unknown parameter, also place
		 * its partial derivative into j_matrix with opposite
		 * sign, effectively setting them equal.  Known or not,
		 * subract the contribution to the residual from k_vector.
		 */
		coefficient = 1.0 / _vnacal_get_correlated_sigma(vpmrp1,
			frequency);
		vnprp2 = vnprp1->vnpr_correlate;
		pindex1 = vnprp1->vnpr_unknown_index;
		j_matrix[j_row][pindex1] = coefficient;	/* partial derivative */
		k_vector[j_row] += coefficient *
		    ssp->ss_p_vector[pindex1][findex];	/* contr. to residual */
		if (vnprp2->vnpr_unknown) {
		    int pindex2 = vnprp2->vnpr_unknown_index;
		    j_matrix[j_row][pindex2] = -coefficient; /* partial drvtv */
		    k_vector[j_row] -= coefficient *
			ssp->ss_p_vector[pindex2][findex];   /* resid */
		} else { /* known parameter value */
		    k_vector[j_row] -= coefficient *	/* contr. to resid */
			_vnacal_get_parameter_value_i(vnprp2->vnpr_parameter,
				frequency);
		}
		++j_row;
	    }
	    assert(j_row == j_rows);
	}

	/*
	 * Solve the j_matrix, k_vector system to create d_vector, the
	 * Gauss-Newton correction to ss_p_vector.
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
	for (int i = 0; i < p_length; ++i) {
	    (void)printf("# d[%2d] = %13.6e %+13.6ej\n", i,
		    creal(d_vector[i]),
		    cimag(d_vector[i]));
	}
#endif /* DEBUG */

	/*
	 * Calculate the squared magnitude of d_vector.
	 */
	sum_d_squared = 0.0;
	for (int i = 0; i < p_length; ++i) {
	    sum_d_squared += creal(d_vector[i] * conj(d_vector[i]));
	}
#ifdef DEBUG
	(void)printf("# sum_d_squared      %13.6e\n", sum_d_squared);
	(void)printf("# best_sum_d_squared %13.6e\n", best_sum_d_squared);
#endif /* DEBUG */

	/*
	 * If we have the best solution so far (or the first solution),
	 * add the new correction to ss_p_vector and remember the solution.
	 */
	if (sum_d_squared < best_sum_d_squared) {
	    double sum_p_squared = 0.0;

#ifdef DEBUG
	    (void)printf("# best\n");
#endif /* DEBUG */
	    /*
	     * Limit the magnitude of d_vector to keep it smaller than
	     * the magnitude of ss_p_vector (or smaller than one if
	     * ss_p_vector is less than one).  This improves stability
	     * at the cost of slowing convergence.  It makes it less
	     * likely that we jump into an adjacent basin of attraction.
	     * We use one over the golden ratio as the maximum norm of
	     * d_vector relative to ss_p_vector (or one).
	     */
	    for (int i = 0; i < p_length; ++i) {
		double complex p = ssp->ss_p_vector[i][findex];

		sum_p_squared += creal(p * conj(p));
	    }
	    if (sum_p_squared < 1.0) {	/* if less than 1, make it 1 */
		sum_p_squared = 1.0;
	    }
	    if (sum_d_squared > sum_p_squared * PHI_INV2) {
		double scale = sqrt(sum_p_squared / sum_d_squared) * PHI_INV;

#ifdef DEBUG
		(void)printf("# scaling d_vector by: %f\n", scale);
#endif /* DEBUG */
		for (int i = 0; i < p_length; ++i) {
		    d_vector[i] *= scale;
		}
	    }

	    /*
	     * Remember this solution.
	     */
	    (void)memcpy((void *)best_x_vector, (void *)x_vector,
		    x_length * sizeof(double complex));
	    for (int i = 0; i < p_length; ++i) {
		best_p_vector[i] = ssp->ss_p_vector[i][findex];
		best_d_vector[i] = d_vector[i];
	    }
	    if (w_vector != NULL) {
		(void)memcpy((void *)best_w_vector, (void *)w_vector,
			equations * sizeof(double));
	    }
	    best_sum_d_squared = sum_d_squared;

	    /*
	     * Update the weight vector.
	     */
	    if (w_vector != NULL) {
		calc_weights(ssp, x_vector, w_vector);
#ifdef DEBUG
		for (int i = 0; i < equations; ++i) {
		    (void)printf("# w[%2d] = %f\n", i, w_vector[i]);
		}
#endif /* DEBUG */
	    }

	    /*
	     * Apply d_vector to ss_p_vector.
	     */
	    for (int i = 0; i < p_length; ++i) {
		ssp->ss_p_vector[i][findex] += d_vector[i];
	    }
#ifdef DEBUG
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("# p[%2d] = %13.6e %+13.6ej\n", i,
			creal(ssp->ss_p_vector[i][findex]),
			cimag(ssp->ss_p_vector[i][findex]));
	    }
#endif /* DEBUG */
	    backtrack_count = 0;

	/*
	 * If the squared error is nearly zero, stop.
	 */
	} else if (cabs(sum_d_squared) / p_length <= DBL_EPSILON) {
#ifdef DEBUG
	    (void)printf("# stop: converged\n");
#endif /* DEBUG */
	    break;

	/*
	 * The new solution is worse: we must have over-corrected.  Use a
	 * backtracking line search that keeps dividing d_vector in half
	 * and retrying from the best solution.
	 */
	} else {
	    if (++backtrack_count > 6) {
#ifdef DEBUG
	    (void)printf("# stop: backtrack count\n");
#endif /* DEBUG */
		break;
	    }
#ifdef DEBUG
	    (void)printf("# retry with half d_vector\n");
#endif /* DEBUG */
	    for (int i = 0; i < p_length; ++i) {
		best_d_vector[i] *= 0.5;
#ifdef DEBUG
	    (void)printf("# half-d[%2d] = %13.6e %+13.6ej\n", i,
		    creal(best_d_vector[i]),
		    cimag(best_d_vector[i]));
#endif /* DEBUG */
		ssp->ss_p_vector[i][findex] =
		    best_p_vector[i] + best_d_vector[i];
	    }
#ifdef DEBUG
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("# p[%2d] = %13.6e %+13.6ej\n", i,
			creal(ssp->ss_p_vector[i][findex]),
			cimag(ssp->ss_p_vector[i][findex]));
	    }
#endif /* DEBUG */
	    if (w_vector != NULL) {
		(void)memcpy((void *)w_vector, (void *)best_w_vector,
			equations * sizeof(double));
		calc_weights(ssp, x_vector, w_vector);
#ifdef DEBUG
		for (int i = 0; i < equations; ++i) {
		    (void)printf("# w[%2d] = %f\n", i, w_vector[i]);
		}
#endif /* DEBUG */
	    }
	}

	/*
	 * Limit the number of iterations.
	 *
	 * TODO: instead of failing here, just return what we have so
	 * far and let calc_rms_error check if it's close enough.
	 */
	if (iteration >= 50) {
	    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
		    "system failed to converge at %e Hz", frequency);
	    goto out;
	}
    }

    /*
     * Load the best solution.
     */
    (void)memcpy((void *)x_vector, (void *)best_x_vector,
	    x_length * sizeof(double complex));
    for (int i = 0; i < p_length; ++i) {
	ssp->ss_p_vector[i][findex] = best_p_vector[i];
    }
success:
    rv = 0;
    /*FALLTHROUGH*/

out:
    free((void *)w_vector);
    free((void *)best_w_vector);
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
    const double frequency = vnp->vn_frequency_vector[findex];
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
     * TODO: need to re-work the way we're weighting the equations
     * Consider basing the weights on the initial guesses only to avoid
     * the situation where the choice of weights effectively eliminates
     * equations and makes the system underdetermined.
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
		double sigma;

		v = ssp->ss_p_vector[vnprp1->vnpr_unknown_index][findex];
		if (vnprp2->vnpr_unknown) {
		    v -= ssp->ss_p_vector[vnprp2->vnpr_unknown_index][findex];
		} else {
		    v -= _vnacal_get_parameter_value_i(vnprp2->vnpr_parameter,
				frequency);
		}
		sigma = _vnacal_get_correlated_sigma(vpmrp1, frequency);
		squared_error += creal(v * conj(v)) / (sigma * sigma);
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
	 * Solve the system.  If there are no unknown parameters and
	 * no errors given on the measurements, we can use a simpler
	 * analytic solution using LU or QR decomposition.  If there
	 * are unknown parameters, measurement errors were given, or
	 * then the system is nonlinear then use the iterative solver.
	 */
	if (unknown_parameters == 0 && vnp->vn_m_error_vector == NULL) {
	
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
#ifdef DEBUG
	    (void)printf("# findex %d error %13.6e\n", findex, error);
#endif /* DEBUG */
	    if (error > 6.0) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"too much error");
		goto out;
	    }
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
