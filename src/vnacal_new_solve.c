/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include "vnacal_new_internal.h"

/* #define DEBUG */

/*
 * PHI_INV: inverse of the golden ratio
 * PHI_INV2: inverse of the golden ratio squared
 */
#define PHI_INV		0.61803398874989484820
#define PHI_INV2	0.38196601125010515180

/*
 * vnacal_new_leakage_term_t: leakage term outside the linear system
 */
typedef struct vnacal_new_leakage_term {
    /* sum of the samples */
    double complex vnlt_sum;

    /* sum of squared magnitudes of the samples */
    double vnlt_sumsq;

    /* count of accumulated samples */
    int vnlt_count;

} vnacal_new_leakage_term_t;

/*
 * vnacal_new_iterator_state_t: coefficient iterator states
 */
typedef enum {
    VNACAL_NI_INIT,			/* not started */
    VNACAL_NI_SYSTEM,			/* in system */
    VNACAL_NI_EQUATION,			/* in equation */
    VNACAL_NI_COEFFICIENT,		/* in coefficient list */
    VNACAL_NI_END_COEFFICIENTS,		/* no remaining coefficients */
    VNACAL_NI_END_EQUATIONS		/* no remaining equations */
} vnacal_new_iterator_state_t;

/*
 * vnacal_new_ms_matrices_t: a measured standard for solve
 */
typedef struct vnacal_new_ms_matrices {
    /* correponding measured standard */
    struct vnacal_new_measurement *vnmm_vnmp;

    /* matrix of measured values for the current frequency */
    double complex *vnmm_m_matrix;

    /* matrix of values of the standard for the current frequency */
    double complex *vnmm_s_matrix;

} vnacal_new_ms_matrices_t;

/*
 * vnacal_new_solve_state: iterate over vnacal_new_coefficient_t's
 */
typedef struct vnacal_new_solve_state {
    /* new calibration structure */
    vnacal_new_t *vnss_vnp;

    /* current frequency index */
    int vnss_findex;

    /* vector of structures correponding to each measured standard */
    vnacal_new_ms_matrices_t *vnss_ms_matrices;

    /* serialized matrix of pointers to leakage term structures */
    vnacal_new_leakage_term_t **vnss_leakage_matrix;

    /* vector of vector of unknown parameter values [index][findex] */
    double complex **vnss_p_vector;

    /* equation iterator state */
    vnacal_new_iterator_state_t vnss_iterator_state;

    /* current system in iterator */
    int vnss_sindex;

    /* current equation in iterator */
    vnacal_new_equation_t *vnss_vnep;

    /* current coefficient in iterator */
    vnacal_new_coefficient_t *vnss_vncp;

} vnacal_new_solve_state_t;

static void vs_free(vnacal_new_solve_state_t *vnssp);	/* forward */

/*
 * vs_init: initialize the solve state structure
 *   @vnssp: solve state structure
 *   @vnp:   associated vnacal_new_t structure
 */
static int vs_init(vnacal_new_solve_state_t *vnssp, vnacal_new_t *vnp)
{
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);

    /*
     * Start initializing the solve state structure.
     */
    (void)memset((void *)vnssp, 0, sizeof(*vnssp));
    vnssp->vnss_vnp = vnp;
    vnssp->vnss_findex = -1;

    /*
     * Allocate a vector of vnacal_new_ms_matrices_t structures,
     * each corresponding to the measured standard with same index.
     */
    if ((vnssp->vnss_ms_matrices = calloc(vnp->vn_measurement_count,
		    sizeof(vnacal_new_ms_matrices_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	return -1;
    }
    for (vnacal_new_measurement_t *vnmp = vnp->vn_measurement_list;
	    vnmp != NULL; vnmp = vnmp->vnm_next) {
	int index = vnmp->vnm_index;
	vnacal_new_ms_matrices_t *vnmmp = &vnssp->vnss_ms_matrices[index];

	/*
	 * Set the back pointer.
	 */
	vnmmp->vnmm_vnmp = vnmp;

	/*
	 * Allocate the temporary M and S matrices.
	 */
	if ((vnmmp->vnmm_m_matrix = calloc(m_rows * m_columns,
			sizeof(double complex))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    vs_free(vnssp);
	    return -1;
	}
	if ((vnmmp->vnmm_s_matrix = calloc(s_rows * s_columns,
			sizeof(double complex))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    vs_free(vnssp);
	    return -1;
	}
    }

    /*
     * If the error term type has leakage terms outside of the linear
     * system, allocate a matrix of pointers to vnacal_new_leakage_term_t
     * structures and populate the off-diagonal elements.
     */
    if (VNACAL_HAS_OUTSIDE_LEAKAGE_TERMS(vlp->vl_type)) {
	if ((vnssp->vnss_leakage_matrix = calloc(m_rows * m_columns,
			sizeof(vnacal_new_leakage_term_t *))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    vs_free(vnssp);
	    return -1;
	}
	for (int r = 0; r < m_rows; ++r) {
	    for (int c = 0; c < m_columns; ++c) {
		if (r != c) {
		    const int cell = r * m_columns + c;
		    vnacal_new_leakage_term_t *vnltp;

		    vnltp = malloc(sizeof(vnacal_new_leakage_term_t));
		    if (vnltp == NULL) {
			_vnacal_error(vcp, VNAERR_SYSTEM,
				"calloc: %s", strerror(errno));
			vs_free(vnssp);
			return -1;
		    }
		    (void)memset((void *)vnltp, 0, sizeof(*vnltp));
		    vnssp->vnss_leakage_matrix[cell] = vnltp;
		}
	    }
	}
    }

    /*
     * If there are unknown parameters, allocate a vector, one entry per
     * frequency, of vectors of unknown parameter values.
     */
    if (vnp->vn_unknown_parameters != 0) {
	if ((vnssp->vnss_p_vector = calloc(vnp->vn_unknown_parameters,
			sizeof(double complex *))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    vs_free(vnssp);
	    return -1;
	}
	for (int i = 0; i < vnp->vn_unknown_parameters; ++i) {
	    if ((vnssp->vnss_p_vector[i] = calloc(vnp->vn_frequencies,
			    sizeof(double complex))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s",
			strerror(errno));
		vs_free(vnssp);
		return -1;
	    }
	}
    }

    /*
     * Init the coefficient iterator state.
     */
    vnssp->vnss_iterator_state = VNACAL_NI_END_EQUATIONS;

    return 0;
}

/*
 * vs_start_frequency: start a new frequency
 *   @vnssp:  solve state structure
 *   @findex: frequency index
 */
static int vs_start_frequency(vnacal_new_solve_state_t *vnssp, int findex)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const double frequency = vnp->vn_frequency_vector[findex];
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);

    /*
     * Set the new frequency index.
     */
    vnssp->vnss_findex = findex;

    /*
     * If the error term type uses leakage terms outside of the linear
     * system, find the sum, sum of squared magnitude and count for
     * each term.
     */
    if (vnssp->vnss_leakage_matrix != NULL) {
	vnacal_new_measurement_t *vnmp;

	for (int i = 0; i < m_rows * m_columns; ++i) {
	    vnacal_new_leakage_term_t *vnltp = vnssp->vnss_leakage_matrix[i];

	    if (vnltp != NULL) {
		vnltp->vnlt_sum   = 0.0;
		vnltp->vnlt_sumsq = 0.0;
		vnltp->vnlt_count = 0;
	    }
	}
	for (vnmp = vnp->vn_measurement_list; vnmp != NULL;
		vnmp = vnmp->vnm_next) {
	    for (int row = 0; row < m_rows; ++row) {
		for (int column = 0; column < m_columns; ++column) {
		    const int m_cell = row * m_columns + column;
		    const int s_cell = row * s_columns + column;
		    vnacal_new_leakage_term_t *vnltp;
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
		    vnltp = vnssp->vnss_leakage_matrix[m_cell];
		    vnltp->vnlt_sum   += m;
		    vnltp->vnlt_sumsq += creal(m * conj(m));
		    ++vnltp->vnlt_count;
		}
	    }
	}
    }

    /*
     * Initialize the unknown parameter values.
     */
    for (vnacal_new_parameter_t *vnprp = vnp->vn_unknown_parameter_list;
	    vnprp != NULL; vnprp = vnprp->vnpr_next_unknown) {
	vnssp->vnss_p_vector[vnprp->vnpr_unknown_index][findex] =
	    _vnacal_get_parameter_value_i(vnprp->vnpr_parameter, frequency);
    }

    /*
     * For each measured standard...
     */
    for (vnacal_new_measurement_t *vnmp = vnp->vn_measurement_list;
	    vnmp != NULL; vnmp = vnmp->vnm_next) {
	vnacal_new_ms_matrices_t *vnmmp;

	/*
	 * Fill vnmm_m_matrix, subtracting out off-diagonal leakage
	 * terms if present.
	 */
	vnmmp = &vnssp->vnss_ms_matrices[vnmp->vnm_index];
	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		const int m_cell = m_row * m_columns + m_column;
		vnacal_new_leakage_term_t *vnltp;
		double complex value;

		if (vnmp->vnm_m_matrix[m_cell] == NULL) {
#ifdef NAN
		    value = NAN;
#else
		    value = 0.0;
#endif
		} else {
		    value = vnmp->vnm_m_matrix[m_cell][findex];
		    if (vnssp->vnss_leakage_matrix != NULL) {
			vnltp = vnssp->vnss_leakage_matrix[m_cell];
			if (vnltp != NULL && vnltp->vnlt_count > 0) {
			    value -= vnltp->vnlt_sum / vnltp->vnlt_count;
			}
		    }
		}
		vnmmp->vnmm_m_matrix[m_cell] = value;
	    }
	}

	/*
	 * Fill vnmm_s_matrix, interpolating between frequency points
	 * as necessary.  The _vnacal_get_parameter_value_i function
	 * returns initial guesses for unknown parameters.
	 */
	for (int s_row = 0; s_row < s_rows; ++s_row) {
	    for (int s_column = 0; s_column < s_columns; ++s_column) {
		const int s_cell = s_row * s_columns + s_column;
		vnacal_new_parameter_t *vnprp;
		double complex value;

		if ((vnprp = vnmp->vnm_s_matrix[s_cell]) == NULL) {
#ifdef NAN
		    value = NAN;
#else
		    value = 0.0;
#endif

		} else if (vnprp->vnpr_unknown) {
		    int uindex = vnprp->vnpr_unknown_index;

		    value = vnssp->vnss_p_vector[uindex][findex];

		} else {
		    value = _vnacal_get_parameter_value_i(vnprp->vnpr_parameter,
			    frequency);
		}
		vnmmp->vnmm_s_matrix[s_cell] = value;
	    }
	}
    }
    vnssp->vnss_iterator_state = VNACAL_NI_INIT;
    return 0;
}

/*
 * vs_start_system: prepare equation iterator for new system
 *   @vnssp: solve state structure
 *   @vnsp:  vnacal_new_system_t structure
 */
static inline void vs_start_system(vnacal_new_solve_state_t *vnssp, int sindex)
{
    vnssp->vnss_iterator_state = VNACAL_NI_SYSTEM;
    vnssp->vnss_sindex = sindex;
    vnssp->vnss_vnep   = NULL;
    vnssp->vnss_vncp   = NULL;
}

/*
 * vs_next_equation: move to the next equation in the system
 *   @vnssp: solve state structure
 */
static bool vs_next_equation(vnacal_new_solve_state_t *vnssp)
{
    /*
     * If there are equations remaining, move to the first / next.
     */
    switch (vnssp->vnss_iterator_state) {
    case VNACAL_NI_INIT:
	abort();	/* must call _vnacal_new_ss_start_system first */

    /*
     * If we're starting a new system, set vnss_vnep to the first
     * equation.
     */
    case VNACAL_NI_SYSTEM:
	{
	    vnacal_new_t *vnp = vnssp->vnss_vnp;
	    vnacal_new_system_t *vnsp;

	    vnsp = &vnp->vn_system_vector[vnssp->vnss_sindex];
	    vnssp->vnss_vnep = vnsp->vns_equation_list;
	}
	break;

    /*
     * If we're already started, advance to the next equation.  It's
     * permitted to advance to the next equation even if we haven't
     * started or completed iteration through the coefficients.
     */
    case VNACAL_NI_EQUATION:
    case VNACAL_NI_COEFFICIENT:
    case VNACAL_NI_END_COEFFICIENTS:
	vnssp->vnss_vnep = vnssp->vnss_vnep->vne_next;
	vnssp->vnss_vncp = NULL;
	break;

    /*
     * If we're at the end of the equations, keep returning false.
     */
    case VNACAL_NI_END_EQUATIONS:
	return false;
    }

    /*
     * If there are no remaining equations, set the state to end and
     * return false.  Otherwise, set the state to in equation, ready to
     * begin iterating over the coefficients.
     */
    if (vnssp->vnss_vnep == NULL) {
	vnssp->vnss_iterator_state = VNACAL_NI_END_EQUATIONS;
	return false;
    }
    vnssp->vnss_iterator_state = VNACAL_NI_EQUATION;

    return true;
}

/*
 * vs_next_coefficient: move to the next coefficient
 *   @vnssp: solve state structure
 */
static bool vs_next_coefficient(vnacal_new_solve_state_t *vnssp)
{
    /*
     * If we're starting a new equation, set vnss_vncp to the first
     * term.  If we're already in the term list, advance to the next
     * term.  If we're at the end of the equation, keep returning false;
     */
    switch (vnssp->vnss_iterator_state) {
    case VNACAL_NI_INIT:
    case VNACAL_NI_SYSTEM:
	abort();		/* must call vs_next_equation first */
	break;

    case VNACAL_NI_EQUATION:
	vnssp->vnss_vncp = vnssp->vnss_vnep->vne_coefficient_list;
	vnssp->vnss_iterator_state = VNACAL_NI_COEFFICIENT;
	break;

    case VNACAL_NI_COEFFICIENT:
	vnssp->vnss_vncp = vnssp->vnss_vncp->vnc_next;
	break;

    case VNACAL_NI_END_COEFFICIENTS:
    case VNACAL_NI_END_EQUATIONS:
	return false;
    }
    if (vnssp->vnss_vncp == NULL) {
	vnssp->vnss_iterator_state = VNACAL_NI_END_COEFFICIENTS;
	return false;
    }
    return true;
}

/*
 * vs_get_coefficient: return the current coefficient index or -1 for RHS
 *   @vnssp: solve state structure
 */
static inline int vs_get_coefficient(const vnacal_new_solve_state_t *vnssp)
{
    assert(vnssp->vnss_iterator_state == VNACAL_NI_COEFFICIENT);
    return vnssp->vnss_vncp->vnc_coefficient;
}

/*
 * vs_get_negative: test if the currnet coefficient has a minus sign
 *   @vnssp: solve state structure
 */
static inline bool vs_get_negative(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vncp->vnc_negative;
}

/*
 * vs_have_m: test if the current coefficient has an m factor
 *   @vnssp: solve state structure
 */
static inline bool vs_have_m(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vncp->vnc_m_cell >= 0;
}

/*
 * vs_get_m: return the m value for the current coefficient
 *   @vnssp: solve state structure
 */
static inline double complex vs_get_m(const vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_coefficient_t *vncp = vnssp->vnss_vncp;
    int m_cell;
    vnacal_new_measurement_t *vnmp;
    vnacal_new_ms_matrices_t *vnmmp;

    assert(vnssp->vnss_iterator_state == VNACAL_NI_COEFFICIENT);
    m_cell = vncp->vnc_m_cell;
    assert(m_cell >= 0);
    vnmp = vnssp->vnss_vnep->vne_vnmp;
    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
    vnmmp = &vnssp->vnss_ms_matrices[vnmp->vnm_index];

    return vnmmp->vnmm_m_matrix[m_cell];
}

/*
 * vs_get_m_cell: get the index in the m matrix for the current coefficient
 *   @vnssp: solve state structure
 */
static inline int vs_get_m_cell(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vncp->vnc_m_cell;
}

/*
 * vs_have_s: test if the current coefficient has an s factor
 *   @vnssp: solve state structure
 */
static inline bool vs_have_s(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vncp->vnc_s_cell >= 0;
}

/*
 * vs_get_s: return the s value for the current coefficient
 *   @vnssp: solve state structure
 */
static inline double complex vs_get_s(const vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_coefficient_t *vncp = vnssp->vnss_vncp;
    int s_cell;
    vnacal_new_measurement_t *vnmp;
    vnacal_new_ms_matrices_t *vnmmp;

    assert(vnssp->vnss_iterator_state == VNACAL_NI_COEFFICIENT);
    s_cell = vncp->vnc_s_cell;
    assert(s_cell >= 0);
    vnmp = vnssp->vnss_vnep->vne_vnmp;
    assert(vnmp->vnm_s_matrix[s_cell] != NULL);
    vnmmp = &vnssp->vnss_ms_matrices[vnmp->vnm_index];

    return vnmmp->vnmm_s_matrix[s_cell];
}

/*
 * vs_get_s_cell: get the index in the s matrix for the current coefficient
 *   @vnssp: solve state structure
 */
static inline int vs_get_s_cell(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vncp->vnc_s_cell;
}

/*
 * vs_update_s_matrices: update unknown parameters in s matrices
 *   @vnssp: solve state structure
 */
static void vs_update_s_matrices(vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int s_rows = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    int findex = vnssp->vnss_findex;

    /*
     * For each measured standard...
     */
    for (vnacal_new_measurement_t *vnmp = vnp->vn_measurement_list;
	    vnmp != NULL; vnmp = vnmp->vnm_next) {
	vnacal_new_ms_matrices_t *vnmmp =
	    &vnssp->vnss_ms_matrices[vnmp->vnm_index];
	double complex *s_matrix = vnmmp->vnmm_s_matrix;

	/*
	 * Patch s_matrix with the current value of the unknown parameters.
	 */
	for (int s_row = 0; s_row < s_rows; ++s_row) {
	    for (int s_column = 0; s_column < s_columns; ++s_column) {
		const int s_cell = s_row * s_columns + s_column;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];
		int uindex = vnprp->vnpr_unknown_index;

		if (vnprp != NULL && vnprp->vnpr_unknown) {
		    s_matrix[s_cell] = vnssp->vnss_p_vector[uindex][findex];
		}
	    }
	}
    }
}

/*
 * vs_free: free resources held by the solve state structure
 *   @vnssp: solve state structure
 */
static void vs_free(vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);

    if (vnssp->vnss_p_vector != NULL) {
	for (int i = vnp->vn_unknown_parameters - 1; i >= 0; --i) {
	    free((void *)vnssp->vnss_p_vector[i]);
	}
	free((void *)vnssp->vnss_p_vector);
	vnssp->vnss_p_vector = NULL;
    }
    if (vnssp->vnss_leakage_matrix != NULL) {
	for (int cell = m_rows * m_columns - 1; cell >= 0; --cell) {
	    free((void *)vnssp->vnss_leakage_matrix[cell]);
	}
	free((void *)vnssp->vnss_leakage_matrix);
	vnssp->vnss_leakage_matrix = NULL;
    }
    if (vnssp->vnss_ms_matrices != NULL) {
	for (int i = vnp->vn_measurement_count - 1; i >= 0; --i) {
	    vnacal_new_ms_matrices_t *vnmmp = &vnssp->vnss_ms_matrices[i];

	    free((void *)vnmmp->vnmm_s_matrix);
	    free((void *)vnmmp->vnmm_m_matrix);
	}
	free((void *)vnssp->vnss_ms_matrices);
	vnssp->vnss_ms_matrices = NULL;
    }
}

#ifdef DEBUG
/*
 * print_rmatrix: print an m by n serialized real matrix in octave form
 *   @name: name of matrix
 *   @a: pointer to first element of matrix
 *   @m: number of rows
 *   @n: number of columns
 */
static void print_rmatrix(const char *name, double *a, int m, int n)
{
    (void)printf("%s = [\n", name);
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    (void)printf(" %9.5f", a[i * n + j]);
	}
	(void)printf("\n");
    }
    (void)printf("]\n");
}

/*
 * print_cmatrix: print an m by n serialized complex matrix in octave form
 *   @name: name of matrix
 *   @a: pointer to first element of matrix
 *   @m: number of rows
 *   @n: number of columns
 */
static void print_cmatrix(const char *name, double complex *a, int m, int n)
{
    (void)printf("%s = [\n", name);
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    double complex v = a[i * n + j];

	    (void)printf(" %9.5f%+9.5fj", creal(v), cimag(v));
	}
	(void)printf("\n");
    }
    (void)printf("]\n");
}
#endif

/*
 * _vnacal_new_solve_simple: solve error terms where all s-parameters are known
 *   @vnssp: solve state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 */
static int _vnacal_new_solve_simple(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int unknowns = vlp->vl_t_terms - 1;

    /*
     * For each system of equations...
     */
    assert(x_length == vnp->vn_systems * unknowns);
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
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
	vs_start_system(vnssp, sindex);
	while (vs_next_equation(vnssp)) {
	    while (vs_next_coefficient(vnssp)) {
		int coefficient = vs_get_coefficient(vnssp);
		double complex value = vs_get_negative(vnssp) ? -1.0 : 1.0;

		if (vs_have_m(vnssp)) {
		    value *= vs_get_m(vnssp);
		}
		if (vs_have_s(vnssp)) {
		    value *= vs_get_s(vnssp);
		}
		if (coefficient == -1) {
		    b_vector[eq_count] = value;
		} else {
		    a_matrix[eq_count][coefficient] = value;
		}
	    }
	    ++eq_count;
	}
#ifdef DEBUG
	print_cmatrix("a", &a_matrix[0][0], eq_count, x_length);
	print_cmatrix("b", b_vector, eq_count, 1);
#endif /* DEBUG */

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
 *   @vnssp: solve state structure
 *   @x_vector: current error parameter solution
 *   @w_vector: new weight vector
 *
 * TODO: we're not calculating weights according to the more elegant
 * solution given in the Van hamme paper.  We need to do some rework of
 * the way we represent equations in order to do it right. This is kind
 * of close, though.
 */
static void calc_weights(vnacal_new_solve_state_t *vnssp,
	const double complex *x_vector, double *w_vector)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const int findex = vnssp->vnss_findex;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_cells = VL_M_COLUMNS(vlp) * VL_M_ROWS(vlp);
    vnacal_new_m_error_t *m_error_vector = vnp->vn_m_error_vector;
    const double noise = m_error_vector[findex].vnme_noise;
    const double tracking = m_error_vector[findex].vnme_tracking;
    double complex m_weight_vector[m_cells];
    int equation = 0;

    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	int offset = sindex * (vlp->vl_t_terms - 1);

	vs_start_system(vnssp, sindex);
	while (vs_next_equation(vnssp)) {
	    vnacal_new_measurement_t *vnmp;
	    vnacal_new_ms_matrices_t *vnmmp;
	    double u = 0.0;

	    vnmp = vnssp->vnss_vnep->vne_vnmp;
	    vnmmp = &vnssp->vnss_ms_matrices[vnmp->vnm_index];
	    for (int m_cell = 0; m_cell < m_cells; ++m_cell) {
		m_weight_vector[m_cell] = 0.0;
	    }
	    while (vs_next_coefficient(vnssp)) {
		int coefficient = vs_get_coefficient(vnssp);
		int m_cell = vs_get_m_cell(vnssp);

		if (m_cell >= 0) {
		    double complex v = 1.0;

		    if (vs_get_negative(vnssp)) {
			v *= -1.0;
		    }
		    if (vs_have_s(vnssp)) {
			v *= vs_get_s(vnssp);
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
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    m = vnmmp->vnmm_m_matrix[m_cell];
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
 * _vnacal_new_solve_auto: solve for both error terms and unknown s-parameters
 *   @vnssp: solve state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 *
 * This implementation is based on the algorithm described in H. Van Hamme
 * and M. Vanden Bossche, "Flexible vector network analyzer calibration
 * with accuracy bounds using an 8-term or a 16-term error correction
 * model," in IEEE Transactions on Microwave Theory and Techniques,
 * vol. 42, no. 6, pp. 976-987, June 1994, doi: 10.1109/22.293566.  There
 * are a few differences, however.  For example, instead of calculating
 * the error bounds on the error parameters, we simply test if the data
 * are consistent with the given linear model and error model.  Also,
 * we do the equation weighting wrong instead of using the "V" matrices.
 */
static int _vnacal_new_solve_auto(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length)
{
    /* pointer to vnacal_new_t structore */
    vnacal_new_t *vnp = vnssp->vnss_vnp;

    /* current frequency index */
    const int findex = vnssp->vnss_findex;

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
     * Test that we have at least as many equations as unknowns.
     */
    equations = vnp->vn_equations;
    assert(x_length == vnp->vn_systems * (vlp->vl_t_terms - 1));
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
     * Iterate using Gauss-Newton to find the unknown parameters, vnss_p_vector.
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

	/* sum of squared magnitudes of the elements of d_vector */
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
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    /*
	     * The vs_start_system, vs_next_equation and vs_next_coefficient
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
	    vs_start_system(vnssp, sindex);
	    while (vs_next_equation(vnssp)) {
		while (vs_next_coefficient(vnssp)) {
		    int coefficient = vs_get_coefficient(vnssp);
		    double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;

		    if (vs_have_m(vnssp)) {
			v *= vs_get_m(vnssp);
		    }
		    if (vs_have_s(vnssp)) {
			v *= vs_get_s(vnssp);
		    }
		    if (w_vector != NULL) {
			v *= w_vector[equation];
		    }
		    if (coefficient == -1) {
			b_vector[equation] = v;
		    } else {
			a_matrix[equation][offset + coefficient] = v;
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
	 * of vnss_p_vector, and restart the loop from the top.
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
	    calc_weights(vnssp, x_vector, w_vector);
#ifdef DEBUG
	    print_rmatrix("w", w_vector, equations, 1);
#endif /* DEBUG */
	    --iteration;
	    continue;
	}

#ifdef DEBUG
    (void)printf("p = [\n");
    for (int i = 0; i < p_length; ++i) {
	(void)printf("  %f%+fj\n",
		creal(vnssp->vnss_p_vector[i][findex]),
		cimag(vnssp->vnss_p_vector[i][findex]));
    }
    (void)printf("]\n\n");
#endif /* DEBUG */

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
	    int offset = sindex * (vlp->vl_t_terms - 1);

	    vs_start_system(vnssp, sindex);
	    while (vs_next_equation(vnssp)) {
		vnacal_new_equation_t *vnep = vnssp->vnss_vnep;
		vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;

		while (vs_next_coefficient(vnssp)) {
		    int coefficient = vs_get_coefficient(vnssp);
		    int s_cell = vs_get_s_cell(vnssp);
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
			double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;
			int unknown = vnprp->vnpr_unknown_index;

			if (vs_have_m(vnssp)) {
			    v *= vs_get_m(vnssp);
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
#ifdef DEBUG
	print_cmatrix("a", &a_matrix[0][0], equations, x_length);
	print_cmatrix("b", b_vector, equations, 1);
#endif /* DEBUG */

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
		    vnssp->vnss_p_vector[pindex1][findex];/*contr. to residual*/
		if (vnprp2->vnpr_unknown) {
		    int pindex2 = vnprp2->vnpr_unknown_index;
		    j_matrix[j_row][pindex2] = -coefficient; /* partial drvtv */
		    k_vector[j_row] -= coefficient *
			vnssp->vnss_p_vector[pindex2][findex];   /* resid */
		} else { /* known parameter value */
		    k_vector[j_row] -= coefficient *	/* contr. to resid */
			_vnacal_get_parameter_value_i(vnprp2->vnpr_parameter,
				frequency);
		}
		++j_row;
	    }
	    assert(j_row == j_rows);
	}
#ifdef DEBUG
	print_cmatrix("j", &j_matrix[0][0], j_rows, p_length);
	print_cmatrix("k", k_vector, j_rows, 1);
#endif /* DEBUG */

	/*
	 * Solve the j_matrix, k_vector system to create d_vector, the
	 * Gauss-Newton correction to vnss_p_vector.
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
	print_cmatrix("d", d_vector, p_length, 1);
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
	 * If the error is within the target tolerance, stop.
	 */
	if (sum_d_squared / (double)p_length <= vnp->vn_p_tolerance *
						vnp->vn_p_tolerance) {
#ifdef DEBUG
	    (void)printf("# stop: converged\n");
#endif /* DEBUG */
	    break;
	}

	/*
	 * If we have the best solution so far (or the first solution),
	 * add the new correction to vnss_p_vector and remember the solution.
	 */
	if (sum_d_squared < best_sum_d_squared) {
	    double sum_p_squared = 0.0;

#ifdef DEBUG
	    (void)printf("# best\n");
#endif /* DEBUG */
	    /*
	     * Limit the magnitude of d_vector to keep it smaller than
	     * the magnitude of vnss_p_vector (or smaller than one if
	     * vnss_p_vector is less than one).  This improves stability
	     * at the cost of slowing convergence.  It makes it less
	     * likely that we jump into an adjacent basin of attraction.
	     * We use one over the golden ratio as the maximum norm of
	     * d_vector relative to vnss_p_vector (or one).
	     */
	    for (int i = 0; i < p_length; ++i) {
		double complex p = vnssp->vnss_p_vector[i][findex];

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
		best_p_vector[i] = vnssp->vnss_p_vector[i][findex];
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
		calc_weights(vnssp, x_vector, w_vector);
#ifdef DEBUG
		for (int i = 0; i < equations; ++i) {
		    (void)printf("# w[%2d] = %f\n", i, w_vector[i]);
		}
#endif /* DEBUG */
	    }

	    /*
	     * Apply d_vector to vnss_p_vector.
	     */
	    for (int i = 0; i < p_length; ++i) {
		vnssp->vnss_p_vector[i][findex] += d_vector[i];
	    }
#ifdef DEBUG
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("# p[%2d] = %13.6e %+13.6ej\n", i,
			creal(vnssp->vnss_p_vector[i][findex]),
			cimag(vnssp->vnss_p_vector[i][findex]));
	    }
#endif /* DEBUG */
	    backtrack_count = 0;

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
		vnssp->vnss_p_vector[i][findex] =
		    best_p_vector[i] + best_d_vector[i];
	    }
#ifdef DEBUG
	    for (int i = 0; i < p_length; ++i) {
		(void)printf("# p[%2d] = %13.6e %+13.6ej\n", i,
			creal(vnssp->vnss_p_vector[i][findex]),
			cimag(vnssp->vnss_p_vector[i][findex]));
	    }
#endif /* DEBUG */
	    if (w_vector != NULL) {
		(void)memcpy((void *)w_vector, (void *)best_w_vector,
			equations * sizeof(double));
		calc_weights(vnssp, x_vector, w_vector);
#ifdef DEBUG
		print_rmatrix("w", w_vector, equations, 1);
#endif /* DEBUG */
	    }
	}
	vs_update_s_matrices(vnssp);

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
	vnssp->vnss_p_vector[i][findex] = best_p_vector[i];
    }
    vs_update_s_matrices(vnssp);

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
 *   @vnssp: solve state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 */
static double calc_rms_error(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const int findex = vnssp->vnss_findex;
    const double frequency = vnp->vn_frequency_vector[findex];
    vnacal_new_leakage_term_t **leakage_matrix = vnssp->vnss_leakage_matrix;
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
	int offset = sindex * (vlp->vl_t_terms - 1);

	vs_start_system(vnssp, sindex);
	while (vs_next_equation(vnssp)) {
	    double complex residual = 0.0;
	    double u = 0.0;

	    for (int i = 0; i < w_terms; ++i) {
		w_term_vector[i] = 0.0;
	    }
	    while (vs_next_coefficient(vnssp)) {
		int coefficient = vs_get_coefficient(vnssp);
		double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;
		int m_cell = vs_get_m_cell(vnssp);

		if (coefficient >= 0) {
		    assert(offset + coefficient < x_length);
		    v *= x_vector[offset + coefficient];
		} else {
		    v *= -1.0;
		}
		if (vs_have_s(vnssp)) {
		    v *= vs_get_s(vnssp);
		}
		if (vs_have_m(vnssp)) {
		    double complex m = vs_get_m(vnssp);
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

		v = vnssp->vnss_p_vector[vnprp1->vnpr_unknown_index][findex];
		if (vnprp2->vnpr_unknown) {
		    v -= vnssp->vnss_p_vector[vnprp2->vnpr_unknown_index][findex];
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
    if (leakage_matrix != NULL) {
	for (int row = 0; row < m_rows; ++row) {
	    for (int column = 0; column < m_columns; ++column) {
		if (row != column) {
		    const int m_cell = row * m_columns + column;
		    const vnacal_new_leakage_term_t *vnltp;
		    double value;

		    vnltp = leakage_matrix[m_cell];
		    if (vnltp->vnlt_count > 1) {
			double complex sum_x = vnltp->vnlt_sum;
			const int n = vnltp->vnlt_count;
			double n_mean_squared;

			n_mean_squared = creal(sum_x * conj(sum_x)) / n;
			value = (vnltp->vnlt_sumsq - n_mean_squared);
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
	const vnacal_layout_t *vlp_in, const vnacal_layout_t *vlp_out)
{
    const int m_rows     = VL_M_ROWS(vlp_in);
    const int m_columns  = VL_M_COLUMNS(vlp_in);
    const int e12_terms  = VL_ERROR_TERMS(vlp_out);
    const double complex *el_in = &e[VL_EL_OFFSET(vlp_in)];
    double complex e_out[e12_terms];
    int el_map[m_rows * m_columns];

    /*
     * The el_in vector contains only the off-diagonal terms.  Construct
     * a map from m_cell to index of el_in so can easily find them.
     */
    {
	int index = 0;

	for (int row = 0; row < m_rows; ++row) {
	    for (int column = 0; column < m_columns; ++column) {
		const int m_cell = row * m_columns + column;

		el_map[m_cell] = (row == column) ? -1 : index++;
	    }
	}
    }

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
		const int m_cell = m_row * m_columns + m_column;

		el[m_row] = el_in[el_map[m_cell]];
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
    const int m_rows    = VL_M_ROWS(vlp_in);
    const int m_columns = VL_M_COLUMNS(vlp_in);
    const int error_terms_in = VL_ERROR_TERMS(vlp_in);
    int error_terms_out = error_terms_in;
    vnacal_new_solve_state_t vnss;
    vnacal_calibration_t *calp = NULL;
    int rc = -1;

    /*
     * Make sure the frequency vector was given.
     */
    if (!vnp->vn_frequencies_valid) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_solve: "
		"calibration frequency " "vector must be given");
	return -1;
    }

    /*
     * Init the state structure.
     */
    if (vs_init(&vnss, vnp) == -1) {
	return -1;
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
	const int x_length = vnp->vn_systems * (vlp_in->vl_t_terms - 1);
	double complex x_vector[x_length];
	double complex e_vector[error_terms_out];
	int eterm_index = 0;

	/*
	 * Prepare the state structure for a new frequency.
	 */
	vs_start_frequency(&vnss, findex);

	/*
	 * Solve the system.  If there are no unknown parameters,
	 * the the system is linear and we can use a simple analytic
	 * method to solve it.  Otherwise, the system is non-linear
	 * and we use an iterative gauss-newton.
	 */
	if (unknown_parameters == 0) {
	    if (_vnacal_new_solve_simple(&vnss, x_vector, x_length) == -1) {
		goto out;
	    }
	} else {
	    if (_vnacal_new_solve_auto(&vnss, x_vector, x_length) == -1) {
		goto out;
	    }
	}

	/*
	 * If the mesurement error was given, calculate the RMS
	 * error of the solution and fail if it's too high.
	 */
	if (vnp->vn_m_error_vector != NULL) {
	    double error;

	    error = calc_rms_error(&vnss, x_vector, x_length);
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
	 * If there are leakage terms outside of the linear system, add them.
	 */
	if (VL_HAS_OUTSIDE_LEAKAGE_TERMS(vlp_in)) {
	    for (int row = 0; row < m_rows; ++row) {
		for (int column = 0; column < m_columns; ++column) {
		    if (row != column) {
			const int cell = row * m_columns + column;
			vnacal_new_leakage_term_t *vnltp;
			double complex term;

			vnltp = vnss.vnss_leakage_matrix[cell];
			if (vnltp->vnlt_count != 0) {
			    term = vnltp->vnlt_sum / vnltp->vnlt_count;
			} else {
			    term = 0.0;
			}
			e_vector[eterm_index++] = term;
		    }
		}
	    }
	}

	/*
	 * If _VNACAL_E12_UE14, convert to VNACAL_E12.
	 */
	if (type_in == _VNACAL_E12_UE14) {
	    rc = convert_ue14_to_e12(e_vector, vlp_in, vlp_out);
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
	assert(vnss.vnss_p_vector[index] != NULL);
	vpmrp->vpmr_gamma_vector = vnss.vnss_p_vector[index];
	vnss.vnss_p_vector[index] = NULL;
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
    vs_free(&vnss);
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
