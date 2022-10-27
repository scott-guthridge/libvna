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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"


/*
 * _vnacal_new_solve_init: initialize the solve state structure
 *   @vnssp: solve state structure
 *   @vnp:   associated vnacal_new_t structure
 */
int _vnacal_new_solve_init(vnacal_new_solve_state_t *vnssp, vnacal_new_t *vnp)
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
     * Init the equation iterator state.
     */
    vnssp->vnss_iterator_state = VNACAL_NI_END_EQUATIONS;

    return 0;
}

/*
 * _vnacal_new_solve_start_frequency: start a new frequency
 *   @vnssp:  solve state structure
 *   @findex: frequency index
 */
int _vnacal_new_solve_start_frequency(vnacal_new_solve_state_t *vnssp,
	int findex)
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
		    if (vnmp->vnm_connectivity_matrix[s_cell]) {
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
 * _vnacal_new_solve_next_equation: move to the next equation in the system
 *   @vnssp: solve state structure
 */
bool _vnacal_new_solve_next_equation(vnacal_new_solve_state_t *vnssp)
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
     * started or completed iteration through the terms.
     */
    case VNACAL_NI_EQUATION:
    case VNACAL_NI_TERM:
    case VNACAL_NI_END_TERMS:
	vnssp->vnss_vnep = vnssp->vnss_vnep->vne_next;
	vnssp->vnss_vntp = NULL;
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
     * begin iterating over the terms.
     */
    if (vnssp->vnss_vnep == NULL) {
	vnssp->vnss_iterator_state = VNACAL_NI_END_EQUATIONS;
	return false;
    }
    vnssp->vnss_iterator_state = VNACAL_NI_EQUATION;

    return true;
}

/*
 * _vnacal_new_solve_next_term: move to the next term
 *   @vnssp: solve state structure
 */
bool _vnacal_new_solve_next_term(vnacal_new_solve_state_t *vnssp)
{
    /*
     * If we're starting a new equation, set vnss_vntp to the first
     * term.  If we're already in the term list, advance to the next
     * term.  If we're at the end of the equation, keep returning false;
     */
    switch (vnssp->vnss_iterator_state) {
    case VNACAL_NI_INIT:
    case VNACAL_NI_SYSTEM:
	abort();		/* must call vs_next_equation first */
	break;

    case VNACAL_NI_EQUATION:
	vnssp->vnss_vntp = vnssp->vnss_vnep->vne_term_list_no_v;
	vnssp->vnss_iterator_state = VNACAL_NI_TERM;
	break;

    case VNACAL_NI_TERM:
	vnssp->vnss_vntp = vnssp->vnss_vntp->vnt_next_no_v;
	break;

    case VNACAL_NI_END_TERMS:
    case VNACAL_NI_END_EQUATIONS:
	return false;
    }
    if (vnssp->vnss_vntp == NULL) {
	vnssp->vnss_iterator_state = VNACAL_NI_END_TERMS;
	return false;
    }
    return true;
}

/*
 * _vnacal_new_solve_update_s_matrices: update unknown parameters in s matrices
 *   @vnssp: solve state structure
 */
void _vnacal_new_solve_update_s_matrices(vnacal_new_solve_state_t *vnssp)
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
 * _vnacal_new_solve_free: free resources held by the solve state structure
 *   @vnssp: solve state structure
 */
void _vnacal_new_solve_free(vnacal_new_solve_state_t *vnssp)
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
	    while (vs_next_term(vnssp)) {
		int xindex = vs_get_xindex(vnssp);
		double complex v = vs_get_negative(vnssp) ? -1.0 : 1.0;
		int m_cell = vs_get_m_cell(vnssp);

		if (xindex >= 0) {
		    assert(offset + xindex < x_length);
		    v *= x_vector[offset + xindex];
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
