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

/*
 * vnacal_new_leakage_term_t
 */
typedef struct vnacal_new_leakage_term {
    /* sum of the samples */
    double complex ncl_sx;

    /* sum of the squared magnitudes of the samples */
    double ncl_sx2;

    /* count of samples */
    int ncl_count;

} vnacal_new_leakage_term_t;

/*
 * TODO: in the next phase, we need to build the Jacobian matrix for
 * unknown parameters.  Note that the squared residual of the leakage
 * terms is:
 *       sum(x[i] * conj(x[i])) - 1/N sum(x[i]) * conj(sum(x[i])), or
 *       ncl_sx2 - ncl_sx * conj(ncl_sx) / N
 */

/*
 * get_measurement: get a measurement value adjusted for leakage
 *   @vnmp: measurement of a calibration standard
 *   @leakage_term_vector: vector of leakage term structures
 *   @leakage_term_map: map from m_cell to leakage term
 *   @findex: frequency index
 *   @m_cell: linear offset into measurement matrix
 */
static double complex get_measurement(const vnacal_new_measurement_t *vnmp,
    const vnacal_new_leakage_term_t *leakage_term_vector,
    const int *leakage_term_map, int findex, int m_cell)
{
    const vnacal_new_t *vnp = vnmp->vnm_ncp;
    const vnacal_layout_t *const vlp = &vnp->vn_layout;
    const int m_row = m_cell / VL_M_COLUMNS(vlp);
    const int m_column = m_cell % VL_M_COLUMNS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);

    if (m_row != m_column && leakage_term_vector != NULL &&
	    m_row < s_rows && m_column < s_columns) {
	const int m_cell = m_row * m_columns + m_column;
	const int term = leakage_term_map[m_cell];
	const vnacal_new_leakage_term_t *nclp = &leakage_term_vector[term];

	if (nclp->ncl_count > 0) {
	    return vnmp->vnm_m_matrix[m_cell][findex] -
		nclp->ncl_sx / nclp->ncl_count;
	}
    }
    return vnmp->vnm_m_matrix[m_cell][findex];
}

/*
 * convert_ue14_to_e12: convert UE14 error terms to E12 error terms
 *   @e: input and output vector
 *   @vlp_in: input layout
 *   @vlp_out: output layout
 *   @leakage_term_map: map m_cell to external leakage term
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
	const int *leakage_term_map)
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
		el[m_row] = el_in[leakage_term_map[m_cell]];
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
    const vnacal_layout_t *const vlp_in = &vnp->vn_layout;
    const vnacal_layout_t *vlp_out = vlp_in;
    vnacal_layout_t vl_e12;
    const vnacal_type_t type_in = VL_TYPE(vlp_in);
    vnacal_type_t type_out = type_in;
    const int m_rows = VL_M_ROWS(vlp_in);
    const int m_columns = VL_M_COLUMNS(vlp_in);
    const int s_columns = VL_S_COLUMNS(vlp_in);
    const int frequencies = vnp->vn_frequencies;
    const int max_equations = vnp->vn_max_equations;
    const int coefficients = vlp_in->vl_t_terms - 1;
    const int leakage_terms = VL_EL_TERMS(vlp_in);
    int leakage_term_map[m_rows * m_columns];
    const int solved_error_terms = VL_ERROR_TERMS(vlp_in);
    int output_error_terms = solved_error_terms;
    double complex **unknown_parameter_vector = NULL;
    vnacal_calibration_t *calp = NULL;
    vnacal_new_leakage_term_t *leakage_term_vector = NULL;
    double complex *a_matrix = NULL;
    double complex *b_vector = NULL;
    double complex *x_vector = NULL;
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
     * If the type is _VNACAL_E12_UE14, set up a different output layout
     * for automatically converting to VNACAL_E12.
     */
    if (type_in == _VNACAL_E12_UE14) {
	_vnacal_layout(&vl_e12, VNACAL_E12, m_rows, m_columns);
	vlp_out = &vl_e12;
	type_out = VNACAL_E12;
	output_error_terms = VL_ERROR_TERMS(vlp_out);
    }

    /*
     * Create the vnacal_calibration_t structure.
     */
    if ((calp = _vnacal_calibration_alloc(vcp, type_out, m_rows, m_columns,
		    frequencies, output_error_terms)) == NULL) {
	goto out;
    }
    (void)memcpy((void *)calp->cal_frequency_vector,
	    (void *)vnp->vn_frequency_vector,
	    frequencies * sizeof(double));
    calp->cal_z0 = vnp->vn_z0;

    /*
     * If there are unknown parameters, allocate a vector, one entry per
     * frequency, of vectors of unknown parameter values.
     */
    if (vnp->vn_unknown_parameters != 0) {
	unknown_parameter_vector = vnp->vn_unknown_parameter_vector;
	if (unknown_parameter_vector == NULL) {
	    if ((unknown_parameter_vector = calloc(frequencies,
			    sizeof(double complex *))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"calloc: %s", strerror(errno));
		goto out;
	    }
	    vnp->vn_unknown_parameter_vector = unknown_parameter_vector;
	}
	for (int findex = 0; findex < frequencies; ++findex) {
	    free((void *)unknown_parameter_vector[findex]);
	    if ((unknown_parameter_vector[findex] =
			calloc(vnp->vn_unknown_parameters,
			    sizeof(double complex))) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"calloc: %s", strerror(errno));
		goto out;
	    }
	}
    }

    /*
     * Allocate the coefficient matrix and right-hand side vector.
     */
    if ((a_matrix = calloc(max_equations * coefficients,
		    sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto out;
    }
    if ((b_vector = calloc(max_equations,
		    sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto out;
    }
    if ((x_vector = calloc(coefficients,
		    sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto out;
    }

    /*
     * If there are leakage terms outside of the linear system, allocate
     * a vector of vnacal_new_leakage_term_t structures to accumulate
     * them and build a map from m_cell to leakage term.
     */
    if (leakage_terms != 0) {
	int term = 0;

	if ((leakage_term_vector = calloc(leakage_terms,
			sizeof(vnacal_new_leakage_term_t))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    goto out;
	}
	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		const int m_cell = m_row * m_columns + m_column;

		if (m_row == m_column) {
		    leakage_term_map[m_cell] = -1;
		} else {
		    leakage_term_map[m_cell] = term++;
		}
	    }
	}
    }

    /*
     * For each frequency, solve for the error parameters.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	const double frequency = vnp->vn_frequency_vector[findex];
	vnacal_new_parameter_t *vnprp;
	double complex error_term_vector[MAX(solved_error_terms,
					     output_error_terms)];
	int eterm_index = 0;

	/*
	 * Initialize the unknown parameter values, if any.
	 */
	for (vnprp = vnp->vn_unknown_parameter_list; vnprp != NULL;
		vnprp = vnprp->vnpr_next_unknown) {
	    const int i = vnprp->vnpr_unknown_index;

	    unknown_parameter_vector[findex][i] =
		_vnacal_get_parameter_value(vnprp->vnpr_parameter, frequency);
	}

	/*
	 * If TE10, UE10 or UE14, compute leakage terms outside of the
	 * linear system.
	 */
	if (leakage_terms != 0) {
	    vnacal_new_leakage_term_t *nclp;
	    vnacal_new_measurement_t *vnmp;

	    for (int term = 0; term < leakage_terms; ++term) {
		nclp = &leakage_term_vector[term];
		nclp->ncl_sx    = 0.0;
		nclp->ncl_sx2   = 0.0;
		nclp->ncl_count = 0;
	    }
	    for (vnmp = vnp->vn_measurement_list; vnmp != NULL;
		    vnmp = vnmp->vnm_next) {
		for (int row = 0; row < m_rows; ++row) {
		    for (int column = 0; column < m_columns; ++column) {
			const int m_cell = row * m_columns + column;
			const int s_cell = row * s_columns + column;
			const int term = leakage_term_map[m_cell];
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
			nclp = &leakage_term_vector[term];
			nclp->ncl_sx  += m;
			nclp->ncl_sx2 += creal(m * conj(m));
			++nclp->ncl_count;
		    }
		}
	    }
	    for (int term = 0; term < leakage_terms; ++term) {
		nclp = &leakage_term_vector[term];

		if (nclp->ncl_count == 0) {
		    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			    "leakage term system is singular");
		    goto out;
		}
	    }
	}

	/*
	 * For each system of equations...
	 */
	for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	    vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	    int unity_index = _vl_unity_offset(vlp_in, sindex);
	    int eq_count = 0;

	    /*
	     * Build a_matrix and b_vector.
	     */
	    for (int a_cell = 0; a_cell < max_equations * coefficients;
		    ++a_cell) {
		a_matrix[a_cell] = 0.0;
	    }
	    for (int b_row = 0; b_row < max_equations; ++b_row) {
		b_vector[b_row] = 0.0;
	    }
	    /*
	     * TODO: consider building per vnacal_new_measurement_t s and
	     * m matrices here to avoid repeated calls to get_measurement
	     * and _vnacal_get_parameter_value.
	     */
	    for (vnacal_new_equation_t *vnep = vnsp->vns_equation_list;
		    vnep != NULL; vnep = vnep->vne_next) {
		vnacal_new_measurement_t *vnmp = vnep->vne_vcsp;
		const int eq_row = vnep->vne_row;
		const int eq_column = vnep->vne_column;
		const int s_cell = eq_row * s_columns + eq_column;

		if (type_in != VNACAL_T16 && type_in != VNACAL_U16 &&
			eq_row != eq_column &&
			!vnmp->vnm_reachability_matrix[s_cell]) {
		    continue;
		}
		for (vnacal_new_term_t *vntp = vnep->vne_term_list;
			vntp != NULL; vntp = vntp->vnt_next) {
		    double complex coefficient = vntp->vnt_negative ?
			-1.0 : 1.0;

		    if (vntp->vnt_measurement >= 0) {
			coefficient *= get_measurement(vnmp,
				leakage_term_vector, leakage_term_map,
				findex, vntp->vnt_measurement);
		    }
		    if (vntp->vnt_sparameter >= 0) {
			vnacal_new_parameter_t *vnprp =
			    vnmp->vnm_s_matrix[vntp->vnt_sparameter];

			assert(vnprp != NULL);
			if (vnprp->vnpr_unknown) {
			    coefficient *=
				unknown_parameter_vector[findex][
				vnprp->vnpr_unknown_index];
			} else {
			    coefficient *=
				_vnacal_get_parameter_value(
					vnprp->vnpr_parameter, frequency);
			}
		    }
		    if (vntp->vnt_coefficient == -1) {
			b_vector[eq_count] = coefficient;
		    } else {
			a_matrix[eq_count * coefficients +
			    vntp->vnt_coefficient] = coefficient;
		    }
		}
		++eq_count;
	    }

	    // if measurement error given, scale rows (do it always?)

	    /*
	     * Solve the linear system.
	     */
	    if (eq_count < coefficients) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"insufficient number of equations to solve "
			"error terms");
		goto out;
	    }
	    if (eq_count == coefficients) {
		double complex determinant;

		determinant = _vnacommon_mldivide(x_vector, a_matrix, b_vector,
			coefficients, 1);
		if (determinant == 0.0 || !isfinite(cabs(determinant))) {
		    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			    "singular linear system");
		    goto out;
		}
	    } else {
		int rank;

		rank = _vnacommon_qrsolve(x_vector, a_matrix, b_vector,
			eq_count, coefficients, 1);
		if (rank < coefficients) {
		    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			    "singular linear system");
		    goto out;
		}
	    }

	    /*
	     * Copy from x_vector to error_term_vector, inserting the
	     * unity term.
	     */
	    for (int k = 0; k < unity_index; ++k) {
		error_term_vector[eterm_index++] = x_vector[k];
	    }
	    error_term_vector[eterm_index++] = 1.0;
	    for (int k = unity_index; k < coefficients; ++k) {
		error_term_vector[eterm_index++] = x_vector[k];
	    }
	}

	/*
	 * Add the leakage terms.
	 */
	for (int term = 0; term < leakage_terms; ++term) {
	    vnacal_new_leakage_term_t *nclp = &leakage_term_vector[term];
	    double complex el;

	    nclp = &leakage_term_vector[term];
	    if (nclp->ncl_count != 0) {
		el = nclp->ncl_sx / nclp->ncl_count;
	    } else {
		el = 0.0;
	    }
	    error_term_vector[eterm_index++] = el;
	}
	assert(eterm_index == solved_error_terms);

	/*
	 * If _VNACAL_E12_UE14, convert to VNACAL_E12.
	 */
	if (type_in == _VNACAL_E12_UE14) {
	    rc = convert_ue14_to_e12(error_term_vector, vlp_in, vlp_out,
		    leakage_term_map);
	    if (rc == -1) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"singular system");
		goto out;
	    }
	}

	/*
	 * TODO: handle unknown parameters here
	 */

	/*
	 * Copy the error terms the calibration structure.
	 */
	for (int term = 0; term < output_error_terms; ++term) {
	    calp->cal_error_term_vector[term][findex] = error_term_vector[term];
	}
    }
    _vnacal_calibration_free(vnp->vn_calibration);
    vnp->vn_calibration = calp;
    calp = NULL;	/* no longer ours to free */
    rc = 0;

out:
    free((void *)leakage_term_vector);
    free((void *)x_vector);
    free((void *)b_vector);
    free((void *)a_matrix);
    _vnacal_calibration_free(calp);
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
