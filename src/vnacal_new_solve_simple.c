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

/* #define DEBUG */

/*
 * init_prev_v_vector: initialize prev_v_vector to identity matrices
 *   @vnssp: pointer to state structure
 *   @prev_v_vector: vector filled in by this function
 */
static double complex *alloc_prev_v_vector(vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int v_rows    = VL_V_ROWS(vlp);
    const int v_columns = VL_V_COLUMNS(vlp);
    const int v_length = vnp->vn_measurement_count * v_rows * v_columns;
    double complex *prev_v_vector = NULL;
    int v_index = 0;

    if ((prev_v_vector = calloc(v_length, sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	return NULL;
    }
    for (int m_index = 0; m_index < vnp->vn_measurement_count; ++m_index) {
	for (int v_row = 0; v_row < v_rows; ++v_row) {
	    for (int v_column = 0; v_column < v_columns; ++v_column) {
		prev_v_vector[v_index++] = v_row == v_column ? 1.0 : 0.0;
	    }
	}
    }
    assert(v_index == v_length);
    return prev_v_vector;
}

/*
 * update_prev_v_vector: update prev_v_vector and return mean of squared diff.
 *   @vnssp: pointer to state structure
 *   @prev_v_vector: vector filled in by this function
 *   @system: index of current linear system
 */
static double update_prev_v_vector(vnacal_new_solve_state_t *vnssp,
	double complex *prev_v_vector, int sindex)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int v_rows    = VL_V_ROWS(vlp);
    const int v_columns = VL_V_COLUMNS(vlp);
    int v_index = 0;
    double sum = 0.0;

    for (int m_index = 0; m_index < vnp->vn_measurement_count; ++m_index) {
	vnacal_new_msv_matrices_t *vnmmp = &vnssp->vnss_msv_matrices[m_index];
	double complex *v_matrix = vnmmp->vnsm_v_matrices[sindex];

	for (int v_cell = 0; v_cell < v_rows * v_columns; ++v_cell) {
	    double complex d = v_matrix[v_cell] - prev_v_vector[v_index];

	    sum += _vnacommon_cabs2(d);
	    prev_v_vector[v_index++] = v_matrix[v_cell];
	}
    }
    if (v_index != 0) {
	sum /= (double)v_index;
    }
    return sum;
}

/*
 * _vnacal_new_solve_simple: solve when all s-parameters are known
 *   @vnssp: solve state structure
 *   @x_vector: vector of unknowns filled in by this function
 *   @x_length: length of x_vector
 */
int _vnacal_new_solve_simple(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length)
{
    vnacal_new_t *vnp = vnssp->vnss_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int findex = vnssp->vnss_findex;
    const int unknowns = vlp->vl_t_terms - 1;
    double frequency = vnp->vn_frequency_vector[findex];
    double *w_vector = NULL;
    double complex *prev_v_vector = NULL;
    int rv = -1;

    /*
     * If a measurement error vector was given, calculates weights
     * for each measurement and allocate the prev_v_vector.
     */
    if (vnp->vn_m_error_vector != NULL) {
	if ((w_vector = vs_calc_weights(vnssp)) == NULL) {
	    goto out;
	}
	if ((prev_v_vector = alloc_prev_v_vector(vnssp)) == NULL) {
	    goto out;
	}
    }

    /*
     * For each system of equations...
     */
    assert(x_length == vnp->vn_systems * unknowns);
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	const int offset = sindex * unknowns;
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	const int equations = vnsp->vns_equation_count;
	int iteration = 0;

	/*
	 * For each iteration on the V matrices (if in use)...
	 */
	for (;;) {
	    double complex a_matrix[equations][unknowns];
	    double complex b_vector[equations];
	    double v_mean_squares;
	    int eq_count = 0;

	    /*
	     * Build the coefficient matrix (a) and right-hand side vector (b).
	     */
	    for (int i = 0; i < equations; ++i) {
		for (int j = 0; j < unknowns; ++j) {
		    a_matrix[i][j] = 0.0;
		}
		b_vector[i] = 0.0;
	    }
	    vs_start_system(vnssp, sindex);
	    while (vs_next_equation(vnssp)) {
		while (vs_next_term(vnssp)) {
		    double complex value = vs_get_negative(vnssp) ? -1.0 : 1.0;
		    const int xindex = vs_get_xindex(vnssp);

		    if (vs_have_m(vnssp)) {
			value *= vs_get_m(vnssp);
		    }
		    if (vs_have_s(vnssp)) {
			value *= vs_get_s(vnssp);
		    }
		    if (vs_have_v(vnssp)) {
			value *= vs_get_v(vnssp);
		    }
		    if (w_vector != NULL) {
			value *= w_vector[eq_count];
		    }
		    if (xindex == -1) {
			b_vector[eq_count] += value;
		    } else {
			a_matrix[eq_count][xindex] += value;
		    }
		}
		++eq_count;
	    }
	    assert(eq_count == equations);

	    /*
	     * Solve for the unknowns using LU decomposition if a_matrix
	     * is square, or QR decomposition if the system is overdetermined.
	     */
	    if (equations < unknowns) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"insufficient number of standards to solve "
			"error terms");
		return -1;
	    }
	    if (equations == unknowns) {
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
	    if (vnp->vn_m_error_vector != NULL) {
		if (vs_update_v_matrices("vnacal_new_solve", vnssp, sindex,
			    &x_vector[offset], unknowns) == -1) {
		    return -1;
		}
	    }
	    if (!vs_have_v(vnssp)) {
		break;
	    }
	    v_mean_squares = update_prev_v_vector(vnssp, prev_v_vector, sindex);
	    if (v_mean_squares <= vnp->vn_v_tolerance * vnp->vn_v_tolerance) {
		break;
	    }
	    if (++iteration >= 50) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"measurement error model failed to converge at %e Hz",
			frequency);
		goto out;
	    }
	}
    }
    rv = 0;

out:
    free((void *)prev_v_vector);
    free((void *)w_vector);
    return rv;
}
