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
    double complex *prev_x_vector = NULL;
    int rv = -1;

    /*
     * If a measurement error vector was given, calculates weights
     * for each measurement and allocate the prev_v_vector.
     */
    if (vnp->vn_m_error_vector != NULL) {
	if ((w_vector = vs_calc_weights(vnssp)) == NULL) {
	    goto out;
	}
	prev_x_vector = malloc(x_length * sizeof(double complex));
	if (prev_x_vector == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	    goto out;
	}
	_vnacal_new_solve_init_x_vector(vnssp, prev_x_vector, x_length);
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
	    double sum_dx_squared;
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
	    /*
	     * If measurement errors were not given or if this particular
	     * system is not over-determined, then we don't need to iterate.
	     */
	    if (!vs_have_v(vnssp)) {
		break;
	    }

	    /*
	     * Here, x_vector depends on the V matrices and the V matrices
	     * depend on x_vector.  Iterate until they converge.
	     */
	    if (vs_update_v_matrices("vnacal_new_solve", vnssp, sindex,
			&x_vector[offset], unknowns) == -1) {
		return -1;
	    }
	    sum_dx_squared = 0.0;
	    for (int i = 0; i < x_length; ++i) {
		double complex d = x_vector[i] - prev_x_vector[i];

		sum_dx_squared += _vnacommon_cabs2(d);
	    }
	    if (sum_dx_squared / (double)x_length <= vnp->vn_et_tolerance *
						     vnp->vn_et_tolerance) {
		break;
	    }
	    if (++iteration >= vnp->vn_iteration_limit) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"measurement error model failed to converge at %e Hz",
			frequency);
		goto out;
	    }
	    (void)memcpy((void *)prev_x_vector, (void *)x_vector,
		    x_length * sizeof(double complex));
	}
    }
    rv = 0;

out:
    free((void *)prev_x_vector);
    free((void *)w_vector);
    return rv;
}
