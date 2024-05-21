/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2023 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#define DEBUG

#ifdef DEBUG
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

	    (void)printf(" %+.6f%+.6fj", creal(v), cimag(v));
	}
	(void)printf("\n");
    }
    (void)printf("]\n");
}
#endif

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
    const int unknowns = vlp->vl_t_terms - 1;	/* unknowns per system */
    double frequency = vnp->vn_frequency_vector[findex];
    double *w_vector = NULL;
    double complex *prev_x_segment = NULL;
    int rv = -1;

    /*
     * If a measurement error vector was given, calculates weights
     * for each measurement and allocate the prev_x_segment.
     */
    if (vnp->vn_m_error_vector != NULL) {
	if ((w_vector = vs_calc_weights(vnssp)) == NULL) {
	    goto out;
	}
	prev_x_segment = malloc(unknowns * sizeof(double complex));
	if (prev_x_segment == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	    goto out;
	}
	_vnacal_new_solve_init_x_vector(vnssp, x_vector, x_length);
    }

    /*
     * For each system of equations...
     */
    assert(x_length == vnp->vn_systems * unknowns);
    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
	const int offset = sindex * unknowns;
	double complex *x_segment = &x_vector[offset];
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[sindex];
	const int equations = vnsp->vns_equation_count;
	int iteration = 0;

	/*
	 * Initialize prev_x_segment if in-use.
	 */
	if (prev_x_segment != NULL) {
	    (void)memcpy((void *)prev_x_segment, (void *)x_segment,
		         unknowns * sizeof(double complex));
	}

	/*
	 * For each iteration on the V matrices (if in use)...
	 */
	for (;;) {
	    double complex a_matrix[equations][unknowns];
	    double complex b_vector[equations];
	    double sum_dx_squared;
	    int eq_count = 0;

#ifdef DEBUG
	    if (vs_have_v(vnssp)) {
		int standard = 0;
		int v_rows, v_columns;
		vnacal_new_measurement_t *vnmp;

		if (VL_IS_T(vlp)) {
		    v_rows    = VL_S_COLUMNS(vlp);
		    v_columns = VL_M_COLUMNS(vlp);
		} else {
		    v_rows    = VL_M_ROWS(vlp);
		    v_columns = VL_S_ROWS(vlp);
		}
		for (vnmp = vnp->vn_measurement_list; vnmp != NULL;
			vnmp = vnmp->vnm_next) {
		    vnacal_new_msv_matrices_t *vnmmp;

		    vnmmp = &vnssp->vnss_msv_matrices[vnmp->vnm_index];
		    for (int sindex = 0; sindex < vnp->vn_systems; ++sindex) {
			char name[24];

			(void)sprintf(name, "v%d_%d", standard + 1, sindex + 1);
			print_cmatrix(name, vnmmp->vnsm_v_matrices[sindex],
				v_rows, v_columns);
		    }
		    ++standard;
		}
	    }
#endif /* DEBUG */

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

		determinant = _vnacommon_mldivide(x_segment,
			&a_matrix[0][0], b_vector, unknowns, 1);
		if (determinant == 0.0 || !isnormal(cabs(determinant))) {
		    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			    "singular linear system");
		    return -1;
		}
	    } else {
		int rank;

		rank = _vnacommon_qrsolve(x_segment, &a_matrix[0][0],
			b_vector, eq_count, unknowns, 1);
		if (rank < unknowns) {
		    _vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			    "singular linear system");
		    return -1;
		}
	    }
#ifdef DEBUG
	    {
		char buf[6 * sizeof(int) + 12];

		(void)sprintf(buf, "x[%d..%d]", offset, offset + unknowns - 1);
		print_cmatrix(buf, x_segment, unknowns, 1);
	    }
#endif /* DEBUG */

	    /*
	     * If measurement errors were not given or if this particular
	     * system is not over-determined, then we don't need to iterate.
	     */
	    if (!vs_have_v(vnssp)) {
		break;
	    }

	    /*
	     * x_segment depends on the V matrices and the V matrices
	     * depend on x_segment.  Iterate until they converge.
	     */
	    if (vs_update_v_matrices("vnacal_new_solve", vnssp, sindex,
			x_segment, unknowns) == -1) {
		return -1;
	    }
	    sum_dx_squared = 0.0;
	    for (int i = 0; i < unknowns; ++i) {
		double complex d = x_segment[i] - prev_x_segment[i];

		sum_dx_squared += _vnacommon_cabs2(d);
	    }
#ifdef DEBUG
	    (void)printf("RMS change in x_segment %e\n",
		    sqrt(sum_dx_squared / (double)unknowns));
#endif /* DEBUG */
	    if (sum_dx_squared / (double)unknowns <= vnp->vn_et_tolerance *
						     vnp->vn_et_tolerance) {
#ifdef DEBUG
		(void)printf("stop: converged\n");
#endif /* DEBUG */
		break;
	    }
	    if (++iteration >= vnp->vn_iteration_limit) {
		_vnacal_error(vcp, VNAERR_MATH, "vnacal_new_solve: "
			"measurement error model failed to converge at %e Hz",
			frequency);
		goto out;
	    }
	    (void)memcpy((void *)prev_x_segment, (void *)x_segment,
		    unknowns * sizeof(double complex));
	}
    }
    rv = 0;

out:
    free((void *)prev_x_segment);
    free((void *)w_vector);
    return rv;
}
