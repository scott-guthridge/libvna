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

	    (void)printf(" %9.5f%+9.5fj", creal(v), cimag(v));
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
