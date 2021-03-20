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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacommon_internal.h"

/*
 * _vnacommon_qrsolve_q: solve the system A X = B and return Q
 *  @x:	nxo result matrix
 *  @a: mxn serialized coefficient matrix (destroyed)
 *  @b: mxo constant term matrix (destroyed)
 *  @q: mxm matrix to receive Q
 *  @m: number of rows in A and B
 *  @n: number of columns in A, and rows in X
 *  @o: number of columns in B and X
 *
 * Solves the system of equations using QR decomposition and returns the
 * Q matrix.  If A has more columns than rows (underdetermined case),
 * the function finds a solution with excess variables set to zero.
 * If A has more rows than columns, (overdetermined case), the function
 * finds a solution that that minimizes error in a least-squares sense.
 *
 * Note: both a and b are destroyed!
 *
 * Returns the rank.
 */
int _vnacommon_qrsolve_q(complex double *x, complex double *a,
	complex double *b, complex double *q, int m, int n, int o)
{
    int diagonals = MIN(m, n);
    int rank;

#define Q(i, j)		((q)[(i) * m + (j)])
#define A(i, j)		((a)[(i) * n + (j)])

    /*
     * Use _vnacommon_qrsolve to solve the system.  On return, the
     * lower triangle of the "a" matrix contains the vectors we need
     * to build Q.
     */
    rank = _vnacommon_qrsolve(x, a, b, m, n, o);

    /*
     * Initialize q to the identity matrix.
     */
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < m; ++j) {
	    Q(i, j) = (i == j) ? 1.0 : 0.0;
	}
    }

    /*
     * Form Q.
     */
    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
	for (int i = 0; i < m; ++i) {
	    double complex s = 0.0;

	    for (int j = diagonal; j < m; ++j) {
		s += Q(i, j) * A(j, diagonal);
	    }
	    for (int j = diagonal; j < m; ++j) {
		Q(i, j) -= 2.0 * s * conj(A(j, diagonal));
	    }
	}
    }

    return rank;
}
#undef A
#undef Q
