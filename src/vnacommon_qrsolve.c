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
 * _vnacommon_qrsolve: solve the system A X = B
 *  @a:       mxn serialized coefficient matrix (destroyed)
 *  @b:       mxo constant term matrix (destroyed)
 *  @x:	      nxo result matrix
 *  @m: number of rows in A and B
 *  @n: number of columns in A, and rows in X
 *  @o: number of columns in B and X
 *
 * Solves the system of equations using QR decomposition, where A doesn't
 * have to be square.  If A has more columns than rows (underdetermined
 * case), the function finds a solution with excess variables set to zero.
 * If A has more rows than columns, (overdetermined case), the function
 * returns a solution that that minimizes error in a least-squares sense.
 *
 * Note: both a and b are destroyed!
 */
void _vnacommon_qrsolve(complex double *a, complex double *b,
	complex double *x, int m, int n, int o)
{
    int diagonals = MIN(m, n);
    complex double d[diagonals];

#define A(i, j)		((a)[(i) * n + (j)])
#define X(i, j)		((x)[(i) * o + (j)])
#define B(i, j)		((b)[(i) * o + (j)])

    /*
     * Find the QR decomposition of A.  On return the lower triangle
     * of A is replaced with the v_i vectors used to construct Q, the
     * upper triangle (above the diagonal) contains R, and d contains
     * the diagonal terms of R.
     */
    _vnacommon_qrd(a, d, m, n);

    /*
     * For each column in X,B...
     */
    for (int k = 0; k < o; ++k) {
	/*
	 * Multiply B_{i,k} on the left by Q_1, Q_2, Q_3, ... Q_diagonals
	 * in order where Q_n = I - 2 * v_n * v_n'
	 */
	for (int i = 0; i < diagonals; ++i) {
	    double complex s = 0.0;

	    for (int j = i; j < m; ++j) {
		s += conj(A(j, i)) * B(j, k);
	    }
	    for (int j = i; j < m; ++j) {
		B(j, k) -= 2 * s * A(j, i);
	    }
	}

	/*
	 * If there are more unknowns than equations (underdetermined
	 * case), set the excess unknowns to zero.
	 */
	for (int i = n - 1; i >= diagonals; --i) {
	    X(i, k) = 0.0;
	}

	/*
	 * Use back substitution to find X, where R * X = Q' * B.
	 */
	for (int i = diagonals - 1; i >= 0; --i) { /* foreach diagonal... */
	    double complex s = 0.0;

	    for (int j = i + 1; j < diagonals; ++j) { /* foreach solved var */
		s += A(i, j) * X(j, k);
	    }
	    X(i, k) = (B(i, k) - s) / d[i];
	}
    }
}
#undef B
#undef X
#undef A
