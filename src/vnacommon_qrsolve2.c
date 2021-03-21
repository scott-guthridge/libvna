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
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacommon_internal.h"

/*
 * _vnacommon_qrsolve2: solve the system Q R X = B
 *  @x:	nxo result matrix
 *  @q: mxm orthogonal matrix
 *  @r: mxn upper-triangular matrix
 *  @b: mxo constant term matrix
 *  @m: number of rows and columns in Q, number of rows in R, B
 *  @n: number of columns in A, R, and number of rows in X
 *  @o: number of columns in B and X
 *
 * Solves the system of equations given the QR decomposition of the
 * coefficient matrix.  If R has more columns than rows (underdetermined
 * case), the function finds a solution with excess variables set to zero.
 * If R has more rows than columns, (overdetermined case), the function
 * finds a solution that that minimizes error in a least-squares sense.
 */
void _vnacommon_qrsolve2(double complex *x, const double complex *q,
	const double complex *r, const double complex *b,
	const int m, const int n, const int o)
{
    int diagonals = MIN(m, n);

#define Q(i, j)		((q)[(i) * m + (j)])
#define R(i, j)		((r)[(i) * n + (j)])
#define X(i, j)		((x)[(i) * o + (j)])
#define B(i, j)		((b)[(i) * o + (j)])

    /*
     * For each column in X,B, use back-substitution to find X,
     * where R X = T.
     */
    for (int j = 0; j < o; ++j) {
	for (int i = diagonals - 1; i >= 0; --i) {
	    double complex s = 0.0;

	    /*
	     * Find T_ij where T = Q' B.
	     */
	    for (int k = 0; k < m; ++k) {
		s += conj(Q(k, i)) * B(k, j);
	    }

	    /*
	     * Use back subtitution to find X, where R X = T.
	     */
	    for (int k = i + 1; k < diagonals; ++k) {
		s -= R(i, k) * X(k, j);
	    }
	    X(i, j) = s / R(i, i);
	}
    }

    /*
     * If the system is undetermined, set the remaining X's to
     * zero.
     */
    for (int i = m; i < n; ++i) {
	for (int k = 0; k < o; ++k) {
	    X(i, k) = 0.0;
	}
    }
}
#undef B
#undef X
#undef R
#undef Q
