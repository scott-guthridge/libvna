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
 * _vnacommon_mldivide: find X = A^-1 * B
 *  @x: serialized result matrix (m x n)
 *  @a: serialized A matrix (m x m), destroyed on return
 *  @b: serialized B matrix (m x n)
 *  @m: dimensions of A, number of rows in X and B
 *  @n: number of columns in X and B
 *
 * Divide matrix B by A from the left, storing the result in X.  Matrix A
 * is destroyed.
 *
 * Returns the determinant of A.
 */
double complex _vnacommon_mldivide(complex double *x, complex double *a,
	const double complex *b, int m, int n)
{
    int row_index[m];
    double complex d;

#define A(i, j)		(a[(i) * m + (j)])
#define B(i, j)		(b[(i) * n + (j)])
#define X(i, j)		(x[(i) * n + (j)])

    /*
     * Replace A with its in-place LU decomposition.
     */
    d = _vnacommon_lu(a, row_index, m);

    /*
     * For each column...
     */
    for (int j = 0; j < n; ++j) {
	/*
	 * Use forward substitution to find the intermediate X' such that
	 * L X' = B.
	 */
	for (int i = 0; i < m; ++i) {
	    double complex s = B(row_index[i], j);

	    for (int k = 0; k < i; ++k) {
		s -= A(i, k) * X(k, j);
		assert(i > k);
	    }
	    X(i, j) = s;
	}
	/*
	 * Use back substitution to find the result X, such that U X = X'.
	 */
	for (int i = m - 1; i >= 0; --i) {
	    double complex s = X(i, j);

	    for (int k = i + 1; k < m; ++k) {
		s -= A(i, k) * X(k, j);
	    }
	    X(i, j) = s / A(i, i);
	}
    }
    return d;
}
#undef A
#undef B
#undef X
