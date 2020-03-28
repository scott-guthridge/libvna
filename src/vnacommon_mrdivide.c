/*
 * Electrical Network Parameter Conversion Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * _vnacommon_mrdivide: find X = B * A^-1
 *  @x: serialized result matrix (m x n)
 *  @b: serialized B matrix (m x n)
 *  @a: serialized A matrix (n x n), destroyed on return
 *  @m: dimensions of A, number of rows in X and B
 *  @n: number of columns in X and B
 *
 * Divide matrix B by A from the right, storing the result in X.
 * Matrix A is destroyed.
 *
 * Returns the determinant of A.
 */
double complex _vnacommon_mrdivide(complex double *x, const complex double *b,
	double complex *a, int m, int n)
{
    int row_index[n];
    double complex d;

#define A(i, j)		(a[(i) * n + (j)])
#define B(i, j)		(b[(i) * n + (j)])
#define X(i, j)		(x[(i) * n + (j)])
    /*
     * Replace A with its in-place LU decomposition.
     */
    d = _vnacommon_lu(a, row_index, n);

    /*
     * For each row...
     */
    for (int i = 0; i < m; ++i) {
	/*
	 * Use back substitution to find the intermediate X' such that
	 * X' U = B.
	 */
	for (int j = 0; j < n; ++j) {
	    double complex s = B(i, j);

	    for (int k = 0; k < j; ++k) {
		s -= A(k, j) * X(i, row_index[k]);
	    }
	    X(i, row_index[j]) = s / A(j, j);
	}
	/*
	 * Use forward substitution to find X such that X * U = X'.
	 */
	for (int j = n - 1; j >= 0; --j) {
	    double complex s = X(i, row_index[j]);

	    for (int k = j + 1; k < n; ++k) {
		s -= A(k, j) * X(i, row_index[k]);
	    }
	    X(i, row_index[j]) = s;
	}
    }
    return d;
}
#undef A
#undef B
#undef X
