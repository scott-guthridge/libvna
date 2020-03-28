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
 * _vnacommon_mldivide: find X = A^-1
 *  @x: serialized result matrix (n x n)
 *  @a: serialized A matrix (n x n), destroyed on return
 *  @n: dimension of a and x
 *
 * Invert A storing the result in X.  Matrix A is destroyed.
 *
 * Returns the determinant.
 */
double complex _vnacommon_minverse(complex double *x, complex double *a, int n)
{
    int row_index[n];
    double complex d;

#define A(i, j)		(a[(i) * n + (j)])
#define X(i, j)		(x[(i) * n + (j)])

    /*
     * Replace A with its in-place LU decomposition.
     */
    d = _vnacommon_lu(a, row_index, n);

    /*
     * For each column...
     */
    for (int j = 0; j < n; ++j) {
	/*
	 * Use forward substitution to find the intermediate X' such that
	 * L X' = I.
	 */
	for (int i = 0; i < n; ++i) {
	    double complex s = (row_index[i] == j) ? 1.0 : 0.0;

	    for (int k = 0; k < i; ++k) {
		s -= A(i, k) * X(k, j);
		assert(i > k);
	    }
	    X(i, j) = s;
	}
	/*
	 * Use back substitution to find the result X, such that U X = X'.
	 */
	for (int i = n - 1; i >= 0; --i) {
	    double complex s = X(i, j);

	    for (int k = i + 1; k < n; ++k) {
		s -= A(i, k) * X(k, j);
	    }
	    X(i, j) = s / A(i, i);
	}
    }
    return d;
}
#undef A
#undef X
