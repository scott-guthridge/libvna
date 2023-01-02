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
 * _vnacommon_qrd: find the QR decomposition of A, destroying A
 *  @a:       mxn serialized input matrix (also output matrix)
 *  @d:       MIN(m, n) length vector to receive major diagonal of R
 *  @rows:    number of rows in A
 *  @columns: number of columns in A
 *
 * On return, the lower triangle of A (including the major diagonal)
 * contains the v_i vectors that can be used to construct Q:
 *
 *   Q = (I - 2 v1 v1')' * (I - 2 v2 v2')' + ... (I - 2 vn vn')'
 *
 * Since v1, v2, v3, etc. become progressively shorter, one must
 * treat the rows above as containing zeros.  We don't actually
 * calculate Q in this function, though -- only the v_j vectors.
 *
 * The upper triangle of A not including the diagonal contains R,
 * with the diagonal terms placed in the d vector.
 */
void _vnacommon_qrd(complex double *a, complex double *d, int rows, int columns)
{
    int diagonals = MIN(rows, columns);

#define A(i, j)		(a[(i) * columns + (j)])

    /*
     * For each diagonal element in A...
     */
    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
	double norm;
	double subdot = 0.0;	/* common subexpression factored */
	double complex alpha;	/* new diagonal value */

	/*
	 * Let v be a(diagonal:columns-1, diagonal).
	 *
	 * Calculate the sub-expression: v(2:)' * v(2:), needed twice
	 * below.
	 */
	for (int row = diagonal + 1; row < rows; ++row) {
	    subdot += _vnacommon_cabs2(A(row, diagonal));
	}

	/*
	 * Compute alpha with norm(v) and angle the opposite direction
	 * as v(1).  This will be the next diagonal term of R.
	 */
	alpha = -cexp(I * carg(A(diagonal, diagonal))) *
	    sqrt(_vnacommon_cabs2(A(diagonal, diagonal)) + subdot);
	d[diagonal] = alpha;

	/*
	 * v = v - alpha * [1, 0, 0, ... 0]'.
	 */
	A(diagonal, diagonal) -= alpha;

	/*
	 * v = v / norm(v).
	 */
	norm = sqrt(_vnacommon_cabs2(A(diagonal, diagonal)) + subdot);
	for (int row = diagonal; row < rows; ++row) {
	    A(row, diagonal) /= norm;
	}

	/*
	 * Multiply R on the left by I - 2 v v'.
	 */
	for (int column = diagonal + 1; column < columns; ++column) {
	    double complex temp = 0.0;

	    for (int row = diagonal; row < rows; ++row) {
		temp += conj(A(row, diagonal)) * A(row, column);
	    }
	    for (int row = diagonal; row < rows; ++row) {
		A(row, column) -= 2.0 * temp * A(row, diagonal);
	    }
	}
    }
}
#undef A
