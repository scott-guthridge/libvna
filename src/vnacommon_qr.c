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
 * _vnacommon_qr: find the QR decomposition of A
 *  @a: mxn serialized coefficient matrix (destroyed)
 *  @q: mxm matrix to receive Q
 *  @r: mxn matrix to receive R
 *  @m: number of rows in A
 *  @n: number of columns in A
 *
 * Note: a is destroyed!
 *
 * Returns the rank.
 */
int _vnacommon_qr(complex double *a, complex double *q, complex double *r,
	int m, int n)
{
    int diagonals = MIN(m, n);
    complex double d[diagonals];
    int rank = 0;

#define Q(i, j)		((q)[(i) * m + (j)])
#define R(i, j)		((r)[(i) * n + (j)])
#define A(i, j)		((a)[(i) * n + (j)])

    /*
     * Find the QR decomposition of A.  On return the lower triangle
     * of A is replaced with the v_i vectors used to construct Q, the
     * upper triangle (above the diagonal) contains the portion of R
     * above the major diagonal, and d contains the major diagonal
     * of R.
     */
    _vnacommon_qrd(a, d, m, n);

    /*
     * Initialize q to the identity matrix, and r to zero.
     */
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < m; ++j) {
	    Q(i, j) = (i == j) ? 1.0 : 0.0;
	}
	for (int j = 0; j < n; ++j) {
	    R(i, j) = 0.0;
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

    /*
     * Form R.
     */
    for (int i = 0; i < diagonals; ++i) {
	R(i, i) = d[i];
	(void)memcpy((void *)&R(i, i + 1), (void *)&A(i, i + 1),
		(n - i - 1) * sizeof(double complex));
    }

    /*
     * Find the rank.
     */
    for (int i = 0; i < diagonals; ++i) {
	if (isnormal(cabs(d[i])) && d[i] != 0.0) {
	    ++rank;
	}
    }
    return rank;
}
#undef A
#undef R
#undef Q
