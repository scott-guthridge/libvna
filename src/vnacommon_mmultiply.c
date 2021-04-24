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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacommon_internal.h"


/*
 * _vnacommon_mmultiply: find C = A x B
 *   @c: serialized result matrix, m x o
 *   @a: serialized A matrix, m x n
 *   @b: serialized B matrix, n x o
 *   @m: first dimension of C and A
 *   @n: second dimension of A, first dimension of B
 *   @o: second dimension of C and B
 */
void _vnacommon_mmultiply(double complex *c, const double complex *a,
	const double complex *b, int m, int n, int o)
{
#define A(i, j) (a[(i) * n + (j)])
#define B(i, j) (b[(i) * o + (j)])
#define C(i, j) (c[(i) * o + (j)])

    for (int i = 0; i < m; ++i) {
	for (int k = 0; k < o; ++k) {
	    double complex s = 0.0;

	    for (int j = 0; j < n; ++j) {
		s += A(i, j) * B(j, k);
	    }
	    C(i, k) = s;
	}
    }
}
#undef A
#undef B
#undef C
