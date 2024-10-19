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

#include <complex.h>
#include <math.h>
#include "vnacommon_internal.h"
#include "vnaconv_internal.h"


/*
 * vnaconv_ztozin: calculate the two-port input port impedances
 *   @z:  serialized z matrix in (n x n)
 *   @zi: zin vector out (length n)
 *   @z0: reference impedance vector in (length n)
 *   @n:  length
 */
void vnaconv_ztozin(const double complex *z, double complex *zi,
	const double complex *z0, int n)
{
    if (n <= 0)
	return;

    double complex a[n * n];
    double complex x[n * n];
#define A(i, j)		(a[(i) * n + (j)])
#define X(i, j)		(x[(i) * n + (j)])
#define Z(i, j)		(z[(i) * n + (j)])

    /*
     * Find x = (z + diag(z0))^-1
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    A(i, j) = Z(i, j);
	}
	A(i, i) += z0[i];
    }
    _vnacommon_minverse(x, a, n);

    /*
     * Find zi[i] = 1 / x[i][i] - z0[i]
     */
    for (int i = 0; i < n; ++i) {
	zi[i] = 1.0 / X(i, i) - z0[i];
    }
}
#undef Z
#undef X
#undef A
