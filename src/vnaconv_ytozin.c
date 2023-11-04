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
 * vnaconv_ytozin: calculate the two-port input port impedances
 *   @y:  serialized z matrix in (n x n)
 *   @zi: zin vector out (length n)
 *   @z0: system impedance vector in (length n)
 *   @n:  length
 */
void vnaconv_ytozin(const double complex *y, double complex *zi,
	const double complex *z0, int n)
{
    if (n <= 0)
	return;

    double complex a[n * n];
    double complex b[n * n];
    double complex s[n * n];
#define A(i, j)         (a[(i) * n + (j)])
#define B(i, j)         (b[(i) * n + (j)])
#define S(i, j)         (s[(i) * n + (j)])
#define Y(i, j)         (y[(i) * n + (j)])

    /*
     * Find:
     *   b  = I - z0* y
     *   a  = I + z0  y
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    B(i, j) = -conj(z0[i]) * Y(i, j);
	    A(i, j) = z0[i] * Y(i, j);
	}
	B(i, i) += 1.0;
	A(i, i) += 1.0;
    }

    /*
     * Find s = b a^-1.  Note that the result isn't quite "s", since
     * we didn't multiply (ki . s * k), but we want only the major
     * diagonal, so we don't need to do that extra step.
     */
    _vnacommon_mrdivide(s, b, a, n, n);

    /*
     * Find zi from the major diagonal of s.
     */
    for (int i = 0; i < n; ++i) {
	const double complex sii = S(i, i);

	zi[i] = (sii * z0[i] + conj(z0[i])) / (1.0 - sii);
    }
}
