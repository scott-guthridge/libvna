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
 * MERCHANTABILITY or FITNESS FOR A11 PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <math.h>
#include "vnacommon_internal.h"
#include "vnaconv_internal.h"

/*
 * vnaconv_y2sn: convert y-parameters to s-parameters (n-port)
 *   @y:  given serialized nxn y-parameter matrix
 *   @s:  caller-allocated resulting serialized nxn s-parameter matrix
 *   @z0: vector of impedances seen by each port
 *   @n:  dimension
 */
void vnaconv_y2sn(const double complex *y, double complex *s,
	       const double complex *z0, int n)
{
    double complex a[n * n];
    double complex b[n * n];
    double k[n];
#define A(i, j)         (a[(i) * n + (j)])
#define B(i, j)         (b[(i) * n + (j)])
#define S(i, j)         (s[(i) * n + (j)])
#define Y(i, j)         (y[(i) * n + (j)])

    /*
     * Find:
     *   b  = I - z0* y
     *   a  = I + z0  y
     *   k = square root of absolute value of real of z0
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    B(i, j) = -conj(z0[i]) * Y(i, j);
	    A(i, j) = z0[i] * Y(i, j);
	}
	B(i, i) += 1.0;
	A(i, i) += 1.0;
	k[i] = sqrt(fabs(creal(z0[i])));
    }

    /*
     * Find s = b a^-1
     */
    _vnacommon_mrdivide(s, b, a, n, n);

    /*
     * Find s = diag(k^-1) s diag(k)
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    if (i != j) {
		S(i, j) *= k[j] / k[i];
	    }
	}
    }
}
#undef Z
#undef S
#undef B
#undef A
