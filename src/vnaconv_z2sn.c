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

#include <math.h>
#include "vnacommon_internal.h"
#include "vnaconv_internal.h"


/*
 * vnaconv_z2sn: convert z-parameters to s-parameters (n-port)
 *   @z:  given serialized nxn z-parameter matrix
 *   @s:  caller-allocated resulting serialized nxn s-parameter matrix
 *   @z0: vector of impedances seen by each port
 *   @n:  dimension
 */
void vnaconv_z2sn(const double complex *z, double complex *s,
	       const double complex *z0, int n)
{
    double complex a[n * n];
    double complex b[n * n];
    double ki[n];
#define A(i, j)         (a[(i) * n + (j)])
#define B(i, j)         (b[(i) * n + (j)])
#define S(i, j)         (s[(i) * n + (j)])
#define Z(i, j)         (z[(i) * n + (j)])

    /*
     * Find:
     *   b  = z - z0*
     *   a  = z + z0
     *   ki = square root of absolute value of real of z0
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    B(i, j) = Z(i, j);
	    A(i, j) = Z(i, j);
	}
	B(i, i) -= conj(z0[i]);
	A(i, i) += z0[i];
	ki[i] = sqrt(fabs(creal(z0[i])));
    }

    /*
     * Find s = b a^-1
     */
    _vnacommon_mrdivide(s, b, a, n, n);

    /*
     * Find s = diag(ki^-1) s diag(ki)
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    if (i != j) {
		S(i, j) *= ki[j] / ki[i];
	    }
	}
    }
}
#undef Y
#undef S
#undef B
#undef A
