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
 * vnaconv_stozn: convert s-parameters to z-parameters (n-port)
 *   @s:  given serialized nxn s-parameter matrix
 *   @z:  caller-allocated resulting serialized nxn z-parameter matrix
 *   @z0: vector of impedances seen by each port
 *   @n:  dimension
 */
void vnaconv_stozn(const double complex *s, double complex *z,
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
     *   a  = I - S
     *   b  = diag(z0*) + s diag(z0)
     *   ki = square root of absolute value of real of z0
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    A(i, j) = -S(i, j);
	    B(i, j) = S(i, j) * z0[j];
	}
	A(i, i) += 1.0;
	B(i, i) += conj(z0[i]);
	ki[i] = sqrt(fabs(creal(z0[i])));
    }

    /*
     * Find z = a^-1 b
     */
    _vnacommon_mldivide(z, a, b, n, n);

    /*
     * Find z = diag(ki) z diag(ki^-1)
     */
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
	    if (i != j) {
		Z(i, j) *= ki[i] / ki[j];
	    }
	}
    }
}
#undef Z
#undef S
#undef B
#undef A
