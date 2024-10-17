/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2024 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * vnaconv_stosrn: renormalize s-parameters (n-port)
 *   @si:  given serialized nxn s-parameter matrix
 *   @so:  caller-allocated resulting serialized nxn s-parameter matrix
 *   @z1: vector of initial impedances for each port
 *   @z2: vector of final impedances for each port
 *   @n:  dimension
 */
void vnaconv_stosrn(const double complex *si, double complex *so,
	const double complex *z1, const double complex *z2, int n)
{
    if (n <= 0)
	return;

    double complex a[n * n];
    double complex b[n * n];
#define A(i, j)         (a[(i) * n + (j)])
#define B(i, j)         (b[(i) * n + (j)])
#define SI(i, j)        (si[(i) * n + (j)])
#define SO(i, j)        (zo[(i) * n + (j)])

    /*
     * Find:
     *   K = diag(0.5 sqrt(fabs(creal(z1) ./ creal(z2))) ./ creal(z1))
     *   VM = diag(conj(z1)) + SI * diag(z1)
     *   IM = I - SI
     *   A = K * (VM + diag(z2) * IM)
     *   B = K * (VM - diag(conj(z2)) * IM)
     */
    for (int r = 0; r < n; ++r) {
	double rz1 = creal(z1[r]);
	double rz2 = creal(z2[r]);
	double k = 0.5 * sqrt(fabs(rz1 / rz2)) / rz1;

	for (int c = 0; c < n; ++c) {
	    double complex v, i;

	    if (r == c) {
		v = conj(z1[r]);
		i = 1.0;
	    } else {
		v = 0.0;
		i = 0.0;
	    }
	    v += z1[r] * SI(r, c);
	    i -= SI(r, c);
	    A(r, c) = k * (v + z2[r] * i);
	    B(r, c) = k * (v - conj(z2[r]) * i);
	}
    }

    /*
     * Find SO = B A^-1
     */
    _vnacommon_mrdivide(so, b, a, n, n);
}
#undef SO
#undef SI
#undef B
#undef A
