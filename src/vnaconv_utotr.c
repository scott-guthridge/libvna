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
#include "vnaconv_internal.h"


/*
 * vnaconv_utotr: convert u-parameters to t-parameters, renormalizing
 */
void vnaconv_utotr(const double complex (*u)[2], double complex (*t)[2],
	           const double complex *z1, const double complex *z2)
{
    const double complex u11 = u[0][0];
    const double complex u12 = u[0][1];
    const double complex u21 = u[1][0];
    const double complex u22 = u[1][1];
    const double complex z11 = z1[0];
    const double complex z12 = z1[1];
    const double complex z21 = z2[0];
    const double complex z22 = z2[1];
    const double complex z11c = conj(z11);
    const double complex z12c = conj(z12);
    const double complex z21c = conj(z21);
    const double complex z22c = conj(z22);
    const double kx = sqrt(fabs(creal(z12) * creal(z21) /
			       (creal(z11) * creal(z22))));
    const double complex z11pz21c  = z11  + z21c;
    const double complex z11cpz21  = z11c + z21;
    const double complex z11cmz21c = z11c - z21c;
    const double complex z12mz22   = z12  - z22;
    const double complex z12pz22c  = z12  + z22c;
    const double complex z12cpz22  = z12c + z22;
    const double complex z12cmz22c = z12c - z22c;
    const double complex d = -4.0 * kx * creal(z11) * creal(z22) *
                             (u11 * u22 - u12 * u21);

#define T11	t[0][0]
#define T12	t[0][1]
#define T21	t[1][0]
#define T22	t[1][1]

    T11 = (z12cmz22c  * (z11cmz21c * u11 - z11pz21c * u12)
           + z12pz22c * (z11cmz21c * u21 - z11pz21c * u22)) / d;
    T12 = (-z12cpz22 * (z11cmz21c * u11 - z11pz21c * u12)
           - z12mz22 * (z11cmz21c * u21 - z11pz21c * u22)) / d;
    T21 = (z12cmz22c  * (z11cpz21 * u11 - (z11 - z21) * u12)
           + z12pz22c * (z11cpz21 * u21 - (z11 - z21) * u22)) / d;
    T22 = (-z12cpz22 * (z11cpz21 * u11 - (z11 - z21) * u12)
           - z12mz22 * (z11cpz21 * u21 - (z11 - z21) * u22)) / d;
}
#undef T22
#undef T21
#undef T12
#undef T11
