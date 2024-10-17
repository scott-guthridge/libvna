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
 * vnaconv_utosr: convert u-parameters to s-parameters, renormalizing
 */
void vnaconv_utosr(const double complex (*u)[2], double complex (*s)[2],
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
    const double kx = sqrt(fabs(creal(z11) * creal(z22) /
			       (creal(z12) * creal(z21))));
    const double complex z11mz21   = z11  - z21;
    const double complex z11pz21c  = z11  + z21c;
    const double complex z11cpz21  = z11c + z21;
    const double complex z11cmz21c = z11c - z21c;
    const double complex z12cpz22  = z12c + z22;
    const double complex z12cmz22c = z12c - z22c;
    const double complex z12mz22   = z12  - z22;
    const double complex z12pz22c  = z12  + z22c;
    const double complex z21pz11c  = z21  + z11c;
    const double complex z22pz12c  = z22  + z12c;
    const double complex d = - z22pz12c * (z21pz11c * u11 - z11mz21 * u12)
                             - z12mz22  * (z21pz11c * u21 - z11mz21 * u22);

#define S11	s[0][0]
#define S12	s[0][1]
#define S21	s[1][0]
#define S22	s[1][1]

    S11 = (-z12cpz22 * (z11cmz21c * u11 - z11pz21c * u12)
           - z12mz22 * (z11cmz21c * u21 - z11pz21c * u22)) / d;
    S12 = -4.0 * kx * creal(z12) * creal(z21) / d;
    S21 = -4.0 / kx * creal(z11) * creal(z22) * (u11 * u22 - u12 * u21) / d;
    S22 = (-z12cmz22c * (z11cpz21 * u11 - z11mz21 * u12)
           - z12pz22c * (z11cpz21 * u21 - z11mz21 * u22)) / d;

}
#undef S22
#undef S21
#undef S12
#undef S11
