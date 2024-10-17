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
 * vnaconv_stosr: renormalize s-parameters
 */
void vnaconv_stosr(const double complex (*si)[2], double complex (*so)[2],
	           const double complex *z1, const double complex *z2)
{
    const double complex s11 = si[0][0];
    const double complex s12 = si[0][1];
    const double complex s21 = si[1][0];
    const double complex s22 = si[1][1];
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
    const double complex z12mz22   = z12  - z22;
    const double complex z12pz22c  = z12  + z22c;
    const double complex z12cpz22  = z12c + z22;
    const double complex z12cmz22c = z12c - z22c;

    const double complex d = (z11cpz21 + z11mz21 * s11) *
                             (z12cpz22 + z12mz22 * s22) -
			     z11mz21 * z12mz22 * s12 * s21;

#define S11	so[0][0]
#define S12	so[0][1]
#define S21	so[1][0]
#define S22	so[1][1]

    S11 = ((z11cmz21c + z11pz21c * s11) * (z12cpz22  + z12mz22  * s22)
           - z11pz21c * z12mz22 * s12 * s21) / d;
    S12 = 4.0 * kx * creal(z12) * creal(z21) * s12 / d;
    S21 = 4.0 / kx * creal(z11) * creal(z22) * s21 / d;
    S22 = ((z11cpz21  + z11mz21  * s11) * (z12cmz22c + z12pz22c * s22)
           - z11mz21  * z12pz22c * s12 * s21) / d;
}
#undef S22
#undef S21
#undef S12
#undef S11
