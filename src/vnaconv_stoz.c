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
#include "vnaconv_internal.h"


/*
 * vnaconv_stoz: convert s-parameters to z-parameters
 */
void vnaconv_stoz(const double complex (*s)[2], double complex (*z)[2],
	       const double complex *z0)
{
    const double complex s11 = s[0][0];
    const double complex s12 = s[0][1];
    const double complex s21 = s[1][0];
    const double complex s22 = s[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex dz = s11 * s22 - s12 * s21;
    const double complex d = 1.0 - s11 + dz - s22;
#define Z11	z[0][0]
#define Z12	z[0][1]
#define Z21	z[1][0]
#define Z22	z[1][1]

    Z11 = -((dz - s11) * z1 + (s22 - 1.0) * z1c) / d;
    Z12 =  (k1i / k2i * s12 * (z2 + z2c)) / d;
    Z21 =  (k2i / k1i * s21 * (z1 + z1c)) / d;
    Z22 = -((dz - s22) * z2 + (s11 - 1.0) * z2c) / d;
}
#undef Z22
#undef Z21
#undef Z12
#undef Z11
