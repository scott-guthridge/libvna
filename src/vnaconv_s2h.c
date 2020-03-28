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
#include "vnaconv_internal.h"


/*
 * vnaconv_s2h: convert s-parameters to h-parameters
 */
void vnaconv_s2h(const vnaconv_array2_t *s, complex (*h)[2],
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
    const double complex dh = s11 * s22 - s12 * s21;
    const double complex d = (dh - s22) * z2 + (s11 - 1.0) * z2c;
#define H11	h[0][0]
#define H12	h[0][1]
#define H21	h[1][0]
#define H22	h[1][1]

    H11 = -(dh * z1 * z2 + s11 * z1 * z2c + s22 * z1c * z2 + z1c * z2c) / d;
    H12 = -k1i / k2i * s12 * (z2 + z2c) / d;
    H21 =  k2i / k1i * s21 * (z1 + z1c) / d;
    H22 = -(1.0 + dh - s11 - s22) / d;
}
#undef H22
#undef H21
#undef H12
#undef H11
