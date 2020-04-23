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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_a2t: convert a-parameters to t-parameters
 */
void vnaconv_a2t(const vnaconv_array2_t *a, complex (*t)[2],
	       const double complex *z0)
{
    const double complex a11 = a[0][0];
    const double complex a12 = a[0][1];
    const double complex a21 = a[1][0];
    const double complex a22 = a[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = k1i / k2i * 2.0 * creal(z2);
#define T11	t[0][0]
#define T12	t[0][1]
#define T21	t[1][0]
#define T22	t[1][1]

    T11 = -(a12 - a22 * z1c - a11 * z2c + a21 * z1c * z2c) / d;
    T12 =  (a12 - a22 * z1c + a11 * z2  - a21 * z1c * z2)  / d;
    T21 = -(a12 + a22 * z1  - a11 * z2c - a21 * z1  * z2c) / d;
    T22 =  (a12 + a22 * z1  + a11 * z2  + a21 * z1  * z2)  / d;
}
#undef T22
#undef T21
#undef T12
#undef T11
