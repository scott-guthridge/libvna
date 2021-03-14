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
#include "vnaconv_internal.h"


/*
 * vnaconv_ztot: convert z-parameters to t-parameters
 */
void vnaconv_ztot(const vnaconv_array2_t *z, complex (*t)[2],
	       const double complex *z0)
{
    const double complex z11 = z[0][0];
    const double complex z12 = z[0][1];
    const double complex z21 = z[1][0];
    const double complex z22 = z[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = k1i / k2i * 2.0 * creal(z2) * z21;
#define T11	t[0][0]
#define T12	t[0][1]
#define T21	t[1][0]
#define T22	t[1][1]

    T11 =  (z12 * z21 - (z11 - z1c) * (z22 - z2c)) / d;
    T12 = -(z12 * z21 - (z11 - z1c) * (z22 + z2))  / d;
    T21 =  (z12 * z21 - (z11 + z1)  * (z22 - z2c)) / d;
    T22 = -(z12 * z21 - (z11 + z1)  * (z22 + z2))  / d;
}
#undef T22
#undef T21
#undef T12
#undef T11
