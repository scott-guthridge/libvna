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
 * vnaconv_ttoz: convert t-parameters to z-parameters
 */
void vnaconv_ttoz(const vnaconv_array2_t *t, complex (*z)[2],
	       const double complex *z0)
{
    const double complex t11 = t[0][0];
    const double complex t12 = t[0][1];
    const double complex t21 = t[1][0];
    const double complex t22 = t[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = -t11 - t12 + t21 + t22;
#define Z11	z[0][0]
#define Z12	z[0][1]
#define Z21	z[1][0]
#define Z22	z[1][1]

    Z11 =             ((t11 + t12) * z1 + (t21 + t22) * z1c)    / d;
    Z12 = k1i / k2i * (t11 * t22 - t12 * t21) * 2.0 * creal(z2) / d;
    Z21 = k2i / k1i *                           2.0 * creal(z1) / d;
    Z22 =             ((t11 - t21) * z2 + (t22 - t12) * z2c)    / d;
}
#undef Z22
#undef Z21
#undef Z12
#undef Z11
