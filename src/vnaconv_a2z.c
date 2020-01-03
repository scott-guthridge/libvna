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
 * vnaconv_a2z: convert a-parameters to z-parameters
 */
void vnaconv_a2z(const vnaconv_array2_t *a, double complex (*z)[2])
{
    const double complex a11 = a[0][0];
    const double complex a12 = a[0][1];
    const double complex a21 = a[1][0];
    const double complex a22 = a[1][1];
#define Z11	z[0][0]
#define Z12	z[0][1]
#define Z21	z[1][0]
#define Z22	z[1][1]

    Z11 = a11 / a21;
    Z12 = (a11 * a22 - a12 * a21) / a21;
    Z21 = 1.0 / a21;
    Z22 = a22 / a21;
}
#undef Z22
#undef Z21
#undef Z12
#undef Z11
