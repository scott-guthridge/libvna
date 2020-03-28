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
 * vnaconv_a2h: convert a-parameters to h-parameters
 */
void vnaconv_a2h(const vnaconv_array2_t *a, double complex (*h)[2])
{
    const double complex a11 = a[0][0];
    const double complex a12 = a[0][1];
    const double complex a21 = a[1][0];
    const double complex a22 = a[1][1];
#define H11	h[0][0]
#define H12	h[0][1]
#define H21	h[1][0]
#define H22	h[1][1]

    H11 =  a12 / a22;
    H12 =  (a11 * a22 - a12 * a21) / a22;
    H21 = -1.0 / a22;
    H22 =  a21 / a22;
}
#undef H22
#undef H21
#undef H12
#undef H11
