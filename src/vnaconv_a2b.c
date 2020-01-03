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
 * vnaconv_a2b: convert a-parameters to b-parameters
 */
void vnaconv_a2b(const vnaconv_array2_t *a, double complex (*b)[2])
{
    const double complex a11 = a[0][0];
    const double complex a12 = a[0][1];
    const double complex a21 = a[1][0];
    const double complex a22 = a[1][1];
    const double complex d = a11 * a22 - a12 * a21;
#define B11	b[0][0]
#define B12	b[0][1]
#define B21	b[1][0]
#define B22	b[1][1]

    B11 =  a22 / d;
    B12 = -a12 / d;
    B21 = -a21 / d;
    B22 =  a11 / d;
}
#undef B22
#undef B21
#undef B12
#undef B11
