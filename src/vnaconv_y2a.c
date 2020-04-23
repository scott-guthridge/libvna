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
 * vnaconv_y2a: convert y-parameters to a-parameters
 */
void vnaconv_y2a(const vnaconv_array2_t *y, double complex (*a)[2])
{
    const double complex y11 = y[0][0];
    const double complex y12 = y[0][1];
    const double complex y21 = y[1][0];
    const double complex y22 = y[1][1];
#define A11	a[0][0]
#define A12	a[0][1]
#define A21	a[1][0]
#define A22	a[1][1]

    A11 = -y22 / y21;
    A12 = -1.0 / y21;
    A21 = (y12 * y21 - y11 * y22) / y21;
    A22 = -y11 / y21;
}
#undef A22
#undef A21
#undef A12
#undef A11
