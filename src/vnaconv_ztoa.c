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
 * vnaconv_ztoa: convert z-parameters to a-parameters
 */
void vnaconv_ztoa(const vnaconv_array2_t *z, double complex (*a)[2])
{
    const double complex z11 = z[0][0];
    const double complex z12 = z[0][1];
    const double complex z21 = z[1][0];
    const double complex z22 = z[1][1];
#define A11	a[0][0]
#define A12	a[0][1]
#define A21	a[1][0]
#define A22	a[1][1]

    A11 = z11 / z21;
    A12 = (z11 * z22 - z12 * z21) / z21;
    A21 = 1.0 / z21;
    A22 = z22 / z21;
}
#undef A22
#undef A21
#undef A12
#undef A11
