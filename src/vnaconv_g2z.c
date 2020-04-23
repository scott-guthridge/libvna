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
 * vnaconv_g2z: convert g-parameters to z-parameters
 */
void vnaconv_g2z(const vnaconv_array2_t *g, double complex (*z)[2])
{
    const double complex g11 = g[0][0];
    const double complex g12 = g[0][1];
    const double complex g21 = g[1][0];
    const double complex g22 = g[1][1];
#define Z11	z[0][0]
#define Z12	z[0][1]
#define Z21	z[1][0]
#define Z22	z[1][1]

    Z11 =  1.0 / g11;
    Z12 = -g12 / g11;
    Z21 =  g21 / g11;
    Z22 =  (g11 * g22 - g12 * g21) / g11;
}
#undef Z22
#undef Z21
#undef Z12
#undef Z11
