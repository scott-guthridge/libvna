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
 * vnaconv_g2y: convert g-parameters to y-parameters
 */
void vnaconv_g2y(const vnaconv_array2_t *g, double complex (*y)[2])
{
    const double complex g11 = g[0][0];
    const double complex g12 = g[0][1];
    const double complex g21 = g[1][0];
    const double complex g22 = g[1][1];
#define Y11	y[0][0]
#define Y12	y[0][1]
#define Y21	y[1][0]
#define Y22	y[1][1]

    Y11 =  (g11 * g22 - g12 * g21) / g22;
    Y12 =  g12 / g22;
    Y21 = -g21 / g22;
    Y22 =  1.0 / g22;
}
#undef Y22
#undef Y21
#undef Y12
#undef Y11
