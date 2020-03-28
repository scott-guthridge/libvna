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
 * vnaconv_b2g: convert b-parameters to g-parameters
 */
void vnaconv_b2g(const vnaconv_array2_t *b, double complex (*g)[2])
{
    const double complex b11 = b[0][0];
    const double complex b12 = b[0][1];
    const double complex b21 = b[1][0];
    const double complex b22 = b[1][1];
#define G11	g[0][0]
#define G12	g[0][1]
#define G21	g[1][0]
#define G22	g[1][1]

    G11 = -b21 / b22;
    G12 = -1.0 / b22;
    G21 =  (b11 * b22 - b12 * b21) / b22;
    G22 = -b12 / b22;
}
#undef G22
#undef G21
#undef G12
#undef G11
