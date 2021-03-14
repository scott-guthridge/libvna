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
 * vnaconv_ytog: convert y-parameters to g-parameters
 */
void vnaconv_ytog(const vnaconv_array2_t *y, double complex (*g)[2])
{
    const double complex y11 = y[0][0];
    const double complex y12 = y[0][1];
    const double complex y21 = y[1][0];
    const double complex y22 = y[1][1];
#define G11	g[0][0]
#define G12	g[0][1]
#define G21	g[1][0]
#define G22	g[1][1]

    G11 = (y11 * y22 - y12 * y21) / y22;
    G12 =  y12 / y22;
    G21 = -y21 / y22;
    G22 =  1.0 / y22;
}
#undef G22
#undef G21
#undef G12
#undef G11
