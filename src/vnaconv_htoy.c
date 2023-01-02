/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2023 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * vnaconv_htoy: convert h-parameters to y-parameters
 */
void vnaconv_htoy(const double complex (*h)[2], double complex (*y)[2])
{
    const double complex h11 = h[0][0];
    const double complex h12 = h[0][1];
    const double complex h21 = h[1][0];
    const double complex h22 = h[1][1];
#define Y11	y[0][0]
#define Y12	y[0][1]
#define Y21	y[1][0]
#define Y22	y[1][1]

    Y11 =  1.0 / h11;
    Y12 = -h12 / h11;
    Y21 =  h21 / h11;
    Y22 =  (h11 * h22 - h12 * h21) / h11;
}
#undef Y22
#undef Y21
#undef Y12
#undef Y11
