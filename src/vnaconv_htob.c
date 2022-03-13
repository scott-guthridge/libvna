/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * vnaconv_htoa: convert h-parameters to b-parameters
 */
void vnaconv_htob(const double complex (*h)[2], double complex (*b)[2])
{
    const double complex h11 = h[0][0];
    const double complex h12 = h[0][1];
    const double complex h21 = h[1][0];
    const double complex h22 = h[1][1];
#define B11	b[0][0]
#define B12	b[0][1]
#define B21	b[1][0]
#define B22	b[1][1]

    B11 =  1.0 / h12;
    B12 = -h11 / h12;
    B21 = -h22 / h12;
    B22 =  (h11 * h22 - h12 * h21) / h12;
}
#undef B22
#undef B21
#undef B12
#undef B11
