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

#include <complex.h>
#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_ztob: convert z-parameters to b-parameters
 */
void vnaconv_ztob(const double complex (*z)[2], double complex (*b)[2])
{
    const double complex z11 = z[0][0];
    const double complex z12 = z[0][1];
    const double complex z21 = z[1][0];
    const double complex z22 = z[1][1];
#define B11	b[0][0]
#define B12	b[0][1]
#define B21	b[1][0]
#define B22	b[1][1]

    B11 =  z22 / z12;
    B12 =  (z12 * z21 - z11 * z22) / z12;
    B21 = -1.0 / z12;
    B22 =  z11 / z12;
}
#undef B22
#undef B21
#undef B12
#undef B11
