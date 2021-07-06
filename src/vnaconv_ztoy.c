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

#include <complex.h>
#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_ztoy: convert z-parameters to y-parameters
 */
void vnaconv_ztoy(const double complex (*z)[2], double complex (*y)[2])
{
    const double complex z11 = z[0][0];
    const double complex z12 = z[0][1];
    const double complex z21 = z[1][0];
    const double complex z22 = z[1][1];
    const double complex d = z11 * z22 - z12 * z21;
#define Y11	y[0][0]
#define Y12	y[0][1]
#define Y21	y[1][0]
#define Y22	y[1][1]

    Y11 =  z22 / d;
    Y12 = -z12 / d;
    Y21 = -z21 / d;
    Y22 =  z11 / d;
}
#undef Y22
#undef Y21
#undef Y12
#undef Y11
