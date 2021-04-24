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
 * vnaconv_btoz: convert b-parameters to z-parameters
 */
void vnaconv_btoz(const vnaconv_array2_t *b, double complex (*z)[2])
{
    const double complex b11 = b[0][0];
    const double complex b12 = b[0][1];
    const double complex b21 = b[1][0];
    const double complex b22 = b[1][1];
#define Z11	z[0][0]
#define Z12	z[0][1]
#define Z21	z[1][0]
#define Z22	z[1][1]

    Z11 = -b22 / b21;
    Z12 = -1.0 / b21;
    Z21 = -(b11 * b22 - b12 * b21) / b21;
    Z22 = -b11 / b21;
}
#undef Z22
#undef Z21
#undef Z12
#undef Z11
