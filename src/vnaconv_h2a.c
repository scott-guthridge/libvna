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
 * vnaconv_h2a: convert h-parameters to a-parameters
 */
void vnaconv_h2a(const vnaconv_array2_t *h, double complex (*a)[2])
{
    const double complex h11 = h[0][0];
    const double complex h12 = h[0][1];
    const double complex h21 = h[1][0];
    const double complex h22 = h[1][1];
#define A11	a[0][0]
#define A12	a[0][1]
#define A21	a[1][0]
#define A22	a[1][1]

    A11 =  (h12 * h21 - h11 * h22) / h21;
    A12 = -h11 / h21;
    A21 = -h22 / h21;
    A22 = -1.0 / h21;
}
#undef A22
#undef A21
#undef A12
#undef A11
