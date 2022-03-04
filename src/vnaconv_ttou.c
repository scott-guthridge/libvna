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
 * vnaconv_ttou: convert t-parameters to u-parameters
 */
void vnaconv_ttou(const double complex (*t)[2], double complex (*u)[2])
{
    const double complex t11 = t[0][0];
    const double complex t12 = t[0][1];
    const double complex t21 = t[1][0];
    const double complex t22 = t[1][1];
    const double complex d = t11 * t22 - t12 * t21;
#define U11	u[0][0]
#define U12	u[0][1]
#define U21	u[1][0]
#define U22	u[1][1]

    U11 =  t22 / d;
    U12 = -t12 / d;
    U21 = -t21 / d;
    U22 =  t11 / d;
}
#undef U22
#undef U21
#undef U12
#undef U11
