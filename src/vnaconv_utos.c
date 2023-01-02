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
 * vnaconv_utos: convert t-parameters to s-parameters
 */
void vnaconv_utos(const double complex (*u)[2], double complex (*s)[2])
{
    const double complex u11 = u[0][0];
    const double complex u12 = u[0][1];
    const double complex u21 = u[1][0];
    const double complex u22 = u[1][1];
#define S11	s[0][0]
#define S12	s[0][1]
#define S21	s[1][0]
#define S22	s[1][1]

    S11 =       -u12 / u11;
    S12 =        1.0 / u11;
    S21 = -u12 * u21 / u11 + u22;
    S22 =        u21 / u11;
}
#undef S22
#undef S21
#undef S12
#undef S11
