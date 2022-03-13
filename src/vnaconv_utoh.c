
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

#include <complex.h>
#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_utoh: convert t-parameters to h-parameters
 */
void vnaconv_utoh(const double complex (*u)[2], double complex (*h)[2],
	       const double complex *z0)
{
    const double complex u11 = u[0][0];
    const double complex u12 = u[0][1];
    const double complex u21 = u[1][0];
    const double complex u22 = u[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = (u11 + u12) * z2c + (u21 + u22) * z2;

#define H11	h[0][0]
#define H12	h[0][1]
#define H21	h[1][0]
#define H22	h[1][1]

    H11 = (u11 * z1c * z2c - u12 * z1  * z2c
         + u21 * z1c * z2 - u22 * z1  * z2) / d;
    H12 =  k1i / k2i * (z2 + z2c) / d;
    H21 = -k2i / k1i * (u11 * u22 - u12 * u21) * (z1 + z1c) / d;
    H22 = (u11 + u12 - u21 - u22) / d;
}
#undef H22
#undef H21
#undef H12
#undef H11
