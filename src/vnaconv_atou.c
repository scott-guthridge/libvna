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
 * vnaconv_atou: convert a-parameters to u-parameters
 */
void vnaconv_atou(const double complex (*a)[2], double complex (*u)[2],
	       const double complex *z0)
{
    const double complex a11 = a[0][0];
    const double complex a12 = a[0][1];
    const double complex a21 = a[1][0];
    const double complex a22 = a[1][1];
    const double complex z1  = z0[0];
    const double complex z2  = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = k2i / k1i * (a11 * a22 - a12 * a21) * (z1 + z1c);
#define U11	u[0][0]
#define U12	u[0][1]
#define U21	u[1][0]
#define U22	u[1][1]

    U11 =  (a11 * z2  + a12 + a21 * z1  * z2  + a22 * z1)  / d;
    U12 = (-a11 * z2  - a12 + a21 * z1c * z2  + a22 * z1c) / d;
    U21 = (-a11 * z2c + a12 - a21 * z1  * z2c + a22 * z1)  / d;
    U22 =  (a11 * z2c - a12 - a21 * z1c * z2c + a22 * z1c) / d;
}
#undef U22
#undef U21
#undef U12
#undef U11
