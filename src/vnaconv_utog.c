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
 * vnaconv_utog: convert u-parameters to g-parameters
 */
void vnaconv_utog(const double complex (*u)[2], double complex (*g)[2],
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
    const double complex d = (u11 - u21) * z1c - (u12 - u22) * z1;
#define G11	g[0][0]
#define G12	g[0][1]
#define G21	g[1][0]
#define G22	g[1][1]

    G11 = (u11 + u12 - u21 - u22) / d;
    G12 = -k1i / k2i * (z2 + z2c) / d;
    G21 =  k2i / k1i * (u11 * u22 - u12 * u21) * (z1 + z1c) / d;
    G22 = (u11 * z1c * z2c - u12 * z1  * z2c
	 + u21 * z1c * z2  - u22 * z1  * z2)  / d;
}
#undef G22
#undef G21
#undef G12
#undef G11
