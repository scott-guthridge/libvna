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
 * vnaconv_gtou: convert g-parameters to u-parameters
 */
void vnaconv_gtou(const double complex (*g)[2], double complex (*u)[2],
	       const double complex *z0)
{
    const double complex g11 = g[0][0];
    const double complex g12 = g[0][1];
    const double complex g21 = g[1][0];
    const double complex g22 = g[1][1];
    const double complex z1  = z0[0];
    const double complex z2  = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = k2i / k1i * g12 * (z1 + z1c);
#define U11	u[0][0]
#define U12	u[0][1]
#define U21	u[1][0]
#define U22	u[1][1]

    U11 = (-g11 * z1  * z2  - g11 * g22 * z1  + g12 * g21 * z1  - g22 - z2)
	/ d;
    U12 = (-g11 * z1c * z2  - g11 * g22 * z1c + g12 * g21 * z1c + g22 + z2)
	/ d;
    U21 = (+g11 * z1  * z2c - g11 * g22 * z1  + g12 * g21 * z1  - g22 + z2c)
	/ d;
    U22 = (+g11 * z1c * z2c - g11 * g22 * z1c + g12 * g21 * z1c + g22 - z2c)
	/ d;
}
#undef U22
#undef U21
#undef U12
#undef U11
