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
 * vnaconv_ytou: convert y-parameters to u-parameters
 */
void vnaconv_ytou(const double complex (*y)[2], double complex (*u)[2],
	       const double complex *z0)
{
    const double complex y11 = y[0][0];
    const double complex y12 = y[0][1];
    const double complex y21 = y[1][0];
    const double complex y22 = y[1][1];
    const double complex z1  = z0[0];
    const double complex z2  = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = k2i / k1i * y12 * (z1 + z1c);
#define U11	u[0][0]
#define U12	u[0][1]
#define U21	u[1][0]
#define U22	u[1][1]

    U11 = (-1.0 - y11 * z1  - y11 * y22 * z1  * z2
		+ y12 * y21 * z1  * z2  - y22 * z2)  / d;
    U12 =  (1.0 - y11 * z1c - y11 * y22 * z1c * z2
		+ y12 * y21 * z1c * z2  + y22 * z2)  / d;
    U21 = (-1.0 - y11 * z1  + y11 * y22 * z1  * z2c
		- y12 * y21 * z1  * z2c + y22 * z2c) / d;
    U22 =  (1.0 - y11 * z1c + y11 * y22 * z1c * z2c
		- y12 * y21 * z1c * z2c - y22 * z2c) / d;
}
#undef U22
#undef U21
#undef U12
#undef U11
