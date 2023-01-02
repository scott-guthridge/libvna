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
 * vnaconv_btozi: calculate the two-port input port impedances
 */
void vnaconv_btozi(const double complex (*b)[2], double complex *zi,
	const double complex *z0)
{
    const double complex b11 = b[0][0];
    const double complex b12 = b[0][1];
    const double complex b21 = b[1][0];
    const double complex b22 = b[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];

    zi[0] = (-b12 + b22 * z2) / (b11 - b21 * z2);
    zi[1] =  (b12 - b11 * z1) / (b21 * z1 - b22);
}
