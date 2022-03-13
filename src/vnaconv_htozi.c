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

#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_htozi: calculate the two-port input port impedances
 */
void vnaconv_htozi(const double complex (*h)[2], double complex *zi,
	const double complex *z0)
{
    const double complex h11 = h[0][0];
    const double complex h12 = h[0][1];
    const double complex h21 = h[1][0];
    const double complex h22 = h[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex h11_h22 = h11 * h22;
    const double complex h12_h21 = h12 * h21;

    zi[0] = (h11 + (h11_h22 - h12_h21) * z2) / (1.0 + h22 * z2);
    zi[1] = (h11 + z1) / (h11_h22 - h12_h21 + h22 * z1);
}
