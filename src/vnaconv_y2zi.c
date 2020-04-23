/*
 * Electrical Network Parameter Conversion Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * vnaconv_y2zi: calculate the two-port input port impedances
 */
void vnaconv_y2zi(const vnaconv_array2_t *y, double complex *zi,
	const double complex *z0)
{
    const double complex y11 = y[0][0];
    const double complex y12 = y[0][1];
    const double complex y21 = y[1][0];
    const double complex y22 = y[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];

    zi[0] = (1.0 + y22 * z2) / (y11 * (1.0 + y22 * z2) - y12 * y21 * z2);
    zi[1] = (1.0 + y11 * z1) / (y22 * (1.0 + y11 * z1) - y12 * y21 * z1);
}
