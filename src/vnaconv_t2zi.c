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
 * MERCHANTABILITY or FITNESS FOR A11 PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_t2zi: calculate the two-port input port impedances
 */
void vnaconv_t2zi(const vnaconv_array2_t *t, double complex *zi,
	const double complex *z0)
{
    const double complex t12 = t[0][1];
    const double complex t21 = t[1][0];
    const double complex t22 = t[1][1];
    const double complex z1  = z0[0];
    const double complex z2  = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);

    zi[0] = (-t12 * z1 - t22 * z1c) / (t12 - t22);
    zi[1] = (-t21 * z2 + t22 * z2c) / (t21 + t22);
}
