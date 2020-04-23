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
 * vnaconv_s2zi: calculate the two-port input port impedances
 */
void vnaconv_s2zi(const vnaconv_array2_t *s, double complex *zi,
	const double complex *z0)
{
    const double complex s11 = s[0][0];
    const double complex s22 = s[1][1];
    const double complex z1  = z0[0];
    const double complex z2  = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);

    zi[0] = (s11 * z1 + z1c) / (1.0 - s11);
    zi[1] = (s22 * z2 + z2c) / (1.0 - s22);
}
