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
 * vnaconv_g2zi: calculate the two-port input port impedances
 */
void vnaconv_g2zi(const vnaconv_array2_t *g, double complex *zi,
	const double complex *z0)
{
    const double complex g11 = g[0][0];
    const double complex g12 = g[0][1];
    const double complex g21 = g[1][0];
    const double complex g22 = g[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex g11_g22 = g11 * g22;
    const double complex g12_g21 = g12 * g21;

    zi[0] = (g22 + z2) / (g11_g22 - g12_g21 + g11 * z2);
    zi[1] = (g22 + (g11_g22 - g12_g21) * z1) / (1.0 + g11 * z1);
}
