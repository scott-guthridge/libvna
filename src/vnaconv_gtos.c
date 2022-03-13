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

#include <complex.h>
#include <math.h>
#include "vnaconv_internal.h"


/*
 * vnaconv_gtos: convert g-parameters to s-parameters
 */
void vnaconv_gtos(const double complex (*g)[2], double complex (*s)[2],
	       const double complex *z0)
{
    const double complex g11 = g[0][0];
    const double complex g12 = g[0][1];
    const double complex g21 = g[1][0];
    const double complex g22 = g[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex dg = g11 * g22 - g12 * g21;
    const double complex d = (dg + g11 * z2) * z1 + g22 + z2;
#define S11	s[0][0]
#define S12	s[0][1]
#define S21	s[1][0]
#define S22	s[1][1]

    S11 = -((dg + g11 * z2) * z1c - g22 - z2)  / d;
    S12 = -k2i / k1i * g12 * (z1 + z1c) / d;
    S21 =  k1i / k2i * g21 * (z2 + z2c) / d;
    S22 =  ((dg - g11 * z2c) * z1 + g22 - z2c) / d;
}
#undef S22
#undef S21
#undef S12
#undef S11
