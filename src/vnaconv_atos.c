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
 * vnaconv_atos: convert a-parameters to s-parameters
 */
void vnaconv_atos(const double complex (*a)[2], double complex (*s)[2],
	       const double complex *z0)
{
    const double complex a11 = a[0][0];
    const double complex a12 = a[0][1];
    const double complex a21 = a[1][0];
    const double complex a22 = a[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex d = a11 * z2 + a12 + a21 * z1 * z2 + a22 * z1;
#define S11	s[0][0]
#define S12	s[0][1]
#define S21	s[1][0]
#define S22	s[1][1]

    S11 =  (a11 * z2 + a12 - a21 * z2 * z1c - a22 * z1c)    / d;
    S12 =  k2i / k1i * (a11 * a22 - a12 * a21) * (z1 + z1c) / d;
    S21 =  k1i / k2i * (z2 + z2c)                           / d;
    S22 = -(a11 * z2c - a12 + a21 * z1 * z2c - a22 * z1)    / d;
}
#undef S22
#undef S21
#undef S12
#undef S11
