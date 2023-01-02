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
 * vnaconv_stob: convert s-parameters to b-parameters
 */
void vnaconv_stob(const double complex (*s)[2], double complex (*b)[2],
	       const double complex *z0)
{
    const double complex s11 = s[0][0];
    const double complex s12 = s[0][1];
    const double complex s21 = s[1][0];
    const double complex s22 = s[1][1];
    const double complex z1 = z0[0];
    const double complex z2 = z0[1];
    const double complex z1c = conj(z1);
    const double complex z2c = conj(z2);
    const double k1i = sqrt(fabs(creal(z1)));
    const double k2i = sqrt(fabs(creal(z2)));
    const double complex db = s11 * s22 - s12 * s21;
    const double complex d = k1i / k2i * s12 * (z2 + z2c);
#define B11	b[0][0]
#define B12	b[0][1]
#define B21	b[1][0]
#define B22	b[1][1]

    B11 = -((db - s22) * z2 + (s11 - 1.0) * z2c) / d;
    B12 = -((db * z2 + s11 * z2c) * z1 + (s22 * z2 + z2c) * z1c) / d;
    B21 = -(1.0 + db - s11 - s22) / d;
    B22 = -((db - s11) * z1 + (s22 - 1.0) * z1c) / d;
}
#undef B22
#undef B21
#undef B12
#undef B11
