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
 * vnaconv_s2a: convert s-parameters to a-parameters
 */
void vnaconv_s2a(const vnaconv_array2_t *s, complex (*a)[2],
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
    const double complex da = s11 * s22 - s12 * s21;
    const double complex d = k2i / k1i * s21 * (z1 + z1c);
#define A11	a[0][0]
#define A12	a[0][1]
#define A21	a[1][0]
#define A22	a[1][1]

    A11 = -((da - s11) * z1 + (s22 - 1.0) * z1c) / d;
    A12 =  ((da * z2 + s11 * z2c) * z1 + (s22 * z2 + z2c) * z1c) / d;
    A21 =  (1.0 + da - s11 - s22) / d;
    A22 = -((da - s22) * z2 + (s11 - 1.0) * z2c) / d;
}
#undef A22
#undef A21
#undef A12
#undef A11
