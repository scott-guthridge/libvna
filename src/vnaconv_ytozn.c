/*
 * Vector Network Analyzer Library
 * Copyright © 2020, 2021 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include <string.h>
#include "vnacommon_internal.h"
#include "vnaconv_internal.h"


/*
 * vnaconv_ytoz: convert y-parameters to z-parameters
 *   @y:  serialized y matrix in  (n x n)
 *   @z:  serialized z matrix out (n x n)
 *   @n:  length
 */
void vnaconv_ytozn(const double complex *y, double complex *z, int n)
{
    double complex u[n * n];

    (void)memcpy((void *)u, (void *)y, n * n * sizeof(double complex));
    _vnacommon_minverse(z, u, n);
}
