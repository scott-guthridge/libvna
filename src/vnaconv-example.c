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

#include <stdlib.h>
#include <stdio.h>
#include <vnaconv.h>

static double complex s[2][2] = {
    {  0.000000 + 0.000000 * I,  0.000000 + 0.000000 * I },
    {  10.00000 + 0.000000 * I,  0.000000 + 0.000000 * I }
};

static double complex z0[2] = { 50.0, 50.0 };

int
main(int argc, char **argv)
{
    double complex z[2][2];

    vnaconv_s2z(s, z, z0);
    (void)printf("%7.1f%+7.1fi    %7.1f%+7.1fi\n",
        creal(z[0][0]), cimag(z[0][0]), creal(z[0][1]), cimag(z[0][1]));
    (void)printf("%7.1f%+7.1fi    %7.1f%+7.1fi\n",
        creal(z[1][0]), cimag(z[1][0]), creal(z[1][1]), cimag(z[1][1]));

    exit(0);
}

