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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vnaconv.h>

/* reference impedances */
#define Z1	75.0
#define Z2	50.0

/* resistor values for impedance matching L pad */
#define R1	(sqrt(Z1) * sqrt(Z1 - Z2))
#define R2	(sqrt(Z1) * Z2 / sqrt(Z1 - Z2))

/* reference impedance vector */
static const double complex z0[] = { Z1, Z2 };

int main(int argc, char **argv)
{
    const double complex z[2][2] = { /* Z-parameters of the L pad */
	{ R1+R2, R2 },
	{ R2,    R2 }
    };
    double complex s[2][2];
    double complex zi[2];

    /*
     * Convert to S-parameters.
     */
    vnaconv_ztos(z, s, z0);
    (void)printf("s-parameters:\n");
    (void)printf("  %7.4f%+7.4fi    %7.4f%+7.4fi\n",
	creal(s[0][0]), cimag(s[0][0]), creal(s[0][1]), cimag(s[0][1]));
    (void)printf("  %7.4f%+7.4fi    %7.4f%+7.4fi\n",
	creal(s[1][0]), cimag(s[1][0]), creal(s[1][1]), cimag(s[1][1]));
    (void)printf("\n");

    /*
     * Convert to input impedance at each port.
     */
    vnaconv_stozi(s, zi, z0);
    (void)printf("input-impedances:\n");
    (void)printf("  %7.4f%+7.4fi    %7.4f%+7.4fi\n",
	creal(zi[0]), cimag(zi[0]), creal(zi[1]), cimag(zi[1]));
    (void)printf("\n");

    exit(0);
}
