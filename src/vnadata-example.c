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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vnadata.h>

/* system impedances */
#define Z1	75.0
#define Z2	50.0

/* values for impedance matching L pad */
#define R1	(sqrt(Z1) * sqrt(Z1 - Z2))
#define R2	(sqrt(Z1) * Z2 / sqrt(Z1 - Z2))

/* system impedance vector */
static const double complex z0_vector[] = { Z1, Z2 };

int main(int argc, char **argv)
{
    /* Z-parameters of the L-pad */
    const double complex z[2][2] = {
	{ R1+R2, R2 },
	{ R2,    R2 }
    };
    vnadata_t *vdp;

    /*
     * Set up Z-parameter matrix.
     */
    if ((vdp = vnadata_alloc_and_init(1, 2, 2, VPT_Z)) == NULL) {
	(void)fprintf(stderr, "%s: vnadata_alloc_and_init: %s\n",
		argv[0], strerror(errno));
	exit(1);
    }
    if (vnadata_set_frequency(vdp, 0, 1.0e+6) == -1) {
	(void)fprintf(stderr, "%s: vnadata_set_frequency: %s\n",
		argv[0], strerror(errno));
	exit(2);
    }
    if (vnadata_set_matrix(vdp, 0, &z[0][0]) == -1) {
	(void)fprintf(stderr, "%s: vnadata_set_matrix: %s\n",
		argv[0], strerror(errno));
	exit(3);
    }
    if (vnadata_set_z0_vector(vdp, z0_vector) == -1) {
	(void)fprintf(stderr, "%s: vnadata_set_z0_vector: %s\n",
		argv[0], strerror(errno));
	exit(4);
    }

    /*
     * Convert to S-parameters and print.
     */
    if (vnadata_convert(vdp, vdp, VPT_S) == -1) {
	(void)fprintf(stderr, "%s: vnadata_convert: %s\n",
		argv[0], strerror(errno));
	exit(5);
    }
    (void)printf("s-parameters:\n");
    for (int row = 0; row < 2; ++row) {
	for (int column = 0; column < 2; ++column) {
	    double complex value = vnadata_get_cell(vdp, 0, row, column);

	    (void)printf("  %7.4f%+7.4fj", creal(value), cimag(value));
	}
	(void)printf("\n");
    }
    (void)printf("\n");

    /*
     * Convert to input impedance at each port and print.
     */
    if (vnadata_convert(vdp, vdp, VPT_ZIN) == -1) {
	(void)fprintf(stderr, "%s: vnadata_convert: %s\n",
		argv[0], strerror(errno));
	exit(6);
    }
    (void)printf("input-impedances:\n");
    for (int port = 0; port < 2; ++port) {
	double complex value = vnadata_get_cell(vdp, 0, 0, port);

	(void)printf("  %7.4f%+7.4fj", creal(value), cimag(value));
    }
    (void)printf("\n");

    exit(0);
}
