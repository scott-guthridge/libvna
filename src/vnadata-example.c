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

#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vnadata.h>

#define PI      3.14159265
#define FMIN	100e+3		/* Hz */
#define FMAX	1e+9		/* Hz */
#define N	9		/* number of frequencies */
#define L	796e-9		/* Henries */
#define C	318e-12		/* Farads */

/*
 * error_fn: error printing function for the library
 *   @message: single line error message without a newline
 *   @error_arg: passed through to the error function (unused here)
 *   @category: category of error (ignored here)
 */
static void error_fn(const char *message, void *error_arg,
	vnaerr_category_t category)
{
    (void)fprintf(stderr, "example: %s\n", message);
}

/*
 * main
 */
int main(int argc, char **argv)
{
    vnadata_t *vdp;
    const double fstep = log(FMAX / FMIN) / (double)(N - 1);

    /*
     * Set up Z-parameter matrix for an L-C divider.
     */
    if ((vdp = vnadata_alloc_and_init(error_fn, /*error_arg*/NULL,
		    VPT_Z, 2, 2, N)) == NULL) {
	exit(1);
    }
    for (int findex = 0; findex < N; ++findex) {
	double f = FMIN * exp((double)findex * fstep);
	double complex s = 2 * PI * I * f;
	double complex z[2][2];

	if (vnadata_set_frequency(vdp, findex, f) == -1) {
	    exit(2);
	}
	z[0][0] = 1.0 / (C * s) + L * s;
	z[0][1] = 1.0 / (C * s);
	z[1][0] = z[0][1];
	z[1][1] = z[0][1];
	if (vnadata_set_matrix(vdp, findex, &z[0][0]) == -1) {
	    exit(3);
	}
    }

    /*
     * Save the parameters in Z real-imaginary, S dB, and Zin
     * magnitude-angle formats.
     */
    if (vnadata_set_format(vdp, "Zri,SdB,Zinma") == -1) {
	exit(4);
    }
    if (vnadata_save(vdp, "vnadata-example.npd") == -1) {
	exit(5);
    }

    /*
     * Print the Z parameters.
     */
    (void)printf("z-parameters (real-imaginary)\n");
    (void)printf("-------------------------\n");
    for (int findex = 0; findex < N; ++findex) {
	double f = vnadata_get_frequency(vdp, findex);

	(void)printf("f %7.2f MHz\n", f / 1.0e+6);
	for (int row = 0; row < 2; ++row) {
	    for (int column = 0; column < 2; ++column) {
		double complex value;

		value = vnadata_get_cell(vdp, findex, row, column);
		(void)printf("  %6.1f %6.1f%s",
			creal(value), cimag(value),
			column < 1 ? "," : "");
	    }
	    (void)printf("\n");
	}
	(void)printf("\n");
    }
    (void)printf("\n");


    /*
     * Convert to S-parameters and print.
     */
    if (vnadata_convert(vdp, vdp, VPT_S) == -1) {
	exit(6);
    }
    (void)printf("s-parameters (dB-degrees)\n");
    (void)printf("-------------------------\n");
    for (int findex = 0; findex < N; ++findex) {
	double f = vnadata_get_frequency(vdp, findex);

	(void)printf("f %7.2f MHz\n", f / 1.0e+6);
	for (int row = 0; row < 2; ++row) {
	    for (int column = 0; column < 2; ++column) {
		double complex value;

		value = vnadata_get_cell(vdp, findex, row, column);
		(void)printf("  %5.1f %6.1f%s",
			20 * log10(cabs(value)), 180 / PI * carg(value),
			column < 1 ? "," : "");
	    }
	    (void)printf("\n");
	}
	(void)printf("\n");
    }
    (void)printf("\n");

    /*
     * Convert to impedance into each port and print.
     */
    if (vnadata_convert(vdp, vdp, VPT_ZIN) == -1) {
	exit(7);
    }
    (void)printf("input-impedances (ohms-degrees)\n");
    (void)printf("------------------------------\n");
    for (int findex = 0; findex < N; ++findex) {
	double f = vnadata_get_frequency(vdp, findex);

	(void)printf("f %7.2f MHz\n", f / 1.0e+6);
	for (int port = 0; port < 2; ++port) {
	    double complex value;

	    value = vnadata_get_cell(vdp, findex, 0, port);
	    (void)printf("  %9.2f %6.1f%s",
		    cabs(value), 180 / PI * carg(value),
		    port < 1 ? "," : "");
	}
	(void)printf("\n");
    }
    (void)printf("\n");
    exit(0);
}
