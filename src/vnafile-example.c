/*
 * Electrical Network Parameter Formatting Library
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

#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vnafile.h>

char *progname;

#define R	50.0
#define L	7.95775e-6
#define C	3.1831e-9
#define Z1	50.0
#if 0
#define Z2	75.0
#define Z3	110.0
#else
#define Z2	50.0
#define Z3	50.0
#endif

#define PI	3.14159265358979323846264338327950288419716939937508

#define NROWS		3
#define NCOLUMNS	3
#define NFREQUENCIES	7

static double frequency_vector[] = {
    1.0e+3, 1.0e+4, 1.0e+5, 1.0e+6, 1.0e+7, 1.0e+8, 1.0e+9
};
static double complex z0_vector[] = { Z1, Z2, Z3 };

static void calc_s3x3(vnadata_t *vdp, double r, double l, double c,
	const double complex *z0_vector, int findex, double frequency)
{
    double complex s = 2.0 * PI * I * frequency;
    double complex s11, s12, s13;
    double complex s21, s22, s23;
    double complex s31, s32, s33;
    double complex z1  = z0_vector[0];
    double complex z2  = z0_vector[1];
    double complex z3  = z0_vector[2];
    double complex z1c = conj(z1);
    double complex z2c = conj(z2);
    double complex z3c = conj(z3);

    s11 = (r - z1c + z3 + s*(l + ((r - z1c)*(z2 + z3) + z2*z3)*c +
	    l*c*s*(r - z1c + z2))) /
	  (r + z1 + z3 + s*(l + ((r + z1)*z2 + (r + z1 + z2)*z3)*c +
	    l*c*s*(r + z1 + z2)));
    
    s12 = c*z1*s*(z3 + l*s) /
          (r + z1 + z3 + s*(l + (r + z1)*c*z3 + l*c*s*(r + z1)));

    s13 = (z1 + c*z1*z2*s) /
          (r + z1 + s*(l + (r + z1)*c*z2 + l*c*s*(r + z1 + z2)));

    s21 = c*z2*s*(z3 + l*s) /
          (r + z3 + s*(l + r*c*z2 + (r + z2)*c*z3 + l*c*s*(r + z2)));

    s22 = (r + z1 + z3 + s*(l - (r + z1)*c*z2c + (r + z1 - z2c)*c*z3 +
    		l*c*s*(r + z1 - z2c))) /
	  (r + z1 + z3 + s*(l + (r + z1)*c*z2 + (r + z1 + z2)*c*z3 +
	  	l*c*s*(r + z1 + z2)));

    s23 = (r + z1)*c*z2*s /
          (r + z1 + s*(l + (r + z1)*c*z2 + l*c*s*(r + z1 + z2)));
    
    s31 = (z3 + c*z2*z3*s) /
          (r + z3 + s*(l + r*c*z2 + (r + z2)*c*z3 + l*c*s*(r + z2)));

    s32 = (r + z1)*c*z3*s /
	  (r + z1 + z3 + s*(l + (r + z1)*c*z3 + l*c*s*(r + z1)));

    s33 = (r + z1 - z3c + s*(l + (r + z1)*c*z2 - (r + z1 + z2)*c*z3c +
    		l*c*s*(r + z1 + z2))) /
	  (r + z1 + z3 + s*(l + (r + z1)*c*z2 + (r + z1 + z2)*c*z3 +
	  	l*c*s*(r + z1 + z2)));

    vnadata_set_cell(vdp, findex, 0, 0, s11);
    vnadata_set_cell(vdp, findex, 0, 1, s12);
    vnadata_set_cell(vdp, findex, 0, 2, s13);
    vnadata_set_cell(vdp, findex, 1, 0, s21);
    vnadata_set_cell(vdp, findex, 1, 1, s22);
    vnadata_set_cell(vdp, findex, 1, 2, s23);
    vnadata_set_cell(vdp, findex, 2, 0, s31);
    vnadata_set_cell(vdp, findex, 2, 1, s32);
    vnadata_set_cell(vdp, findex, 2, 2, s33);
}

/*
 * error_fn: error printing function for the library
 *   @message: single line error message without a newline
 *   @error_arg: passed through to the error function (unused here)
 */
static void error_fn(const char *message, void *error_arg)
{
    (void)fprintf(stderr, "example: %s\n", message);
}

int
main(int argc, char **argv)
{
    vnadata_t *vdp;
    vnafile_t *vfp;

    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    --argc;
    ++argv;

    vdp = vnadata_alloc_and_init(NFREQUENCIES, NROWS, NCOLUMNS, VPT_S);
    vnadata_set_frequency_vector(vdp, frequency_vector);
    for (int findex = 0; findex < NFREQUENCIES; ++findex) {
	calc_s3x3(vdp, R, L, C, z0_vector, findex, frequency_vector[findex]);
    }
    vnadata_set_z0_vector(vdp, z0_vector);

    vfp = vnafile_alloc(error_fn, NULL);
    if (vfp == NULL) {
	(void)fprintf(stderr, "%s: vnafile_create: %s\n",
	    progname, strerror(errno));
	exit(2);
	/*NOTREACHED*/
    }
    vnafile_set_format(vfp, "Sri");
    if (vnafile_save(vfp, "matrix.out", vdp) == -1) {
	(void)fprintf(stderr, "%s: vnafile_save: %s\n",
	    progname, strerror(errno));
	exit(3);
	/*NOTREACHED*/
    }

    exit(0);
    /*NOTREACHED*/
}
