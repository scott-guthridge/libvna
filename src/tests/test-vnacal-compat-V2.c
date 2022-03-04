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

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnacal_internal.h"
#include "libt.h"
#include "libt_vnacal.h"


#define NTRIALS		67

/*
 * Command Line Options
 */
char *progname;
static const char options[] = "av";
static const char *const usage[] = {
    "[-av]",
    NULL
};
static const char *const help[] = {
    "-a	 abort on data miscompare",
    "-v	 show verbose output",
    NULL
};
bool opt_a = false;
int  opt_v = 0;

/*
 * filename: test file name
 */
static const char file[] = "compat-V2.vnacal";
static const char *pathname = file;

/*
 * error_fn: error reporting function
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)fprintf(stderr, "%s: %s\n", progname, message);
}

/*
 * CV2_F: number of freuquencies for the compat_V2 test
 */
#define CV2_F	11

/*
 * compat_V2_frequency_vector: frequency vector for compat_V2 test
 */
static const double compat_V2_frequency_vector[CV2_F] = {
    1.000000e+05,
    1.584893e+05,
    2.511886e+05,
    3.981072e+05,
    6.309573e+05,
    1.000000e+06,
    1.584893e+06,
    2.511886e+06,
    3.981072e+06,
    6.309573e+06,
    1.000000e+07
};

/*
 * compat_V2_measured: "measured" s-parameters for old "VNACAL 2.0" format
 *
 * These tables were generated using the E12 example with 11 calibration
 * and 11 measurement points from 100 kHz to 10 MHz.
 */
static const double complex compat_V2_measured[4][CV2_F] = {
    {	/* s11 */
	-3.926540e-03 + 4.341532e-04*I,
	-9.773344e-03 + 1.726204e-03*I,
	-2.397616e-02 + 6.848069e-03*I,
	-5.649961e-02 + 2.696819e-02*I,
	-1.171773e-01 + 1.031452e-01*I,
	-1.379310e-01 + 3.448276e-01*I,
	+2.440045e-01 + 6.724525e-01*I,
	+8.548239e-01 + 3.846625e-01*I,
	+9.586034e-01 - 2.523617e-01*I,
	+6.399219e-01 - 7.672778e-01*I,
	+1.061320e-01 - 9.942893e-01*I
    },
    {	/* s12 */
	+9.939136e-01 - 1.099960e-01*I,
	+9.845757e-01 - 1.742980e-01*I,
	+9.604276e-01 - 2.759055e-01*I,
	+8.958771e-01 - 4.339283e-01*I,
	+7.175210e-01 - 6.559979e-01*I,
	+2.891602e-01 - 8.052561e-01*I,
	-1.570320e-01 - 5.267873e-01*I,
	-1.809236e-01 - 1.774419e-01*I,
	-9.240888e-02 - 4.711767e-02*I,
	-3.972649e-02 - 1.194356e-02*I,
	-1.625128e-02 - 3.004438e-03*I
    },
    {	/* s21 */
	+9.939350e-01 - 1.098983e-01*I,
	+9.847092e-01 - 1.739230e-01*I,
	+9.612490e-01 - 2.745518e-01*I,
	+9.006954e-01 - 4.299166e-01*I,
	+7.414183e-01 - 6.526327e-01*I,
	+3.448276e-01 - 8.620690e-01*I,
	-2.383455e-01 - 6.568568e-01*I,
	-3.176208e-01 - 1.429263e-01*I,
	-1.275371e-01 + 3.357539e-02*I,
	-2.705835e-02 + 3.244345e-02*I,
	-1.185832e-03 + 1.110938e-02*I
    },
    {	/* s22 */
	+6.013177e-03 - 6.654756e-04*I,
	+1.496251e-02 - 2.648791e-03*I,
	+3.666232e-02 - 1.053212e-02*I,
	+8.590211e-02 - 4.160766e-02*I,
	+1.728184e-01 - 1.580003e-01*I,
	+1.749419e-01 - 4.871800e-01*I,
	-2.386401e-01 - 8.005541e-01*I,
	-6.906384e-01 - 6.773475e-01*I,
	-8.860722e-01 - 4.517927e-01*I,
	-9.568318e-01 - 2.876665e-01*I,
	-9.832025e-01 - 1.817685e-01*I
    }
};

/*
 * compat_V2_m: measurement matrix for vnacal_apply_m
 */
static const double complex *compat_V2_m[4] = {
    compat_V2_measured[0], compat_V2_measured[1],
    compat_V2_measured[2], compat_V2_measured[3],
};

/*
 * compat_V2_expected: expected s-parameters for old "VNACAL 2.0" format
 */
static const double complex compat_V2_expected[4][CV2_F] = {
    {	/* s11 */
	-4.974876e-03 + 4.999875e-04*I,
	-1.239974e-02 + 1.990222e-03*I,
	-3.052222e-02 + 7.916587e-03*I,
	-7.250960e-02 + 3.135099e-02*I,
	-1.533550e-01 + 1.208076e-01*I,
	-2.000000e-01 + 4.000000e-01*I,
	+1.247191e-01 + 7.723058e-01*I,
	+6.206602e-01 + 7.235185e-01*I,
	+8.601119e-01 + 4.945027e-01*I,
	+9.473713e-01 + 3.161807e-01*I,
	+9.796082e-01 + 1.999200e-01*I
    },
    {	/* s12 */
	+9.949751e-01 - 9.999750e-02*I,
	+9.872848e-01 - 1.584643e-01*I,
	+9.674892e-01 - 2.509389e-01*I,
	+9.150093e-01 - 3.956228e-01*I,
	+7.704206e-01 - 6.069102e-01*I,
	+4.000000e-01 - 8.000000e-01*I,
	-9.930313e-02 - 6.149210e-01*I,
	-1.967360e-01 - 2.293399e-01*I,
	-1.085388e-01 - 6.240202e-02*I,
	-4.759378e-02 - 1.588420e-02*I,
	-1.959216e-02 - 3.998401e-03*I
    },
    {	/* s21 */
	+9.949751e-01 - 9.999750e-02*I,
	+9.872848e-01 - 1.584643e-01*I,
	+9.674892e-01 - 2.509389e-01*I,
	+9.150093e-01 - 3.956228e-01*I,
	+7.704206e-01 - 6.069102e-01*I,
	+4.000000e-01 - 8.000000e-01*I,
	-9.930313e-02 - 6.149210e-01*I,
	-1.967360e-01 - 2.293399e-01*I,
	-1.085388e-01 - 6.240202e-02*I,
	-4.759378e-02 - 1.588420e-02*I,
	-1.959216e-02 - 3.998401e-03*I
    },
    {	/* s22 */
	+4.974876e-03 - 4.999875e-04*I,
	+1.239974e-02 - 1.990222e-03*I,
	+3.052222e-02 - 7.916587e-03*I,
	+7.250960e-02 - 3.135099e-02*I,
	+1.533550e-01 - 1.208076e-01*I,
	+2.000000e-01 - 4.000000e-01*I,
	-1.247191e-01 - 7.723058e-01*I,
	-6.206602e-01 - 7.235185e-01*I,
	-8.601119e-01 - 4.945027e-01*I,
	-9.473713e-01 - 3.161807e-01*I,
	-9.796082e-01 - 1.999200e-01*I
    }
};

/*
 * test_vnacal_compat_e2: test compatibility load of old E term format
 */
static libt_result_t test_vnacal_compat_V2()
{
    vnacal_t *vcp = NULL;
    vnadata_t *vdp = NULL;
    libt_result_t result = T_FAIL;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_load VNACAL 2.0 format\n");
    }

    /*
     * Load the old format save file.
     */
    if ((vcp = vnacal_load(pathname, error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_load: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Create a vnadata_t structure to hold the result.
     */
    if ((vdp = vnadata_alloc(error_fn, NULL)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Apply the calibration.
     */
    if (vnacal_apply_m(vcp, 0, compat_V2_frequency_vector, CV2_F,
		(double complex *const *)compat_V2_m, 2, 2, vdp) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Check the result.
     */
    for (int findex = 0; findex < CV2_F; ++findex) {
	double f = compat_V2_frequency_vector[findex];

	if (opt_v >= 2) {
	    (void)printf("findex %d  f %e\n", findex, f);
	    (void)printf("  computed s parameters:\n");
	    for (int s_row = 0; s_row < 2; ++s_row) {
		(void)printf("  ");
		for (int s_column = 0; s_column < 2; ++s_column) {
		    double complex v;

		    v = vnadata_get_cell(vdp, findex, s_row, s_column);
		    (void)printf(" %8.5f%+8.5fj", creal(v), cimag(v));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	}
	for (int s_row = 0; s_row < 2; ++s_row) {
	    for (int s_column = 0; s_column < 2; ++s_column) {
		int s_cell = s_row * 2 + s_column;
		double complex expected, actual;

		expected = compat_V2_expected[s_cell][findex];
		actual = vnadata_get_cell(vdp, findex, s_row, s_column);
		if (!libt_isequal(actual, expected)) {
		    if (opt_a) {
			assert(!"data miscompare");
		    }
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    result = T_PASS;

out:
    vnadata_free(vdp);
    vnacal_free(vcp);
    libt_report(result);
    return result;
}

/*
 * print_usage: print a usage message and exit
 */
static void print_usage()
{
    const char *const *cpp;

    for (cpp = usage; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s: usage %s\n", progname, *cpp);
    }
    for (cpp = help; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s\n", *cpp);
    }
    exit(99);
}

/*
 * main
 */
int
main(int argc, char **argv)
{
    char *srcdir = NULL;

    /*
     * Parse Options
     */
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case -1:
	    break;

	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	default:
	    print_usage();
	}
	break;
    }

    /*
     * If srcdir is defined in the environment, incorporate it into
     * pathname.
     */
    if ((srcdir = getenv("srcdir")) != NULL) {
	char *cp;

	if ((cp = malloc(strlen(srcdir) + 1 +
			sizeof(file))) == NULL) {
	    (void)fprintf(stderr, "%s: malloc: %s\n",
		    progname, strerror(errno));
	    exit(99);
	}
	pathname = cp;
	(void)strcpy(cp, srcdir);
	cp += strlen(cp);
	*cp++ = '/';
	(void)strcpy(cp, file);
    }

    /*
     * Set the compare precision and run the test.
     */
    libt_isequal_init();
    if (libt_isequal_eps < 0.00001) {	/* file has 6 digits precision */
	libt_isequal_eps = 0.00001;
    }
    exit(test_vnacal_compat_V2());
}
