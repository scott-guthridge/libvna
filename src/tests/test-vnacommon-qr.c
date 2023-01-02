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

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnacommon_internal.h"
#include "libt.h"


#define N_MATRIX_TRIALS	100

/*
 * Options
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
int opt_v = 0;

/*
 * test_qr: test QR decomposition
 */
static libt_result_t test_vnacommon_qr()
{
    libt_result_t result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int m = 1; m <= 5; ++m) {
	    for (int n = 1; n <= 5; ++n) {
		double complex a[m][n];
		double complex t[m][n];
		double complex q[m][m];
		double complex r[m][n];

		/*
		 * If -v, print the test header.
		 */
		if (opt_v) {
		    (void)printf("Test vnacommon_qr: trial %3d "
			    "size %d x %d\n", trial, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Fill A with random numbers and copy to T.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			t[i][j] = a[i][j] = libt_crandn();
		    }
		}


		/*
		 * Find the QR decomposition.
		 */
		_vnacommon_qr(*t, *q, *r, m, n);
		if (opt_v) {
		    libt_print_cmatrix("a", *a, m, n);
		    libt_print_cmatrix("q", *q, m, m);
		    libt_print_cmatrix("r", *r, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Test that Q Q' is the identity matrix.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < m; ++j) {
			double complex s = 0.0;

			for (int k = 0; k < m; ++k) {
			    s += q[i][k] * conj(q[j][k]);
			}
			if (!libt_isequal(s, i == j ? 1.0 : 0.0)) {
			    if (opt_a) {
				assert(!"data miscompare");
			    }
			    result = T_FAIL;
			    goto out;
			}
		    }
		}

		/*
		 * Test that R is upper-triangular.
		 */
		for (int i = 1; i < m; ++i) {
		    for (int j = 0; j < MIN(i, n); ++j) {
			if (!libt_isequal(r[i][j], 0.0)) {
			    if (opt_a) {
				assert(!"data miscompare");
			    }
			    result = T_FAIL;
			    goto out;
			}
		    }
		}

		/*
		 * Test that Q R == A
		 */
		_vnacommon_mmultiply(*t, *q, *r, m, m, n);
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			if (!libt_isequal(t[i][j], a[i][j])) {
			    if (opt_a) {
				assert(!"data miscompare");
			    }
			    result = T_FAIL;
			    goto out;
			}
		    }
		}
	    }
	}
    }
    result = T_PASS;

out:
    libt_report(result);;
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
    exit(2);
}

/*
 * main: test program
 */
int
main(int argc, char **argv)
{
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	case -1:
	    break;

	default:
	    print_usage();
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
    }
    libt_isequal_init();
    exit(test_vnacommon_qr());
}
