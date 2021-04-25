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

#include "src/archdep.h"

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
#include "src/vnacommon_internal.h"
#include "test.h"


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
static bool opt_a = false;
static int opt_v = 0;

/*
 * test_qrd: test QR decomposition
 */
static test_result_t test_vnacommon_qrd()
{
    test_result_t result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int m = 1; m <= 5; ++m) {
	    for (int n = 1; n <= 5; ++n) {
		int diagonals = MIN(m, n);
		double complex a[m][n];
		double complex t[m][n];
		double complex d[diagonals];
		double complex q[m][m];
		double complex r[m][n];

		/*
		 * If -v, print the test header.
		 */
		if (opt_v) {
		    (void)printf("Test vnacommon_qrd: trial %3d "
			    "size %d x %d\n", trial, m, n);
		    (void)fflush(stdout);
		}

		/*
		 * Fill A with random numbers and copy to R.
		 * Decompose into factored Q_n and R components.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			r[i][j] = a[i][j] = test_crandn();
		    }
		}
		_vnacommon_qrd(*r, d, m, n);
		if (opt_v) {
		    test_print_cmatrix("a",  *a, m, n);
		    test_print_cmatrix("qr", *r, m, n);
		    test_print_cmatrix("d",  d, 1, diagonals);
		    (void)fflush(stdout);
		}

		/*
		 * Initialize q to the identity matrix.
		 */
		for (int i = 0; i < m; ++i) {
		    for (int j = 0; j < m; ++j) {
			q[i][j] = (i == j) ? 1.0 : 0.0;
		    }
		}

		/*
		 * Form Q.
		 */
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    for (int i = 0; i < m; ++i) {
			double complex s = 0.0;

			for (int j = diagonal; j < m; ++j) {
			    s += q[i][j] * r[j][diagonal];
			}
			for (int j = diagonal; j < m; ++j) {
			    q[i][j] -= 2.0 * s * conj(r[j][diagonal]);
			}
		    }
		}
		if (opt_v) {
		    test_print_cmatrix("q", *q, m, m);
		    (void)fflush(stdout);
		}

		/*
		 * Form R.
		 */
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    r[diagonal][diagonal] = d[diagonal];
		}
		for (int i = 0; i < diagonals; ++i) {
		    for (int j = 0; j < i; ++j) {
			r[i][j] = 0.0;
		    }
		}
		for (int i = n; i < m; ++i) {
		    for (int j = 0; j < n; ++j) {
			r[i][j] = 0.0;
		    }
		}
		if (opt_v) {
		    test_print_cmatrix("r", *r, m, n);
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
			if (!test_isequal(s, i == j ? 1.0 : 0.0)) {
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
			if (!test_isequal(t[i][j], a[i][j])) {
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
    test_report(result);;
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
	    opt_v = true;
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
    test_init_isequal();
    exit(test_vnacommon_qrd());
}
