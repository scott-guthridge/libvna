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
 * test_vnacommon_lu: test LU factorization
 */
static libt_result_t test_vnacommon_lu()
{
    static const int sizes[] = { 1, 2, 3, 10 };
    double complex a[10 * 10];
    double complex t[10 * 10];
    int row_index[10];
    libt_result_t result = T_SKIPPED;

#define A(i, j) (a[(i) * n + (j)])
#define T(i, j) (t[(i) * n + (j)])
    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int size_index = 0; size_index < sizeof(sizes) / sizeof(int);
		++size_index) {
	    int n = sizes[size_index];
	    double complex d;

	    /*
	     * If -v, print the test header.
	     */
	    if (opt_v) {
		(void)printf("Test vnacommon_lu: trial %3d size %d x %d\n",
			trial, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Generate a random complex matrix in T and copy to A.
	     */
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
		    A(i, j) = T(i, j) = libt_crandn();
		}
	    }
	    if (opt_v) {
		libt_print_cmatrix("a", a, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Compute the in-place L U factorization.
	     */
	    errno = 0;
	    d = _vnacommon_lu(a, row_index, n);
	    if (opt_v) {
		libt_print_cmatrix("LU factorization", a, n, n);
		(void)printf("determinant %8.5f%+8.5fj\n", creal(d), cimag(d));
		(void)printf("\n");
		(void)fflush(stdout);
	    }
	    if (cabs(d) < libt_isequal_eps) {
		(void)fprintf(stderr, "%s: test_vnacommon_lu: warning: "
			"skipping nearly singular test matrix\n", progname);
		continue;
	    }

	    /*
	     * Check the result.
	     */
	    for (int i = 0; i < n; ++i) {
		for (int j = 0; j < n; ++j) {
		    double complex s = (i <= j) ? A(i, j) : 0.0;

		    for (int k = 0; k < i; ++k) {
			if (k <= j) {
			    s += A(i, k) * A(k, j);
			}
		    }
		    if (!libt_isequal(s, T(row_index[i], j))) {
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
    result = T_PASS;

out:
    libt_report(result);;
    return result;
}
#undef A
#undef T

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
    exit(test_vnacommon_lu());
}
