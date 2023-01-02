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
 * qrsolve_helper: generate a random system of equations and solve
 */
static int qrsolve_helper(double complex *x, double complex *a,
	double complex *b, int m, int n, int o)
{
    double complex u[m][n];
    double complex v[m][o];

    /*
     * Generate random matrices A and B, and make copies
     * in U and V, respectively.
     */
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    a[i * n + j] = u[i][j] = libt_crandn();
	}
	for (int k = 0; k < o; ++k) {
	    b[i * o + k] = v[i][k] = libt_crandn();
	}
    }

    /*
     * Solve the system.  This call destroys both u and v.
     */
    return _vnacommon_qrsolve(x, &u[0][0], &v[0][0], m, n, o);
}

/*
 * find_axb_error: find the squared error in A * X = B
 */
static double find_axb_error(const double complex *a, const double complex *x,
	const double complex *b, int m, int n, int o)
{
    double squared_error = 0.0;

    for (int k = 0; k < o; ++k) {
	for (int i = 0; i < m; ++i) {
	    double complex s = 0.0;
	    double e;

	    for (int j = 0; j < n; ++j) {
		s += a[i * n + j] * x[j * o + k];
	    }
	    e = cabs(s - b[i * o + k]);
	    squared_error += e * e;
	}
    }
    return squared_error;
}

/*
 * Test vnacommon_qrsolve.
 */
static libt_result_t test_vnacommon_qrsolve()
{
    libt_result_t result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	/*
	 * Test square coefficient matrices.
	 */
	for (int n = 1; n <= 10; ++n) {
	    const int o = 3;
	    double complex a[n][n];
	    double complex b[n][o];
	    double complex x[n][o];
	    int rank;

	    /*
	     * If -v, print the test header.
	     */
	    if (opt_v) {
		(void)printf("Test vnacommon_qrsolve: trial %3d "
			"size %d x %d\n", trial, n, n);
		(void)fflush(stdout);
	    }

	    /*
	     * Generate random matrices A and B, and solve for X.
	     */
	    rank = qrsolve_helper(&x[0][0], &a[0][0], &b[0][0], n, n, o);
	    if (opt_v) {
		libt_print_cmatrix("a", *a, n, n);
		libt_print_cmatrix("b", *b, n, o);
		libt_print_cmatrix("x", *x, n, o);
		(void)printf("rank %d\n", rank);
		(void)fflush(stdout);
	    }

	    /*
	     * Verify A X == B.
	     */
	    for (int k = 0; k < o; ++k) {
		for (int i = 0; i < n; ++i) {
		    double complex s = 0.0;

		    for (int j = 0; j < n; ++j) {
			s += a[i][j] * x[j][k];
		    }
		    if (!libt_isequal(s, b[i][k])) {
			if (opt_a) {
			    assert(!"data miscompare");
			}
			result = T_FAIL;
			goto out;
		    }
		}
	    }

	    /*
	     * Check rank.
	     */
	    if (rank != n) {
		if (opt_a) {
		    assert(!"incorrect rank");
		}
		result = T_FAIL;
		goto out;
	    }
	}

	/*
	 * Test more columns than rows (underdetermined case).
	 */
	for (int m = 1; m <= 4; ++m) {
	    for (int n = m + 1; n <= 5; ++n) {
		for (int o = 1; o <= 2; ++o) {
		    double complex a[m][n];
		    double complex b[m][o];
		    double complex x[n][o];
		    int rank;

		    /*
		     * If -v, print the test header.
		     */
		    if (opt_v) {
			(void)printf("Test vnacommon_qrsolve: trial %3d "
				"A size %d x %d, B size %d x %d\n",
				trial, m, n, n, o);
			(void)fflush(stdout);
		    }

		    /*
		     * Generate random matrices A and B, and solve for X.
		     */
		    rank = qrsolve_helper(&x[0][0], &a[0][0], &b[0][0],
			    m, n, o);
		    if (opt_v) {
			libt_print_cmatrix("a", *a, m, n);
			libt_print_cmatrix("b", *b, n, o);
			libt_print_cmatrix("x", *x, n, o);
			(void)printf("rank %d\n", rank);
			(void)fflush(stdout);
		    }

		    /*
		     * Verify A X == B.
		     */
		    for (int k = 0; k < o; ++k) {
			for (int i = 0; i < m; ++i) {
			    double complex s = 0.0;

			    for (int j = 0; j < n; ++j) {
				s += a[i][j] * x[j][k];
			    }
			    if (!libt_isequal(s, b[i][k])) {
				if (opt_a) {
				    assert(!"data miscompare");
				}
				result = T_FAIL;
				goto out;
			    }
			}
		    }

		    /*
		     * Check rank.
		     */
		    if (rank != m) {
			if (opt_a) {
			    assert(!"incorrect rank");
			}
			result = T_FAIL;
			goto out;
		    }
		}
	    }
	}

	/*
	 * Test more rows than columns (overdetermined case).
	 */
	for (int n = 1; n <= 4; ++n) {
	    for (int m = n + 1; m <= 5; ++m) {
		for (int o = 1; o <= 2; ++o) {
		    double complex a[m][n];
		    double complex b[m][o];
		    double complex x[n][o];
		    int rank;
		    double error0;

		    /*
		     * If -v, print the test header.
		     */
		    if (opt_v) {
			(void)printf("Test vnacommon_qrsolve: trial %3d "
				"A size %d x %d, B size %d x %d\n",
				trial, m, n, n, o);
			(void)fflush(stdout);
		    }

		    /*
		     * Generate random matrices A and B, and solve for X.
		     */
		    rank = qrsolve_helper(&x[0][0], &a[0][0], &b[0][0],
			    m, n, o);
		    if (opt_v) {
			libt_print_cmatrix("a", *a, m, n);
			libt_print_cmatrix("b", *b, m, o);
			libt_print_cmatrix("x", *x, n, o);
			(void)printf("rank %d\n", rank);
			(void)fflush(stdout);
		    }

		    /*
		     * Get the squared error of the result, then perturb
		     * each X_{j,k} value and verify that the error
		     * doesn't decrease when moving away from the result.
		     */
		    error0 = find_axb_error(&a[0][0], &x[0][0], &b[0][0],
			    m, n, o);
		    for (int k = 0; k < o; ++k) {
			for (int j = 0; j < n; ++j) {
			    static const double complex deltas[] =
				{ 0.001, 0.001 * I, -0.001, -0.001 * I };
			    double complex x0 = x[j][k];

			    for (int i = 0; i < 4; ++i) {
				double e;

				x[j][k] = x0 + deltas[i];
				e = find_axb_error(&a[0][0], &x[0][0], &b[0][0],
					m, n, o);
				if (e < error0) {
				    if (opt_a) {
					assert(!"bad result");
				    }
				    result = T_FAIL;
				    goto out;
				}
			    }
			    x[j][k] = x0;	/* restore x[j][k] */
			}
		    }

		    /*
		     * Check rank.
		     */
		    if (rank != n) {
			if (opt_a) {
			    assert(!"incorrect rank");
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
    exit(test_vnacommon_qrsolve());
}
