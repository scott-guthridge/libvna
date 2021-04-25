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
 * qrsolve_q_helper: generate a random system of equations and solve
 */
static int qrsolve_q_helper(double complex *x, double complex *a,
	double complex *b, double complex *q, int m, int n, int o)
{
    double complex u[m][n];
    double complex v[m][o];

    /*
     * Generate random matrices A and B, and make copies
     * in U and V, respectively.
     */
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    a[i * n + j] = u[i][j] = test_crandn();
	}
	for (int k = 0; k < o; ++k) {
	    b[i * o + k] = v[i][k] = test_crandn();
	}
    }

    /*
     * Solve the system.  This call destroys both u and v.
     */
    return _vnacommon_qrsolve_q(x, &u[0][0], &v[0][0], q, m, n, o);
}

/*
 * test_vnacommon_qrsolve_q: test QR decomposition with Q
 */
static test_result_t test_vnacommon_qrsolve_q()
{
    test_result_t result = T_SKIPPED;

    for (int trial = 1; trial <= N_MATRIX_TRIALS; ++trial) {
	for (int m = 1; m <= 5; ++m) {
	    for (int n = 1; n <= 5; ++n) {
		for (int o = 1; o <= 2; ++o) {
		    double complex a[m][n];
		    double complex b[m][o];
		    double complex x[n][o];
		    double complex q[m][m];
		    int diagonals = MIN(m, n);
		    int rank;

		    /*
		     * If -v, print the test header.
		     */
		    if (opt_v) {
			(void)printf("Test vnacommon_qrsolve_q: trial %3d "
				"size %d x %d\n", trial, m, n);
			(void)fflush(stdout);
		    }

		    /*
		     * Generate random matrices A and B, and solve for X.
		     */
		    rank = qrsolve_q_helper(*x, *a, *b, *q, m, n, o);
		    if (opt_v) {
			test_print_cmatrix("a", *a, m, n);
			test_print_cmatrix("b", *b, m, o);
			test_print_cmatrix("x", *x, n, o);
			test_print_cmatrix("q", *q, m, m);
			(void)printf("rank %d\n", rank);
			(void)fflush(stdout);
		    }

		    /*
		     * If m <= n, verify A X == B.  Otherwise, the system
		     * is overdetermined and the equality won't hold.
		     */
		    if (m <= n) {
			for (int k = 0; k < o; ++k) {
			    for (int i = 0; i < m; ++i) {
				double complex s = 0.0;

				for (int j = 0; j < n; ++j) {
				    s += a[i][j] * x[j][k];
				}
				if (!test_isequal(s, b[i][k])) {
				    if (opt_a) {
					assert(!"data miscompare");
				    }
				    result = T_FAIL;
				    goto out;
				}
			    }
			}
		    }

		    /*
		     * Check rank.
		     */
		    if (rank != diagonals) {
			if (opt_a) {
			    assert(!"incorrect rank");
			}
			result = T_FAIL;
			goto out;
		    }

#if 1//ZZ
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
#endif
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
    exit(test_vnacommon_qrsolve_q());
}
