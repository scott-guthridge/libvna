/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2024 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include "vnaconv_internal.h"
#include "libt.h"
#include "libt_crand.h"


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
 * TEST_EQUAL: fail the test if x and y are not equal
 *   Assumes local variable "result" and label "out" are defined.
 */
#define TEST_EQUAL(x, y, label) \
    if (opt_a) { \
	assert(libt_isequal_label((x), (y), (label))); \
    } else { \
	if (!libt_isequal_label((x), (y), (label))) { \
	    result = T_FAIL; \
	    goto out; \
	} \
    } \

/*
 * run_trial: test one trial of vnaconv_stosrn()
 */
static libt_result_t run_trial(int n)
{
    double complex z1[n];
    double complex z2[n];
    double complex a1[n];
    double complex b1[n];
    double complex a2[n];
    double complex b2[n];
    double complex s1[n][n];
    double complex s2[n][n];
    double complex x[n];
    libt_result_t result = T_SKIPPED;

    for (int i = 0; i < n; ++i) {
	z1[i] = libt_crandn();
	z2[i] = libt_crandn();
	a1[i] = libt_crandn();
	for (int j = 0; j < n; ++j) {
	    s1[i][j] = libt_crandn();
	}
    }
    _vnacommon_mmultiply(b1, *s1, a1, n, n, 1);

    for (int i = 0; i < n; ++i) {
	double complex v, c;
	double complex z, zc;
	double zr, ki;

	z = z1[i];
	zc = conj(z);
	zr = creal(z);
	ki = sqrt(fabs(zr));
	v = ki * (zc * a1[i] + z * b1[i]) / zr;
	c = ki * (a1[i] - b1[i]) / zr;

	z = z2[i];
	zc = conj(z);
	zr = creal(z);
	ki = sqrt(fabs(zr));
	a2[i] = 0.5 * (v +  z * c) / ki;
	b2[i] = 0.5 * (v - zc * c) / ki;
    }
    if (opt_v) {
	libt_print_cmatrix("z1", z1, n, 1);
	libt_print_cmatrix("z2", z2, n, 1);
	libt_print_cmatrix("a1", a1, n, 1);
	libt_print_cmatrix("b1", b1, n, 1);
	libt_print_cmatrix("a2", a2, n, 1);
	libt_print_cmatrix("b2", b2, n, 1);
	libt_print_cmatrix("s1", *s1, n, n);
    }

    vnaconv_stosrn(*s1, *s2, z1, z2, n);
    if (opt_v) {
	libt_print_cmatrix("s2", *s2, n, n);
    }
    _vnacommon_mmultiply(x, *s2, a2, n, n, 1);
    for (int i = 0; i < n; ++i) {
	char name[3 * sizeof(int) + 5];

	(void)sprintf(name, "b2[%d]", i + 1);
	TEST_EQUAL(x[i], b2[i], name);
    }
    result = T_PASS;

out:
    return result;
}

/*
 * test_renormalize_NxN: test vnaconv_stosrn()
 */
static libt_result_t test_renormalize_NxN()
{
    libt_result_t result = T_SKIPPED;

    for (int trial = 0; trial < 10000; ++trial) {
	if (opt_v) {
	    (void)printf("Test renormalize NxN: trial %3d\n", trial);
	}
	for (int n = 1; n <= 5; ++n) {
	    result = run_trial(n);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	if (opt_v) {
	    printf("-------------\n");
	}
    }

out:
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
	    opt_a = 1;
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
    exit(test_renormalize_NxN());
}
