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
#include "vnaconv_internal.h"
#include "test.h"


#define Z1	z0[0]
#define Z2	z0[1]
#define Z3	z0[2]

#define S11	s[0][0]
#define S12	s[0][1]
#define S13	s[0][2]
#define S21	s[1][0]
#define S22	s[1][1]
#define S23	s[1][2]
#define S31	s[2][0]
#define S32	s[2][1]
#define S33	s[2][2]

#define Z11	z[0][0]
#define Z12	z[0][1]
#define Z13	z[0][2]
#define Z21	z[1][0]
#define Z22	z[1][1]
#define Z23	z[1][2]
#define Z31	z[2][0]
#define Z32	z[2][1]
#define Z33	z[2][2]


#define Y11	y[0][0]
#define Y12	y[0][1]
#define Y13	y[0][2]
#define Y21	y[1][0]
#define Y22	y[1][1]
#define Y23	y[1][2]
#define Y31	y[2][0]
#define Y32	y[2][1]
#define Y33	y[2][2]

#define U11	u[0][0]
#define U12	u[0][1]
#define U13	u[0][2]
#define U21	u[1][0]
#define U22	u[1][1]
#define U23	u[1][2]
#define U31	u[2][0]
#define U32	u[2][1]
#define U33	u[2][2]


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
 * TEST_EQUAL: fail the test if x and y are not equal
 *   Assumes local variable "result" and label "out" are defined.
 */
#define TEST_EQUAL(x, y, label) \
    if (opt_a) { \
	assert(test_isequal_label((x), (y), (label))); \
    } else { \
	if (!test_isequal_label((x), (y), (label))) { \
	    result = T_FAIL; \
	    goto out; \
	} \
    } \

/*
 * test_conversions_3x3: test parameter conversions
 */
static test_result_t test_conversions_3x3()
{
    test_result_t result = T_SKIPPED;

    for (int trial = 0; trial < 10000; ++trial) {
	double k1i, k2i, k3i;
	double complex z0[3];
	double complex z1c, z2c, z3c;
	double complex a1, a2, a3, b1, b2, b3;
	double complex v1, i1, v2, i2, v3, i3;
	double complex s[3][3];
	double complex z[3][3];
	double complex y[3][3];
	double complex u[3][3];
	double complex v[3];
	double complex zi[3];

	Z1  = test_crandn();
	Z2  = test_crandn();
	Z3  = test_crandn();
	z1c = conj(Z1);
	z2c = conj(Z2);
	z3c = conj(Z3);
	k1i = sqrt(fabs(creal(Z1)));
	k2i = sqrt(fabs(creal(Z2)));
	k3i = sqrt(fabs(creal(Z3)));
	a1  = test_crandn();
	a2  = test_crandn();
	a3  = test_crandn();
	S11 = test_crandn();
	S12 = test_crandn();
	S13 = test_crandn();
	S21 = test_crandn();
	S22 = test_crandn();
	S23 = test_crandn();
	S31 = test_crandn();
	S32 = test_crandn();
	S33 = test_crandn();
	b1 = S11 * a1 + S12 * a2 + S13 * a3;
	b2 = S21 * a1 + S22 * a2 + S23 * a3;
	b3 = S31 * a1 + S32 * a2 + S33 * a3;
	v1 = k1i * (z1c * a1 + Z1 * b1) / creal(Z1);
	v2 = k2i * (z2c * a2 + Z2 * b2) / creal(Z2);
	v3 = k3i * (z3c * a3 + Z3 * b3) / creal(Z3);
	i1 = k1i * (a1 - b1) / creal(Z1);
	i2 = k2i * (a2 - b2) / creal(Z2);
	i3 = k3i * (a3 - b3) / creal(Z3);

	if (opt_v) {
	    (void)printf("Test conversions: trial %3d\n",
		    trial);
	    (void)printf("Z1 %9.5f%+9.5fj  Z2 %9.5f%+9.5fj  Z3 %9.5f%+9.5fj\n",
		creal(Z1), cimag(Z1), creal(Z2), cimag(Z2),
		creal(Z3), cimag(Z3));
	    (void)printf("v1 %9.5f%+9.5fj  i1 %9.5f%+9.5fj\n",
		creal(v1), cimag(v1), creal(i1), cimag(i1));
	    (void)printf("v2 %9.5f%+9.5fj  i2 %9.5f%+9.5fj\n",
		creal(v2), cimag(v2), creal(i2), cimag(i2));
	    (void)printf("v3 %9.5f%+9.5fj  i3 %9.5f%+9.5fj\n",
		creal(v3), cimag(v3), creal(i3), cimag(i3));
	    (void)printf("\n");
	    test_print_cmatrix("s", *s, 3, 3);
	}
	TEST_EQUAL(S11 * a1 + S12 * a2 + S13 * a3, b1, "S11,S12,S13");
	TEST_EQUAL(S21 * a1 + S22 * a2 + S23 * a3, b2, "S21,S22,S23");
	TEST_EQUAL(S31 * a1 + S32 * a2 + S33 * a3, b3, "S31,S32,S33");

	vnaconv_stozn(*s, *z, z0, 3);
	if (opt_v) {
	    test_print_cmatrix("z", *z, 3, 3);
	}
	TEST_EQUAL(Z11 * i1 + Z12 * i2 + Z13 * i3, v1, "stoz: Z11,Z12,Z13");
	TEST_EQUAL(Z21 * i1 + Z22 * i2 + Z23 * i3, v2, "stoz: Z21,Z22,Z23");
	TEST_EQUAL(Z31 * i1 + Z32 * i2 + Z33 * i3, v3, "stoz: Z31,Z32,Z33");

	vnaconv_stoyn(*s, *y, z0, 3);
	if (opt_v) {
	    test_print_cmatrix("y", *y, 3, 3);
	}
	TEST_EQUAL(Y11 * v1 + Y12 * v2 + Y13 * v3, i1, "stoy: Y11,Y12,Y13");
	TEST_EQUAL(Y21 * v1 + Y22 * v2 + Y23 * v3, i2, "stoy: Y21,Y22,Y23");
	TEST_EQUAL(Y31 * v1 + Y32 * v2 + Y33 * v3, i3, "stoy: Y31,Y32,Y33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztosn(*z, *u, z0, 3);
	TEST_EQUAL(U11, S11, "ztos: S11");
	TEST_EQUAL(U12, S12, "ztos: S12");
	TEST_EQUAL(U13, S13, "ztos: S13");
	TEST_EQUAL(U21, S21, "ztos: S21");
	TEST_EQUAL(U22, S22, "ztos: S22");
	TEST_EQUAL(U23, S23, "ztos: S23");
	TEST_EQUAL(U31, S31, "ztos: S31");
	TEST_EQUAL(U32, S32, "ztos: S32");
	TEST_EQUAL(U33, S33, "ztos: S33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztoyn(*z, *u, 3);
	TEST_EQUAL(U11, Y11, "ztoy: Y11");
	TEST_EQUAL(U12, Y12, "ztoy: Y12");
	TEST_EQUAL(U13, Y13, "ztoy: Y13");
	TEST_EQUAL(U21, Y21, "ztoy: Y21");
	TEST_EQUAL(U22, Y22, "ztoy: Y22");
	TEST_EQUAL(U23, Y23, "ztoy: Y23");
	TEST_EQUAL(U31, Y31, "ztoy: Y31");
	TEST_EQUAL(U32, Y32, "ztoy: Y32");
	TEST_EQUAL(U33, Y33, "ztoy: Y33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytosn(*y, *u, z0, 3);
	TEST_EQUAL(U11, S11, "ytos: S11");
	TEST_EQUAL(U12, S12, "ytos: S12");
	TEST_EQUAL(U13, S13, "ytos: S13");
	TEST_EQUAL(U21, S21, "ytos: S21");
	TEST_EQUAL(U22, S22, "ytos: S22");
	TEST_EQUAL(U23, S23, "ytos: S23");
	TEST_EQUAL(U31, S31, "ytos: S31");
	TEST_EQUAL(U32, S32, "ytos: S32");
	TEST_EQUAL(U33, S33, "ytos: S33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytozn(*y, *u, 3);
	TEST_EQUAL(U11, Z11, "ytoz: Z11");
	TEST_EQUAL(U12, Z12, "ytoz: Z12");
	TEST_EQUAL(U13, Z13, "ytoz: Z13");
	TEST_EQUAL(U21, Z21, "ytoz: Z21");
	TEST_EQUAL(U22, Z22, "ytoz: Z22");
	TEST_EQUAL(U23, Z23, "ytoz: Z23");
	TEST_EQUAL(U31, Z31, "ytoz: Z31");
	TEST_EQUAL(U32, Z32, "ytoz: Z32");
	TEST_EQUAL(U33, Z33, "ytoz: Z33");

	vnaconv_stozin(*s, zi, z0, 3);
	if (opt_v) {
	    test_print_cmatrix("zi", zi, 3, 1);
	}

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_stozimn(*s, v, z0, 3, 3);
	TEST_EQUAL(v[0], zi[0], "stozimn: zi0");
	TEST_EQUAL(v[1], zi[1], "stozimn: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_ztozin(*z, v, z0, 3);
	TEST_EQUAL(v[0], zi[0], "ztozi: zi0");
	TEST_EQUAL(v[1], zi[1], "ztozi: zi1");
	TEST_EQUAL(v[2], zi[2], "ztozi: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_ytozin(*y, v, z0, 3);
	TEST_EQUAL(v[0], zi[0], "ytozi: zi0");
	TEST_EQUAL(v[1], zi[1], "ytozi: zi1");
	TEST_EQUAL(v[2], zi[2], "ytozi: zi1");


	if (opt_v) {
	    printf("-------------\n");
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
	    opt_a = 1;
	    continue;

	case 'v':
	    opt_v = 1;
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
    exit(test_conversions_3x3());
}
