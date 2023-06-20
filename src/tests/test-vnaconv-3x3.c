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
#include "vnaconv_internal.h"
#include "libt.h"
#include "libt_crand.h"


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

#define X11	x[0][0]
#define X12	x[0][1]
#define X13	x[0][2]
#define X21	x[1][0]
#define X22	x[1][1]
#define X23	x[1][2]
#define X31	x[2][0]
#define X32	x[2][1]
#define X33	x[2][2]


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
 * test_conversions_3x3: test parameter conversions
 */
static libt_result_t test_conversions_3x3()
{
    libt_result_t result = T_SKIPPED;

    for (int trial = 0; trial < 10000; ++trial) {
	double k1i, k2i, k3i;
	double complex z0[3];
	double complex z1c, z2c, z3c;
	double complex a1, a2, a3, b1, b2, b3;
	double complex v1, i1, v2, i2, v3, i3;
	double complex s[3][3];
	double complex z[3][3];
	double complex y[3][3];
	double complex x[3][3];
	double complex xi[3];
	double complex zi[3];

	/*
	 * Set up test values.
	 */
	Z1  = libt_crandn();
	Z2  = libt_crandn();
	Z3  = libt_crandn();
	z1c = conj(Z1);
	z2c = conj(Z2);
	z3c = conj(Z3);
	k1i = sqrt(fabs(creal(Z1)));
	k2i = sqrt(fabs(creal(Z2)));
	k3i = sqrt(fabs(creal(Z3)));
	a1  = libt_crandn();
	a2  = libt_crandn();
	a3  = libt_crandn();
	S11 = libt_crandn();
	S12 = libt_crandn();
	S13 = libt_crandn();
	S21 = libt_crandn();
	S22 = libt_crandn();
	S23 = libt_crandn();
	S31 = libt_crandn();
	S32 = libt_crandn();
	S33 = libt_crandn();
	b1 = S11 * a1 + S12 * a2 + S13 * a3;
	b2 = S21 * a1 + S22 * a2 + S23 * a3;
	b3 = S31 * a1 + S32 * a2 + S33 * a3;
	v1 = k1i * (z1c * a1 + Z1 * b1) / creal(Z1);
	v2 = k2i * (z2c * a2 + Z2 * b2) / creal(Z2);
	v3 = k3i * (z3c * a3 + Z3 * b3) / creal(Z3);
	i1 = k1i * (a1 - b1) / creal(Z1);
	i2 = k2i * (a2 - b2) / creal(Z2);
	i3 = k3i * (a3 - b3) / creal(Z3);

	/*
	 * Calculate input impedance looking into each DUT port assuming
	 * that the other ports are terminated in their system impedances,
	 * i.e. not driven.  Because of this definition, it's not simply
	 * v1 / i1, v2 / i2, etc. which would be the effective impedance
	 * given that the other ports are also driven.
	 */
	zi[0] = (S11 * Z1 + z1c) / (1.0 - S11);
	zi[1] = (S22 * Z2 + z2c) / (1.0 - S22);
	zi[2] = (S33 * Z3 + z3c) / (1.0 - S33);

	if (opt_v) {
	    (void)printf("Test conversions: trial %3d\n",
		    trial);
	    (void)printf("Z1 %9.5f%+9.5fj  Z2 %9.5f%+9.5fj  Z3 %9.5f%+9.5fj\n",
		creal(Z1), cimag(Z1), creal(Z2), cimag(Z2),
		creal(Z3), cimag(Z3));
	    (void)printf("a1 %9.5f%+9.5fj  b1 %9.5f%+9.5fj\n",
		creal(a1), cimag(a1), creal(b1), cimag(b1));
	    (void)printf("a2 %9.5f%+9.5fj  b2 %9.5f%+9.5fj\n",
		creal(a2), cimag(a2), creal(b2), cimag(b2));
	    (void)printf("a3 %9.5f%+9.5fj  be %9.5f%+9.5fj\n",
		creal(a3), cimag(a3), creal(b3), cimag(b3));
	    (void)printf("v1 %9.5f%+9.5fj  i1 %9.5f%+9.5fj\n",
		creal(v1), cimag(v1), creal(i1), cimag(i1));
	    (void)printf("v2 %9.5f%+9.5fj  i2 %9.5f%+9.5fj\n",
		creal(v2), cimag(v2), creal(i2), cimag(i2));
	    (void)printf("v3 %9.5f%+9.5fj  i3 %9.5f%+9.5fj\n",
		creal(v3), cimag(v3), creal(i3), cimag(i3));
	    libt_print_cmatrix("zi", zi, 3, 1);
	    (void)printf("\n");
	    libt_print_cmatrix("s", *s, 3, 3);
	}
	TEST_EQUAL(S11 * a1 + S12 * a2 + S13 * a3, b1, "S11,S12,S13");
	TEST_EQUAL(S21 * a1 + S22 * a2 + S23 * a3, b2, "S21,S22,S23");
	TEST_EQUAL(S31 * a1 + S32 * a2 + S33 * a3, b3, "S31,S32,S33");

	vnaconv_stozn(*s, *z, z0, 3);
	if (opt_v) {
	    libt_print_cmatrix("z", *z, 3, 3);
	}
	TEST_EQUAL(Z11 * i1 + Z12 * i2 + Z13 * i3, v1, "stoz: Z11,Z12,Z13");
	TEST_EQUAL(Z21 * i1 + Z22 * i2 + Z23 * i3, v2, "stoz: Z21,Z22,Z23");
	TEST_EQUAL(Z31 * i1 + Z32 * i2 + Z33 * i3, v3, "stoz: Z31,Z32,Z33");

	vnaconv_stoyn(*s, *y, z0, 3);
	if (opt_v) {
	    libt_print_cmatrix("y", *y, 3, 3);
	}
	TEST_EQUAL(Y11 * v1 + Y12 * v2 + Y13 * v3, i1, "stoy: Y11,Y12,Y13");
	TEST_EQUAL(Y21 * v1 + Y22 * v2 + Y23 * v3, i2, "stoy: Y21,Y22,Y23");
	TEST_EQUAL(Y31 * v1 + Y32 * v2 + Y33 * v3, i3, "stoy: Y31,Y32,Y33");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ztosn(*z, *x, z0, 3);
	TEST_EQUAL(X11, S11, "ztos: S11");
	TEST_EQUAL(X12, S12, "ztos: S12");
	TEST_EQUAL(X13, S13, "ztos: S13");
	TEST_EQUAL(X21, S21, "ztos: S21");
	TEST_EQUAL(X22, S22, "ztos: S22");
	TEST_EQUAL(X23, S23, "ztos: S23");
	TEST_EQUAL(X31, S31, "ztos: S31");
	TEST_EQUAL(X32, S32, "ztos: S32");
	TEST_EQUAL(X33, S33, "ztos: S33");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ztoyn(*z, *x, 3);
	TEST_EQUAL(X11, Y11, "ztoy: Y11");
	TEST_EQUAL(X12, Y12, "ztoy: Y12");
	TEST_EQUAL(X13, Y13, "ztoy: Y13");
	TEST_EQUAL(X21, Y21, "ztoy: Y21");
	TEST_EQUAL(X22, Y22, "ztoy: Y22");
	TEST_EQUAL(X23, Y23, "ztoy: Y23");
	TEST_EQUAL(X31, Y31, "ztoy: Y31");
	TEST_EQUAL(X32, Y32, "ztoy: Y32");
	TEST_EQUAL(X33, Y33, "ztoy: Y33");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ytosn(*y, *x, z0, 3);
	TEST_EQUAL(X11, S11, "ytos: S11");
	TEST_EQUAL(X12, S12, "ytos: S12");
	TEST_EQUAL(X13, S13, "ytos: S13");
	TEST_EQUAL(X21, S21, "ytos: S21");
	TEST_EQUAL(X22, S22, "ytos: S22");
	TEST_EQUAL(X23, S23, "ytos: S23");
	TEST_EQUAL(X31, S31, "ytos: S31");
	TEST_EQUAL(X32, S32, "ytos: S32");
	TEST_EQUAL(X33, S33, "ytos: S33");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ytozn(*y, *x, 3);
	TEST_EQUAL(X11, Z11, "ytoz: Z11");
	TEST_EQUAL(X12, Z12, "ytoz: Z12");
	TEST_EQUAL(X13, Z13, "ytoz: Z13");
	TEST_EQUAL(X21, Z21, "ytoz: Z21");
	TEST_EQUAL(X22, Z22, "ytoz: Z22");
	TEST_EQUAL(X23, Z23, "ytoz: Z23");
	TEST_EQUAL(X31, Z31, "ytoz: Z31");
	TEST_EQUAL(X32, Z32, "ytoz: Z32");
	TEST_EQUAL(X33, Z33, "ytoz: Z33");

	(void)memset((void *)xi, 0, sizeof(xi));
	vnaconv_stozin(*s, xi, z0, 3);
	TEST_EQUAL(xi[0], zi[0], "stozin: Zi1");
	TEST_EQUAL(xi[1], zi[1], "stozin: Zi2");
	TEST_EQUAL(xi[2], zi[2], "stozin: Zi3");

	(void)memset((void *)xi, 0, sizeof(xi));
	vnaconv_ztozin(*z, xi, z0, 3);
	TEST_EQUAL(xi[0], zi[0], "ztozin: Zi1");
	TEST_EQUAL(xi[1], zi[1], "ztozin: Zi2");
	TEST_EQUAL(xi[2], zi[2], "ztozin: Zi3");

	(void)memset((void *)xi, 0, sizeof(xi));
	vnaconv_ytozin(*y, xi, z0, 3);
	TEST_EQUAL(xi[0], zi[0], "ytozin: Zi1");
	TEST_EQUAL(xi[1], zi[1], "ytozin: Zi2");
	TEST_EQUAL(xi[2], zi[2], "ytozin: Zi3");


	if (opt_v) {
	    printf("-------------\n");
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
    exit(test_conversions_3x3());
}
