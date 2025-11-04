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
#include "vnaconv_internal.h"
#include "libt.h"
#include "libt_crand.h"


#define Z11	z1[0]	/* original z0 */
#define Z12	z1[1]
#define Z21	z2[0]	/* new z0 */
#define Z22	z2[1]

#define S1_11	s1[0][0]
#define S1_12	s1[0][1]
#define S1_21	s1[1][0]
#define S1_22	s1[1][1]

#define S2_11	s2[0][0]
#define S2_12	s2[0][1]
#define S2_21	s2[1][0]
#define S2_22	s2[1][1]

#define T1_11	t1[0][0]
#define T1_12	t1[0][1]
#define T1_21	t1[1][0]
#define T1_22	t1[1][1]

#define T2_11	t2[0][0]
#define T2_12	t2[0][1]
#define T2_21	t2[1][0]
#define T2_22	t2[1][1]

#define U1_11	u1[0][0]
#define U1_12	u1[0][1]
#define U1_21	u1[1][0]
#define U1_22	u1[1][1]

#define U2_11	u2[0][0]
#define U2_12	u2[0][1]
#define U2_21	u2[1][0]
#define U2_22	u2[1][1]

#define X11	x[0][0]
#define X12	x[0][1]
#define X21	x[1][0]
#define X22	x[1][1]


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
    }

/*
 * test_conversions_2x2: test parameter conversions
 */
static libt_result_t test_renormalize_2x2()
{
    libt_result_t result = T_SKIPPED;

    for (int trial = 0; trial < 10000; ++trial) {
	double k11i, k12i;
	double k21i, k22i;
	double complex z1[2];
	double complex z2[2];
	double complex z11c, z12c;
	double complex z21c, z22c;
	double complex a11, a12, b11, b12;
	double complex a21, a22, b21, b22;
	double complex v1, i1, v2, i2;
	double complex s1[2][2], s2[2][2];
	double complex t1[2][2], t2[2][2];
	double complex u1[2][2], u2[2][2];
	double complex x[2][2];

	/*
	 * Set up test values.
	 */
	Z11  = libt_crandn();
	Z12  = libt_crandn();
	Z21  = libt_crandn();
	Z22  = libt_crandn();
	z11c = conj(Z11);
	z12c = conj(Z12);
	z21c = conj(Z21);
	z22c = conj(Z22);
	k11i = sqrt(fabs(creal(Z11)));
	k12i = sqrt(fabs(creal(Z12)));
	k21i = sqrt(fabs(creal(Z21)));
	k22i = sqrt(fabs(creal(Z22)));
	a11  = libt_crandn();
	a12  = libt_crandn();
	S1_11 = libt_crandn();
	S1_12 = libt_crandn();
	S1_21 = libt_crandn();
	S1_22 = libt_crandn();
	b11 = S1_11 * a11 + S1_12 * a12;
	b12 = S1_22 * a12 + S1_21 * a11;
	v1 = k11i * (z11c * a11 + Z11 * b11) / creal(Z11);
	i1 = k11i * (a11 - b11) / creal(Z11);
	v2 = k12i * (z12c * a12 + Z12 * b12) / creal(Z12);
	i2 = k12i * (a12 - b12) / creal(Z12);
	a21 = 0.5 * (v1 + Z21 * i1)  / k21i;
	b21 = 0.5 * (v1 - z21c * i1) / k21i;
	a22 = 0.5 * (v2 + Z22 * i2)  / k22i;
	b22 = 0.5 * (v2 - z22c * i2) / k22i;

	if (opt_v) {
	    (void)printf("Test renormalize 2x2: trial %3d\n",
		    trial);
	    (void)printf("Z11 %9.5f%+9.5fj  Z12 %9.5f%+9.5fj\n",
		creal(Z11), cimag(Z11), creal(Z12), cimag(Z12));
	    (void)printf("Z21 %9.5f%+9.5fj  Z22 %9.5f%+9.5fj\n",
		creal(Z21), cimag(Z21), creal(Z22), cimag(Z22));
	    (void)printf("a11 %9.5f%+9.5fj  b11 %9.5f%+9.5fj\n",
		creal(a11), cimag(a11), creal(b11), cimag(b11));
	    (void)printf("a12 %9.5f%+9.5fj  b12 %9.5f%+9.5fj\n",
		creal(a12), cimag(a12), creal(b12), cimag(b12));
	    (void)printf("v1 %9.5f%+9.5fj  i1 %9.5f%+9.5fj\n",
		creal(v1), cimag(v1), creal(i1), cimag(i1));
	    (void)printf("v2 %9.5f%+9.5fj  i2 %9.5f%+9.5fj\n",
		creal(v2), cimag(v2), creal(i2), cimag(i2));
	    (void)printf("a21 %9.5f%+9.5fj  b21 %9.5f%+9.5fj\n",
		creal(a21), cimag(a21), creal(b21), cimag(b21));
	    (void)printf("a22 %9.5f%+9.5fj  b22 %9.5f%+9.5fj\n",
		creal(a22), cimag(a22), creal(b22), cimag(b22));
	    (void)printf("\n");
	    libt_print_cmatrix("s1", *s1, 2, 2);
	}
	TEST_EQUAL(S1_11 * a11 + S1_12 * a12, b11, "S1_11,S1_12");
	TEST_EQUAL(S1_21 * a11 + S1_22 * a12, b12, "S1_21,S1_22");

	vnaconv_stot(s1, t1);
	if (opt_v) {
	    libt_print_cmatrix("t1", *t1, 2, 2);
	}
	TEST_EQUAL(T1_11 * a12 + T1_12 * b12, b11, "stot: T1_11,T1_12");
	TEST_EQUAL(T1_21 * a12 + T1_22 * b12, a11, "stot: T1_22,T1_22");

	vnaconv_stou(s1, u1);
	if (opt_v) {
	    libt_print_cmatrix("u1", *u1, 2, 2);
	}
	TEST_EQUAL(U1_11 * b11 + U1_12 * a11, a12, "stou: U1_11,U1_12");
	TEST_EQUAL(U1_21 * b11 + U1_22 * a11, b12, "stou: U1_22,U1_22");

	vnaconv_stosr(s1, s2, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("s2", *s2, 2, 2);
	}
	TEST_EQUAL(S2_11 * a21 + S2_12 * a22, b21, "stosr: S2_11,S2_12");
	TEST_EQUAL(S2_21 * a21 + S2_22 * a22, b22, "stosr: S2_21,S2_22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_stosrn(*s1, *x, z1, z2, 2);
	if (opt_v) {
	    libt_print_cmatrix("stosrn", *x, 2, 2);
	}
	TEST_EQUAL(X11, S2_11, "stosrn: X11");
	TEST_EQUAL(X12, S2_12, "stosrn: X12");
	TEST_EQUAL(X21, S2_21, "stosrn: X21");
	TEST_EQUAL(X22, S2_22, "stosrn: X22");

	vnaconv_stotr(s1, t2, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("t2", *t2, 2, 2);
	}
	TEST_EQUAL(T2_11 * a22 + T2_12 * b22, b21, "stotr: T2_11,T2_12");
	TEST_EQUAL(T2_21 * a22 + T2_22 * b22, a21, "stotr: T2_22,T2_22");

	vnaconv_stour(s1, u2, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("u2", *u2, 2, 2);
	}
	TEST_EQUAL(U2_11 * b21 + U2_12 * a21, a22, "stou: U2_11,U2_12");
	TEST_EQUAL(U2_21 * b21 + U2_22 * a21, b22, "stou: U2_22,U2_22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ttosr(t1, x, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("ttosr", *x, 2, 2);
	}
	TEST_EQUAL(X11, S2_11, "ttosr: X11");
	TEST_EQUAL(X12, S2_12, "ttosr: X12");
	TEST_EQUAL(X21, S2_21, "ttosr: X21");
	TEST_EQUAL(X22, S2_22, "ttosr: X22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ttotr(t1, x, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("ttotr", *x, 2, 2);
	}
	TEST_EQUAL(X11, T2_11, "ttotr: X11");
	TEST_EQUAL(X12, T2_12, "ttotr: X12");
	TEST_EQUAL(X21, T2_21, "ttotr: X21");
	TEST_EQUAL(X22, T2_22, "ttotr: X22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_ttour(t1, x, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("ttour", *x, 2, 2);
	}
	TEST_EQUAL(X11, U2_11, "ttour: X11");
	TEST_EQUAL(X12, U2_12, "ttour: X12");
	TEST_EQUAL(X21, U2_21, "ttour: X21");
	TEST_EQUAL(X22, U2_22, "ttour: X22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_utosr(u1, x, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("utosr", *x, 2, 2);
	}
	TEST_EQUAL(X11, S2_11, "utosr: X11");
	TEST_EQUAL(X12, S2_12, "utosr: X12");
	TEST_EQUAL(X21, S2_21, "utosr: X21");
	TEST_EQUAL(X22, S2_22, "utosr: X22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_utotr(u1, x, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("utotr", *x, 2, 2);
	}
	TEST_EQUAL(X11, T2_11, "utotr: X11");
	TEST_EQUAL(X12, T2_12, "utotr: X12");
	TEST_EQUAL(X21, T2_21, "utotr: X21");
	TEST_EQUAL(X22, T2_22, "utotr: X22");

	(void)memset((void *)x, 0, sizeof(x));
	vnaconv_utour(u1, x, z1, z2);
	if (opt_v) {
	    libt_print_cmatrix("utour", *x, 2, 2);
	}
	TEST_EQUAL(X11, U2_11, "utour: X11");
	TEST_EQUAL(X12, U2_12, "utour: X12");
	TEST_EQUAL(X21, U2_21, "utour: X21");
	TEST_EQUAL(X22, U2_22, "utour: X22");

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
    exit(test_renormalize_2x2());
}
