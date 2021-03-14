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
#include <unistd.h>
#include "vnaconv_internal.h"


#define PI	3.1415926535897932384626433832795
#define EPS	1.0e-4

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

#define T11	t[0][0]
#define T12	t[0][1]
#define T21	t[1][0]
#define T22	t[1][1]

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

#define H11	h[0][0]
#define H12	h[0][1]
#define H21	h[1][0]
#define H22	h[1][1]

#define G11	g[0][0]
#define G12	g[0][1]
#define G21	g[1][0]
#define G22	g[1][1]

#define A11	a[0][0]
#define A12	a[0][1]
#define A21	a[1][0]
#define A22	a[1][1]

#define B11	b[0][0]
#define B12	b[0][1]
#define B21	b[1][0]
#define B22	b[1][1]

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
static bool opt_v = false;

/*
 * crandn: generate normally distributed complex numbers
 */
static double complex crandn()
{
    double u1 = (random() + 1.0) / RAND_MAX;
    double u2 = (double)random() / RAND_MAX;
    double r = sqrt(-2.0 * log(u1));
    double a = 2 * PI * u2;

    return r * (cos(a) + I * sin(a));
}

/*
 * cmatrix_print: print an m by n serialized complex matrix
 */
static void cmatrix_print(const char *tag, double complex *a, int m, int n)
{
    (void)printf("%s:\n", tag);
    for (int i = 0; i < m; ++i) {
	for (int j = 0; j < n; ++j) {
	    double complex v = a[i * n + j];

	    (void)printf(" %9.5f%+9.5fj", creal(v), cimag(v));
	}
	(void)printf("\n");
    }
    (void)printf("\n");
}

/*
 * isequal: test if x and y are approximately equal
 */
static int isequal(double complex x, double complex y, const char *label)
{
    int rv;
    double d = cabs(csqrt(x * y));

    if (d < 1.0) {
	d = 1.0;
    }
    rv = cabs(x - y) / d < EPS;
    if (!rv) {
	printf("%s: %9.5f%+9.5fj != %9.5f%+9.5fj\n",
		label, creal(x), cimag(x), creal(y), cimag(y));
	printf("|x-y| = %9.5f\n", cabs(x - y));
    }
    return rv;
}

/*
 * test_result_type
 */
typedef enum test_result {
    T_PASS,
    T_FAIL,
    T_SKIPPED
} test_result_type;

/*
 * TEST_EQUAL: fail the test if x and y are not equal
 *   Assumes local variable "result" and label "out" are defined.
 */
#define TEST_EQUAL(x, y, label) \
    if (opt_a) { \
	assert(isequal((x), (y), (label))); \
    } else { \
	if (!isequal((x), (y), (label))) { \
	    result = T_FAIL; \
	    goto out; \
	} \
    } \

/*
 * test counters
 */
static int test_count = 0;
static int fail_count = 0;

/*
 * report_test_result: report a test result
 */
static void report_test_result(const char *test_name, test_result_type result)
{
    const char *result_name;

    switch (result) {
    case T_PASS:
	result_name = "PASS";
	break;
    case T_FAIL:
	result_name = "FAIL";
	break;
    case T_SKIPPED:
	result_name = "SKIPPED";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)printf("Test %2d: %58s %s\n", ++test_count, test_name, result_name);
    (void)fflush(stdout);
    if (result == T_FAIL) {
	++fail_count;
    }
}

/*
 * test_conversions_2x2: test parameter conversions
 */
static void test_conversions_2x2()
{
    test_result_type result = T_SKIPPED;

    for (int trial = 0; trial < 10000; ++trial) {
	double k1i, k2i;
	double complex z0[2];
	double complex z1c, z2c;
	double complex a1, a2, b1, b2;
	double complex v1, i1, v2, i2;
	double complex s[2][2];
	double complex t[2][2];
	double complex z[2][2];
	double complex y[2][2];
	double complex h[2][2];
	double complex g[2][2];
	double complex a[2][2];
	double complex b[2][2];
	double complex u[2][2];
	double complex v[2];
	double complex zi[2];

	Z1  = crandn();
	Z2  = crandn();
	z1c = conj(Z1);
	z2c = conj(Z2);
	k1i = sqrt(fabs(creal(Z1)));
	k2i = sqrt(fabs(creal(Z2)));
	a1  = crandn();
	a2  = crandn();
	S11 = crandn();
	S12 = crandn();
	S21 = crandn();
	S22 = crandn();
	b1 = S11 * a1 + S12 * a2;
	b2 = S22 * a2 + S21 * a1;
	v1 = k1i * (z1c * a1 + Z1 * b1) / creal(Z1);
	v2 = k2i * (z2c * a2 + Z2 * b2) / creal(Z2);
	i1 = k1i * (a1 - b1) / creal(Z1);
	i2 = k2i * (a2 - b2) / creal(Z2);

	if (opt_v) {
	    (void)printf("Test conversions: trial %3d\n",
		    trial);
	    (void)printf("Z1 %9.5f%+9.5fj  Z2 %9.5f%+9.5fj\n",
		creal(Z1), cimag(Z1), creal(Z2), cimag(Z2));
	    (void)printf("v1 %9.5f%+9.5fj  i1 %9.5f%+9.5fj\n",
		creal(v1), cimag(v1), creal(i1), cimag(i1));
	    (void)printf("v2 %9.5f%+9.5fj  i2 %9.5f%+9.5fj\n",
		creal(v2), cimag(v2), creal(i2), cimag(i2));
	    (void)printf("\n");
	    cmatrix_print("s", *s, 2, 2);
	}
	TEST_EQUAL(S11 * a1 + S12 * a2, b1, "S11,S12");
	TEST_EQUAL(S21 * a1 + S22 * a2, b2, "S21,S22");

	vnaconv_stot(s, t);
	if (opt_v) {
	    cmatrix_print("t", *t, 2, 2);
	}
	TEST_EQUAL(T11 * a2 + T12 * b2, b1, "stot: T11,T12");
	TEST_EQUAL(T21 * a2 + T22 * b2, a1, "stot: T22,T22");

	vnaconv_stoz(s, z, z0);
	if (opt_v) {
	    cmatrix_print("z", *z, 2, 2);
	}
	TEST_EQUAL(Z11 * i1 + Z12 * i2, v1, "stoz: Z11,Z12");
	TEST_EQUAL(Z21 * i1 + Z22 * i2, v2, "stoz: Z21,Z22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_stozn(*s, *u, z0, 2);
	TEST_EQUAL(U11 * i1 + U12 * i2, v1, "stozn: U11,U12");
	TEST_EQUAL(U21 * i1 + U22 * i2, v2, "stozn: U21,U22");

	vnaconv_stoy(s, y, z0);
	if (opt_v) {
	    cmatrix_print("y", *y, 2, 2);
	}
	TEST_EQUAL(Y11 * v1 + Y12 * v2, i1, "stoy: Y11,Y12");
	TEST_EQUAL(Y21 * v1 + Y22 * v2, i2, "stoy: Y21,Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_stoyn(*s, *u, z0, 2);
	TEST_EQUAL(U11 * v1 + U12 * v2, i1, "stoyn: U11,U12");
	TEST_EQUAL(U21 * v1 + U22 * v2, i2, "stoyn: U21,U22");

	vnaconv_stoh(s, h, z0);
	if (opt_v) {
	    cmatrix_print("h", *h, 2, 2);
	}
	TEST_EQUAL(H11 * i1 + H12 * v2, v1, "stoh: H11,H12");
	TEST_EQUAL(H21 * i1 + H22 * v2, i2, "stoh: H21,H22");

	vnaconv_stog(s, g, z0);
	if (opt_v) {
	    cmatrix_print("g", *g, 2, 2);
	}
	TEST_EQUAL(G11 * v1 + G12 * i2, i1, "stog: G11,G12");
	TEST_EQUAL(G21 * v1 + G22 * i2, v2, "stog: G21,G22");

	vnaconv_stoa(s, a, z0);
	if (opt_v) {
	    cmatrix_print("a", *a, 2, 2);
	}
	TEST_EQUAL(A11 * v2 + A12 * -i2, v1, "stoa: A11,A12");
	TEST_EQUAL(A21 * v2 + A22 * -i2, i1, "stob: A21,A22");

	vnaconv_stob(s, b, z0);
	if (opt_v) {
	    cmatrix_print("b", *b, 2, 2);
	}
	TEST_EQUAL(B11 * v1 + B12 * i1,  v2, "stob: B11,B12");
	TEST_EQUAL(B21 * v1 + B22 * i1, -i2, "stob: B21,B22");

	vnaconv_stozi(s, zi, z0);
	if (opt_v) {
	    cmatrix_print("zi", zi, 2, 1);
	}

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_stozin(*s, v, z0, 2);
	TEST_EQUAL(v[0], zi[0], "stozin: zi0");
	TEST_EQUAL(v[1], zi[1], "stozin: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_stozimn(*s, v, z0, 2, 2);
	TEST_EQUAL(v[0], zi[0], "stozimn: zi0");
	TEST_EQUAL(v[1], zi[1], "stozimn: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttos(t, u);
	TEST_EQUAL(U11, S11, "ttos: S11");
	TEST_EQUAL(U12, S12, "ttos: S12");
	TEST_EQUAL(U21, S21, "ttos: S21");
	TEST_EQUAL(U22, S22, "ttos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttoz(t, u, z0);
	TEST_EQUAL(U11, Z11, "ttoz: Z11");
	TEST_EQUAL(U12, Z12, "ttoz: Z12");
	TEST_EQUAL(U21, Z21, "ttoz: Z21");
	TEST_EQUAL(U22, Z22, "ttoz: Z22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttoy(t, u, z0);
	TEST_EQUAL(U11, Y11, "ttoy: Y11");
	TEST_EQUAL(U12, Y12, "ttoy: Y12");
	TEST_EQUAL(U21, Y21, "ttoy: Y21");
	TEST_EQUAL(U22, Y22, "ttoy: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttoh(t, u, z0);
	TEST_EQUAL(U11, H11, "ttoh: H11");
	TEST_EQUAL(U12, H12, "ttoh: H12");
	TEST_EQUAL(U21, H21, "ttoh: H21");
	TEST_EQUAL(U22, H22, "ttoh: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttog(t, u, z0);
	TEST_EQUAL(U11, G11, "ttog: G11");
	TEST_EQUAL(U12, G12, "ttog: G12");
	TEST_EQUAL(U21, G21, "ttog: G21");
	TEST_EQUAL(U22, G22, "ttog: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttoa(t, u, z0);
	TEST_EQUAL(U11, A11, "ttoa: A11");
	TEST_EQUAL(U12, A12, "ttoa: A12");
	TEST_EQUAL(U21, A21, "ttoa: A21");
	TEST_EQUAL(U22, A22, "ttoa: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ttob(t, u, z0);
	TEST_EQUAL(U11, B11, "ttob: B11");
	TEST_EQUAL(U12, B12, "ttob: B12");
	TEST_EQUAL(U21, B21, "ttob: B21");
	TEST_EQUAL(U22, B22, "ttob: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_ttozi(t, v, z0);
	TEST_EQUAL(v[0], zi[0], "ttozi: zi0");
	TEST_EQUAL(v[1], zi[1], "ttozi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztos(z, u, z0);
	TEST_EQUAL(U11, S11, "ztos: S11");
	TEST_EQUAL(U12, S12, "ztos: S12");
	TEST_EQUAL(U21, S21, "ztos: S21");
	TEST_EQUAL(U22, S22, "ztos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztosn(*z, *u, z0, 2);
	TEST_EQUAL(U11, S11, "ztosn: S11");
	TEST_EQUAL(U12, S12, "ztosn: S12");
	TEST_EQUAL(U21, S21, "ztosn: S21");
	TEST_EQUAL(U22, S22, "ztosn: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztot(z, u, z0);
	TEST_EQUAL(U11, T11, "ztot: T11");
	TEST_EQUAL(U12, T12, "ztot: T12");
	TEST_EQUAL(U21, T21, "ztot: T21");
	TEST_EQUAL(U22, T22, "ztot: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztoy(z, u);
	TEST_EQUAL(U11, Y11, "ztoy: Y11");
	TEST_EQUAL(U12, Y12, "ztoy: Y12");
	TEST_EQUAL(U21, Y21, "ztoy: Y21");
	TEST_EQUAL(U22, Y22, "ztoy: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztoyn(*z, *u, 2);
	TEST_EQUAL(U11, Y11, "ztoyn: Y11");
	TEST_EQUAL(U12, Y12, "ztoyn: Y12");
	TEST_EQUAL(U21, Y21, "ztoyn: Y21");
	TEST_EQUAL(U22, Y22, "ztoyn: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztoh(z, u);
	TEST_EQUAL(U11, H11, "ztoh: H11");
	TEST_EQUAL(U12, H12, "ztoh: H12");
	TEST_EQUAL(U21, H21, "ztoh: H21");
	TEST_EQUAL(U22, H22, "ztoh: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztog(z, u);
	TEST_EQUAL(U11, G11, "ztog: G11");
	TEST_EQUAL(U12, G12, "ztog: G12");
	TEST_EQUAL(U21, G21, "ztog: G21");
	TEST_EQUAL(U22, G22, "ztog: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztoa(z, u);
	TEST_EQUAL(U11, A11, "ztoa: A11");
	TEST_EQUAL(U12, A12, "ztoa: A12");
	TEST_EQUAL(U21, A21, "ztoa: A21");
	TEST_EQUAL(U22, A22, "ztoa: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ztob(z, u);
	TEST_EQUAL(U11, B11, "ztob: B11");
	TEST_EQUAL(U12, B12, "ztob: B12");
	TEST_EQUAL(U21, B21, "ztob: B21");
	TEST_EQUAL(U22, B22, "ztob: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_ztozi(z, v, z0);
	TEST_EQUAL(v[0], zi[0], "ztozi: zi0");
	TEST_EQUAL(v[1], zi[1], "ztozi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytos(y, u, z0);
	TEST_EQUAL(U11, S11, "ytos: S11");
	TEST_EQUAL(U12, S12, "ytos: S12");
	TEST_EQUAL(U21, S21, "ytos: S21");
	TEST_EQUAL(U22, S22, "ytos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytot(y, u, z0);
	TEST_EQUAL(U11, T11, "ytot: T11");
	TEST_EQUAL(U12, T12, "ytot: T12");
	TEST_EQUAL(U21, T21, "ytot: T21");
	TEST_EQUAL(U22, T22, "ytot: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytoz(y, u);
	TEST_EQUAL(U11, Z11, "ytoz: Y11");
	TEST_EQUAL(U12, Z12, "ytoz: Y12");
	TEST_EQUAL(U21, Z21, "ytoz: Y21");
	TEST_EQUAL(U22, Z22, "ytoz: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytoh(y, u);
	TEST_EQUAL(U11, H11, "ytoh: H11");
	TEST_EQUAL(U12, H12, "ytoh: H12");
	TEST_EQUAL(U21, H21, "ytoh: H21");
	TEST_EQUAL(U22, H22, "ytoh: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytog(y, u);
	TEST_EQUAL(U11, G11, "ytog: G11");
	TEST_EQUAL(U12, G12, "ytog: G12");
	TEST_EQUAL(U21, G21, "ytog: G21");
	TEST_EQUAL(U22, G22, "ytog: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytoa(y, u);
	TEST_EQUAL(U11, A11, "ytoa: A11");
	TEST_EQUAL(U12, A12, "ytoa: A12");
	TEST_EQUAL(U21, A21, "ytoa: A21");
	TEST_EQUAL(U22, A22, "ytoa: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_ytob(y, u);
	TEST_EQUAL(U11, B11, "ytob: B11");
	TEST_EQUAL(U12, B12, "ytob: B12");
	TEST_EQUAL(U21, B21, "ytob: B21");
	TEST_EQUAL(U22, B22, "ytob: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_ytozi(y, v, z0);
	TEST_EQUAL(v[0], zi[0], "ytozi: zi0");
	TEST_EQUAL(v[1], zi[1], "ytozi: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_ytozin(*y, v, z0, 2);
	TEST_EQUAL(v[0], zi[0], "ytozin: zi0");
	TEST_EQUAL(v[1], zi[1], "ytozin: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htos(h, u, z0);
	TEST_EQUAL(U11, S11, "htos: S11");
	TEST_EQUAL(U12, S12, "htos: S12");
	TEST_EQUAL(U21, S21, "htos: S21");
	TEST_EQUAL(U22, S22, "htos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htot(h, u, z0);
	TEST_EQUAL(U11, T11, "htot: T11");
	TEST_EQUAL(U12, T12, "htot: T12");
	TEST_EQUAL(U21, T21, "htot: T21");
	TEST_EQUAL(U22, T22, "htot: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htoz(h, u);
	TEST_EQUAL(U11, Z11, "htoz: Y11");
	TEST_EQUAL(U12, Z12, "htoz: Y12");
	TEST_EQUAL(U21, Z21, "htoz: Y21");
	TEST_EQUAL(U22, Z22, "htoz: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htoy(h, u);
	TEST_EQUAL(U11, Y11, "htoy: Y11");
	TEST_EQUAL(U12, Y12, "htoy: Y12");
	TEST_EQUAL(U21, Y21, "htoy: Y21");
	TEST_EQUAL(U22, Y22, "htoy: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htog(h, u);
	TEST_EQUAL(U11, G11, "htog: G11");
	TEST_EQUAL(U12, G12, "htog: G12");
	TEST_EQUAL(U21, G21, "htog: G21");
	TEST_EQUAL(U22, G22, "htog: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htoa(h, u);
	TEST_EQUAL(U11, A11, "htoa: A11");
	TEST_EQUAL(U12, A12, "htoa: A12");
	TEST_EQUAL(U21, A21, "htoa: A21");
	TEST_EQUAL(U22, A22, "htoa: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_htob(h, u);
	TEST_EQUAL(U11, B11, "htob: B11");
	TEST_EQUAL(U12, B12, "htob: B12");
	TEST_EQUAL(U21, B21, "htob: B21");
	TEST_EQUAL(U22, B22, "htob: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_htozi(h, v, z0);
	TEST_EQUAL(v[0], zi[0], "htozi: zi0");
	TEST_EQUAL(v[1], zi[1], "htozi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtos(g, u, z0);
	TEST_EQUAL(U11, S11, "gtos: S11");
	TEST_EQUAL(U12, S12, "gtos: S12");
	TEST_EQUAL(U21, S21, "gtos: S21");
	TEST_EQUAL(U22, S22, "gtos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtot(g, u, z0);
	TEST_EQUAL(U11, T11, "gtot: T11");
	TEST_EQUAL(U12, T12, "gtot: T12");
	TEST_EQUAL(U21, T21, "gtot: T21");
	TEST_EQUAL(U22, T22, "gtot: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtoz(g, u);
	TEST_EQUAL(U11, Z11, "gtoz: Z11");
	TEST_EQUAL(U12, Z12, "gtoz: Z12");
	TEST_EQUAL(U21, Z21, "gtoz: Z21");
	TEST_EQUAL(U22, Z22, "gtoz: Z22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtoy(g, u);
	TEST_EQUAL(U11, Y11, "gtoy: Y11");
	TEST_EQUAL(U12, Y12, "gtoy: Y12");
	TEST_EQUAL(U21, Y21, "gtoy: Y21");
	TEST_EQUAL(U22, Y22, "gtoy: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtoh(g, u);
	TEST_EQUAL(U11, H11, "gtoh: H11");
	TEST_EQUAL(U12, H12, "gtoh: H12");
	TEST_EQUAL(U21, H21, "gtoh: H21");
	TEST_EQUAL(U22, H22, "gtoh: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtoa(g, u);
	TEST_EQUAL(U11, A11, "gtoa: A11");
	TEST_EQUAL(U12, A12, "gtoa: A12");
	TEST_EQUAL(U21, A21, "gtoa: A21");
	TEST_EQUAL(U22, A22, "gtoa: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_gtob(g, u);
	TEST_EQUAL(U11, B11, "gtob: B11");
	TEST_EQUAL(U12, B12, "gtob: B12");
	TEST_EQUAL(U21, B21, "gtob: B21");
	TEST_EQUAL(U22, B22, "gtob: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_gtozi(g, v, z0);
	TEST_EQUAL(v[0], zi[0], "gtozi: zi0");
	TEST_EQUAL(v[1], zi[1], "gtozi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atos(a, u, z0);
	TEST_EQUAL(U11, S11, "atos: S11");
	TEST_EQUAL(U12, S12, "atos: S12");
	TEST_EQUAL(U21, S21, "atos: S21");
	TEST_EQUAL(U22, S22, "atos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atot(a, u, z0);
	TEST_EQUAL(U11, T11, "atot: T11");
	TEST_EQUAL(U12, T12, "atot: T12");
	TEST_EQUAL(U21, T21, "atot: T21");
	TEST_EQUAL(U22, T22, "atot: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atoz(a, u);
	TEST_EQUAL(U11, Z11, "atoz: Y11");
	TEST_EQUAL(U12, Z12, "atoz: Y12");
	TEST_EQUAL(U21, Z21, "atoz: Y21");
	TEST_EQUAL(U22, Z22, "atoz: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atoy(a, u);
	TEST_EQUAL(U11, Y11, "atoy: Y11");
	TEST_EQUAL(U12, Y12, "atoy: Y12");
	TEST_EQUAL(U21, Y21, "atoy: Y21");
	TEST_EQUAL(U22, Y22, "atoy: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atoh(a, u);
	TEST_EQUAL(U11, H11, "atoh: H11");
	TEST_EQUAL(U12, H12, "atoh: H12");
	TEST_EQUAL(U21, H21, "atoh: H21");
	TEST_EQUAL(U22, H22, "atoh: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atog(a, u);
	TEST_EQUAL(U11, G11, "atog: G11");
	TEST_EQUAL(U12, G12, "atog: G12");
	TEST_EQUAL(U21, G21, "atog: G21");
	TEST_EQUAL(U22, G22, "atog: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_atob(a, u);
	TEST_EQUAL(U11, B11, "atob: B11");
	TEST_EQUAL(U12, B12, "atob: B12");
	TEST_EQUAL(U21, B21, "atob: B21");
	TEST_EQUAL(U22, B22, "atob: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_atozi(a, v, z0);
	TEST_EQUAL(v[0], zi[0], "atozi: zi0");
	TEST_EQUAL(v[1], zi[1], "atozi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btos(b, u, z0);
	TEST_EQUAL(U11, S11, "btos: S11");
	TEST_EQUAL(U12, S12, "btos: S12");
	TEST_EQUAL(U21, S21, "btos: S21");
	TEST_EQUAL(U22, S22, "btos: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btot(b, u, z0);
	TEST_EQUAL(U11, T11, "btot: T11");
	TEST_EQUAL(U12, T12, "btot: T12");
	TEST_EQUAL(U21, T21, "btot: T21");
	TEST_EQUAL(U22, T22, "btot: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btoz(b, u);
	TEST_EQUAL(U11, Z11, "btoz: Y11");
	TEST_EQUAL(U12, Z12, "btoz: Y12");
	TEST_EQUAL(U21, Z21, "btoz: Y21");
	TEST_EQUAL(U22, Z22, "btoz: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btoy(b, u);
	TEST_EQUAL(U11, Y11, "btoy: Y11");
	TEST_EQUAL(U12, Y12, "btoy: Y12");
	TEST_EQUAL(U21, Y21, "btoy: Y21");
	TEST_EQUAL(U22, Y22, "btoy: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btoh(b, u);
	TEST_EQUAL(U11, H11, "btoh: H11");
	TEST_EQUAL(U12, H12, "btoh: H12");
	TEST_EQUAL(U21, H21, "btoh: H21");
	TEST_EQUAL(U22, H22, "btoh: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btog(b, u);
	TEST_EQUAL(U11, G11, "btog: G11");
	TEST_EQUAL(U12, G12, "btog: G12");
	TEST_EQUAL(U21, G21, "btog: G21");
	TEST_EQUAL(U22, G22, "btog: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_btoa(b, u);
	TEST_EQUAL(U11, A11, "btoa: A11");
	TEST_EQUAL(U12, A12, "btoa: A12");
	TEST_EQUAL(U21, A21, "btoa: A21");
	TEST_EQUAL(U22, A22, "btoa: A22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_btozi(b, v, z0);
	TEST_EQUAL(v[0], zi[0], "btozi: zi0");
	TEST_EQUAL(v[1], zi[1], "btozi: zi1");

	if (opt_v) {
	    printf("-------------\n");
	}
    }
    result = T_PASS;

out:
    report_test_result("2x2 Conversions", result);
}


/*
 * test_conversions_3x3: test parameter conversions
 */
static void test_conversions_3x3()
{
    test_result_type result = T_SKIPPED;

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

	Z1  = crandn();
	Z2  = crandn();
	Z3  = crandn();
	z1c = conj(Z1);
	z2c = conj(Z2);
	z3c = conj(Z3);
	k1i = sqrt(fabs(creal(Z1)));
	k2i = sqrt(fabs(creal(Z2)));
	k3i = sqrt(fabs(creal(Z3)));
	a1  = crandn();
	a2  = crandn();
	a3  = crandn();
	S11 = crandn();
	S12 = crandn();
	S13 = crandn();
	S21 = crandn();
	S22 = crandn();
	S23 = crandn();
	S31 = crandn();
	S32 = crandn();
	S33 = crandn();
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
	    cmatrix_print("s", *s, 3, 3);
	}
	TEST_EQUAL(S11 * a1 + S12 * a2 + S13 * a3, b1, "S11,S12,S13");
	TEST_EQUAL(S21 * a1 + S22 * a2 + S23 * a3, b2, "S21,S22,S23");
	TEST_EQUAL(S31 * a1 + S32 * a2 + S33 * a3, b3, "S31,S32,S33");

	vnaconv_stozn(*s, *z, z0, 3);
	if (opt_v) {
	    cmatrix_print("z", *z, 3, 3);
	}
	TEST_EQUAL(Z11 * i1 + Z12 * i2 + Z13 * i3, v1, "stoz: Z11,Z12,Z13");
	TEST_EQUAL(Z21 * i1 + Z22 * i2 + Z23 * i3, v2, "stoz: Z21,Z22,Z23");
	TEST_EQUAL(Z31 * i1 + Z32 * i2 + Z33 * i3, v3, "stoz: Z31,Z32,Z33");

	vnaconv_stoyn(*s, *y, z0, 3);
	if (opt_v) {
	    cmatrix_print("y", *y, 3, 3);
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
	    cmatrix_print("zi", zi, 3, 1);
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
    report_test_result("3x3 Conversions", result);
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
    /*NOTREACHED*/
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
	    /*NOTREACHED*/
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
	/*NOTREACHED*/
    }
    test_conversions_2x2();
    test_conversions_3x3();

    exit(fail_count != 0);
    /*NOTREACHED*/
}
