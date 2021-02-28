/*
 * Electrical Network Parameter Conversion Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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

    for (int trial = 0; trial < 100000; ++trial) {
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

	vnaconv_s2t(s, t);
	if (opt_v) {
	    cmatrix_print("t", *t, 2, 2);
	}
	TEST_EQUAL(T11 * a2 + T12 * b2, b1, "s2t: T11,T12");
	TEST_EQUAL(T21 * a2 + T22 * b2, a1, "s2t: T22,T22");

	vnaconv_s2z(s, z, z0);
	if (opt_v) {
	    cmatrix_print("z", *z, 2, 2);
	}
	TEST_EQUAL(Z11 * i1 + Z12 * i2, v1, "s2z: Z11,Z12");
	TEST_EQUAL(Z21 * i1 + Z22 * i2, v2, "s2z: Z21,Z22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_s2zn(*s, *u, z0, 2);
	TEST_EQUAL(U11 * i1 + U12 * i2, v1, "s2zn: U11,U12");
	TEST_EQUAL(U21 * i1 + U22 * i2, v2, "s2zn: U21,U22");

	vnaconv_s2y(s, y, z0);
	if (opt_v) {
	    cmatrix_print("y", *y, 2, 2);
	}
	TEST_EQUAL(Y11 * v1 + Y12 * v2, i1, "s2y: Y11,Y12");
	TEST_EQUAL(Y21 * v1 + Y22 * v2, i2, "s2y: Y21,Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_s2yn(*s, *u, z0, 2);
	TEST_EQUAL(U11 * v1 + U12 * v2, i1, "s2yn: U11,U12");
	TEST_EQUAL(U21 * v1 + U22 * v2, i2, "s2yn: U21,U22");

	vnaconv_s2h(s, h, z0);
	if (opt_v) {
	    cmatrix_print("h", *h, 2, 2);
	}
	TEST_EQUAL(H11 * i1 + H12 * v2, v1, "s2h: H11,H12");
	TEST_EQUAL(H21 * i1 + H22 * v2, i2, "s2h: H21,H22");

	vnaconv_s2g(s, g, z0);
	if (opt_v) {
	    cmatrix_print("g", *g, 2, 2);
	}
	TEST_EQUAL(G11 * v1 + G12 * i2, i1, "s2g: G11,G12");
	TEST_EQUAL(G21 * v1 + G22 * i2, v2, "s2g: G21,G22");

	vnaconv_s2a(s, a, z0);
	if (opt_v) {
	    cmatrix_print("a", *a, 2, 2);
	}
	TEST_EQUAL(A11 * v2 + A12 * -i2, v1, "s2a: A11,A12");
	TEST_EQUAL(A21 * v2 + A22 * -i2, i1, "s2b: A21,A22");

	vnaconv_s2b(s, b, z0);
	if (opt_v) {
	    cmatrix_print("b", *b, 2, 2);
	}
	TEST_EQUAL(B11 * v1 + B12 * i1,  v2, "s2b: B11,B12");
	TEST_EQUAL(B21 * v1 + B22 * i1, -i2, "s2b: B21,B22");

	vnaconv_s2zi(s, zi, z0);
	if (opt_v) {
	    cmatrix_print("zi", zi, 2, 1);
	}

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_s2zin(*s, v, z0, 2);
	TEST_EQUAL(v[0], zi[0], "s2zin: zi0");
	TEST_EQUAL(v[1], zi[1], "s2zin: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_s2zimn(*s, v, z0, 2, 2);
	TEST_EQUAL(v[0], zi[0], "s2zimn: zi0");
	TEST_EQUAL(v[1], zi[1], "s2zimn: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2s(t, u);
	TEST_EQUAL(U11, S11, "t2s: S11");
	TEST_EQUAL(U12, S12, "t2s: S12");
	TEST_EQUAL(U21, S21, "t2s: S21");
	TEST_EQUAL(U22, S22, "t2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2z(t, u, z0);
	TEST_EQUAL(U11, Z11, "t2z: Z11");
	TEST_EQUAL(U12, Z12, "t2z: Z12");
	TEST_EQUAL(U21, Z21, "t2z: Z21");
	TEST_EQUAL(U22, Z22, "t2z: Z22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2y(t, u, z0);
	TEST_EQUAL(U11, Y11, "t2y: Y11");
	TEST_EQUAL(U12, Y12, "t2y: Y12");
	TEST_EQUAL(U21, Y21, "t2y: Y21");
	TEST_EQUAL(U22, Y22, "t2y: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2h(t, u, z0);
	TEST_EQUAL(U11, H11, "t2h: H11");
	TEST_EQUAL(U12, H12, "t2h: H12");
	TEST_EQUAL(U21, H21, "t2h: H21");
	TEST_EQUAL(U22, H22, "t2h: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2g(t, u, z0);
	TEST_EQUAL(U11, G11, "t2g: G11");
	TEST_EQUAL(U12, G12, "t2g: G12");
	TEST_EQUAL(U21, G21, "t2g: G21");
	TEST_EQUAL(U22, G22, "t2g: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2a(t, u, z0);
	TEST_EQUAL(U11, A11, "t2a: A11");
	TEST_EQUAL(U12, A12, "t2a: A12");
	TEST_EQUAL(U21, A21, "t2a: A21");
	TEST_EQUAL(U22, A22, "t2a: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_t2b(t, u, z0);
	TEST_EQUAL(U11, B11, "t2b: B11");
	TEST_EQUAL(U12, B12, "t2b: B12");
	TEST_EQUAL(U21, B21, "t2b: B21");
	TEST_EQUAL(U22, B22, "t2b: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_t2zi(t, v, z0);
	TEST_EQUAL(v[0], zi[0], "t2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "t2zi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2s(z, u, z0);
	TEST_EQUAL(U11, S11, "z2s: S11");
	TEST_EQUAL(U12, S12, "z2s: S12");
	TEST_EQUAL(U21, S21, "z2s: S21");
	TEST_EQUAL(U22, S22, "z2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2sn(*z, *u, z0, 2);
	TEST_EQUAL(U11, S11, "z2sn: S11");
	TEST_EQUAL(U12, S12, "z2sn: S12");
	TEST_EQUAL(U21, S21, "z2sn: S21");
	TEST_EQUAL(U22, S22, "z2sn: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2t(z, u, z0);
	TEST_EQUAL(U11, T11, "z2t: T11");
	TEST_EQUAL(U12, T12, "z2t: T12");
	TEST_EQUAL(U21, T21, "z2t: T21");
	TEST_EQUAL(U22, T22, "z2t: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2y(z, u);
	TEST_EQUAL(U11, Y11, "z2y: Y11");
	TEST_EQUAL(U12, Y12, "z2y: Y12");
	TEST_EQUAL(U21, Y21, "z2y: Y21");
	TEST_EQUAL(U22, Y22, "z2y: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2yn(*z, *u, 2);
	TEST_EQUAL(U11, Y11, "z2yn: Y11");
	TEST_EQUAL(U12, Y12, "z2yn: Y12");
	TEST_EQUAL(U21, Y21, "z2yn: Y21");
	TEST_EQUAL(U22, Y22, "z2yn: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2h(z, u);
	TEST_EQUAL(U11, H11, "z2h: H11");
	TEST_EQUAL(U12, H12, "z2h: H12");
	TEST_EQUAL(U21, H21, "z2h: H21");
	TEST_EQUAL(U22, H22, "z2h: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2g(z, u);
	TEST_EQUAL(U11, G11, "z2g: G11");
	TEST_EQUAL(U12, G12, "z2g: G12");
	TEST_EQUAL(U21, G21, "z2g: G21");
	TEST_EQUAL(U22, G22, "z2g: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2a(z, u);
	TEST_EQUAL(U11, A11, "z2a: A11");
	TEST_EQUAL(U12, A12, "z2a: A12");
	TEST_EQUAL(U21, A21, "z2a: A21");
	TEST_EQUAL(U22, A22, "z2a: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2b(z, u);
	TEST_EQUAL(U11, B11, "z2b: B11");
	TEST_EQUAL(U12, B12, "z2b: B12");
	TEST_EQUAL(U21, B21, "z2b: B21");
	TEST_EQUAL(U22, B22, "z2b: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_z2zi(z, v, z0);
	TEST_EQUAL(v[0], zi[0], "z2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "z2zi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2s(y, u, z0);
	TEST_EQUAL(U11, S11, "y2s: S11");
	TEST_EQUAL(U12, S12, "y2s: S12");
	TEST_EQUAL(U21, S21, "y2s: S21");
	TEST_EQUAL(U22, S22, "y2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2t(y, u, z0);
	TEST_EQUAL(U11, T11, "y2t: T11");
	TEST_EQUAL(U12, T12, "y2t: T12");
	TEST_EQUAL(U21, T21, "y2t: T21");
	TEST_EQUAL(U22, T22, "y2t: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2z(y, u);
	TEST_EQUAL(U11, Z11, "y2z: Y11");
	TEST_EQUAL(U12, Z12, "y2z: Y12");
	TEST_EQUAL(U21, Z21, "y2z: Y21");
	TEST_EQUAL(U22, Z22, "y2z: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2h(y, u);
	TEST_EQUAL(U11, H11, "y2h: H11");
	TEST_EQUAL(U12, H12, "y2h: H12");
	TEST_EQUAL(U21, H21, "y2h: H21");
	TEST_EQUAL(U22, H22, "y2h: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2g(y, u);
	TEST_EQUAL(U11, G11, "y2g: G11");
	TEST_EQUAL(U12, G12, "y2g: G12");
	TEST_EQUAL(U21, G21, "y2g: G21");
	TEST_EQUAL(U22, G22, "y2g: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2a(y, u);
	TEST_EQUAL(U11, A11, "y2a: A11");
	TEST_EQUAL(U12, A12, "y2a: A12");
	TEST_EQUAL(U21, A21, "y2a: A21");
	TEST_EQUAL(U22, A22, "y2a: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2b(y, u);
	TEST_EQUAL(U11, B11, "y2b: B11");
	TEST_EQUAL(U12, B12, "y2b: B12");
	TEST_EQUAL(U21, B21, "y2b: B21");
	TEST_EQUAL(U22, B22, "y2b: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_y2zi(y, v, z0);
	TEST_EQUAL(v[0], zi[0], "y2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "y2zi: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_y2zin(*y, v, z0, 2);
	TEST_EQUAL(v[0], zi[0], "y2zin: zi0");
	TEST_EQUAL(v[1], zi[1], "y2zin: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2s(h, u, z0);
	TEST_EQUAL(U11, S11, "h2s: S11");
	TEST_EQUAL(U12, S12, "h2s: S12");
	TEST_EQUAL(U21, S21, "h2s: S21");
	TEST_EQUAL(U22, S22, "h2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2t(h, u, z0);
	TEST_EQUAL(U11, T11, "h2t: T11");
	TEST_EQUAL(U12, T12, "h2t: T12");
	TEST_EQUAL(U21, T21, "h2t: T21");
	TEST_EQUAL(U22, T22, "h2t: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2z(h, u);
	TEST_EQUAL(U11, Z11, "h2z: Y11");
	TEST_EQUAL(U12, Z12, "h2z: Y12");
	TEST_EQUAL(U21, Z21, "h2z: Y21");
	TEST_EQUAL(U22, Z22, "h2z: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2y(h, u);
	TEST_EQUAL(U11, Y11, "h2y: Y11");
	TEST_EQUAL(U12, Y12, "h2y: Y12");
	TEST_EQUAL(U21, Y21, "h2y: Y21");
	TEST_EQUAL(U22, Y22, "h2y: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2g(h, u);
	TEST_EQUAL(U11, G11, "h2g: G11");
	TEST_EQUAL(U12, G12, "h2g: G12");
	TEST_EQUAL(U21, G21, "h2g: G21");
	TEST_EQUAL(U22, G22, "h2g: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2a(h, u);
	TEST_EQUAL(U11, A11, "h2a: A11");
	TEST_EQUAL(U12, A12, "h2a: A12");
	TEST_EQUAL(U21, A21, "h2a: A21");
	TEST_EQUAL(U22, A22, "h2a: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_h2b(h, u);
	TEST_EQUAL(U11, B11, "h2b: B11");
	TEST_EQUAL(U12, B12, "h2b: B12");
	TEST_EQUAL(U21, B21, "h2b: B21");
	TEST_EQUAL(U22, B22, "h2b: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_h2zi(h, v, z0);
	TEST_EQUAL(v[0], zi[0], "h2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "h2zi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2s(g, u, z0);
	TEST_EQUAL(U11, S11, "g2s: S11");
	TEST_EQUAL(U12, S12, "g2s: S12");
	TEST_EQUAL(U21, S21, "g2s: S21");
	TEST_EQUAL(U22, S22, "g2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2t(g, u, z0);
	TEST_EQUAL(U11, T11, "g2t: T11");
	TEST_EQUAL(U12, T12, "g2t: T12");
	TEST_EQUAL(U21, T21, "g2t: T21");
	TEST_EQUAL(U22, T22, "g2t: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2z(g, u);
	TEST_EQUAL(U11, Z11, "g2z: Z11");
	TEST_EQUAL(U12, Z12, "g2z: Z12");
	TEST_EQUAL(U21, Z21, "g2z: Z21");
	TEST_EQUAL(U22, Z22, "g2z: Z22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2y(g, u);
	TEST_EQUAL(U11, Y11, "g2y: Y11");
	TEST_EQUAL(U12, Y12, "g2y: Y12");
	TEST_EQUAL(U21, Y21, "g2y: Y21");
	TEST_EQUAL(U22, Y22, "g2y: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2h(g, u);
	TEST_EQUAL(U11, H11, "g2h: H11");
	TEST_EQUAL(U12, H12, "g2h: H12");
	TEST_EQUAL(U21, H21, "g2h: H21");
	TEST_EQUAL(U22, H22, "g2h: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2a(g, u);
	TEST_EQUAL(U11, A11, "g2a: A11");
	TEST_EQUAL(U12, A12, "g2a: A12");
	TEST_EQUAL(U21, A21, "g2a: A21");
	TEST_EQUAL(U22, A22, "g2a: A22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_g2b(g, u);
	TEST_EQUAL(U11, B11, "g2b: B11");
	TEST_EQUAL(U12, B12, "g2b: B12");
	TEST_EQUAL(U21, B21, "g2b: B21");
	TEST_EQUAL(U22, B22, "g2b: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_g2zi(g, v, z0);
	TEST_EQUAL(v[0], zi[0], "g2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "g2zi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2s(a, u, z0);
	TEST_EQUAL(U11, S11, "a2s: S11");
	TEST_EQUAL(U12, S12, "a2s: S12");
	TEST_EQUAL(U21, S21, "a2s: S21");
	TEST_EQUAL(U22, S22, "a2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2t(a, u, z0);
	TEST_EQUAL(U11, T11, "a2t: T11");
	TEST_EQUAL(U12, T12, "a2t: T12");
	TEST_EQUAL(U21, T21, "a2t: T21");
	TEST_EQUAL(U22, T22, "a2t: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2z(a, u);
	TEST_EQUAL(U11, Z11, "a2z: Y11");
	TEST_EQUAL(U12, Z12, "a2z: Y12");
	TEST_EQUAL(U21, Z21, "a2z: Y21");
	TEST_EQUAL(U22, Z22, "a2z: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2y(a, u);
	TEST_EQUAL(U11, Y11, "a2y: Y11");
	TEST_EQUAL(U12, Y12, "a2y: Y12");
	TEST_EQUAL(U21, Y21, "a2y: Y21");
	TEST_EQUAL(U22, Y22, "a2y: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2h(a, u);
	TEST_EQUAL(U11, H11, "a2h: H11");
	TEST_EQUAL(U12, H12, "a2h: H12");
	TEST_EQUAL(U21, H21, "a2h: H21");
	TEST_EQUAL(U22, H22, "a2h: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2g(a, u);
	TEST_EQUAL(U11, G11, "a2g: G11");
	TEST_EQUAL(U12, G12, "a2g: G12");
	TEST_EQUAL(U21, G21, "a2g: G21");
	TEST_EQUAL(U22, G22, "a2g: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_a2b(a, u);
	TEST_EQUAL(U11, B11, "a2b: B11");
	TEST_EQUAL(U12, B12, "a2b: B12");
	TEST_EQUAL(U21, B21, "a2b: B21");
	TEST_EQUAL(U22, B22, "a2b: B22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_a2zi(a, v, z0);
	TEST_EQUAL(v[0], zi[0], "a2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "a2zi: zi1");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2s(b, u, z0);
	TEST_EQUAL(U11, S11, "b2s: S11");
	TEST_EQUAL(U12, S12, "b2s: S12");
	TEST_EQUAL(U21, S21, "b2s: S21");
	TEST_EQUAL(U22, S22, "b2s: S22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2t(b, u, z0);
	TEST_EQUAL(U11, T11, "b2t: T11");
	TEST_EQUAL(U12, T12, "b2t: T12");
	TEST_EQUAL(U21, T21, "b2t: T21");
	TEST_EQUAL(U22, T22, "b2t: T22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2z(b, u);
	TEST_EQUAL(U11, Z11, "b2z: Y11");
	TEST_EQUAL(U12, Z12, "b2z: Y12");
	TEST_EQUAL(U21, Z21, "b2z: Y21");
	TEST_EQUAL(U22, Z22, "b2z: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2y(b, u);
	TEST_EQUAL(U11, Y11, "b2y: Y11");
	TEST_EQUAL(U12, Y12, "b2y: Y12");
	TEST_EQUAL(U21, Y21, "b2y: Y21");
	TEST_EQUAL(U22, Y22, "b2y: Y22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2h(b, u);
	TEST_EQUAL(U11, H11, "b2h: H11");
	TEST_EQUAL(U12, H12, "b2h: H12");
	TEST_EQUAL(U21, H21, "b2h: H21");
	TEST_EQUAL(U22, H22, "b2h: H22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2g(b, u);
	TEST_EQUAL(U11, G11, "b2g: G11");
	TEST_EQUAL(U12, G12, "b2g: G12");
	TEST_EQUAL(U21, G21, "b2g: G21");
	TEST_EQUAL(U22, G22, "b2g: G22");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_b2a(b, u);
	TEST_EQUAL(U11, A11, "b2a: A11");
	TEST_EQUAL(U12, A12, "b2a: A12");
	TEST_EQUAL(U21, A21, "b2a: A21");
	TEST_EQUAL(U22, A22, "b2a: A22");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_b2zi(b, v, z0);
	TEST_EQUAL(v[0], zi[0], "b2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "b2zi: zi1");

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

    for (int trial = 0; trial < 100000; ++trial) {
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

	vnaconv_s2zn(*s, *z, z0, 3);
	if (opt_v) {
	    cmatrix_print("z", *z, 3, 3);
	}
	TEST_EQUAL(Z11 * i1 + Z12 * i2 + Z13 * i3, v1, "s2z: Z11,Z12,Z13");
	TEST_EQUAL(Z21 * i1 + Z22 * i2 + Z23 * i3, v2, "s2z: Z21,Z22,Z23");
	TEST_EQUAL(Z31 * i1 + Z32 * i2 + Z33 * i3, v3, "s2z: Z31,Z32,Z33");

	vnaconv_s2yn(*s, *y, z0, 3);
	if (opt_v) {
	    cmatrix_print("y", *y, 3, 3);
	}
	TEST_EQUAL(Y11 * v1 + Y12 * v2 + Y13 * v3, i1, "s2y: Y11,Y12,Y13");
	TEST_EQUAL(Y21 * v1 + Y22 * v2 + Y23 * v3, i2, "s2y: Y21,Y22,Y23");
	TEST_EQUAL(Y31 * v1 + Y32 * v2 + Y33 * v3, i3, "s2y: Y31,Y32,Y33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2sn(*z, *u, z0, 3);
	TEST_EQUAL(U11, S11, "z2s: S11");
	TEST_EQUAL(U12, S12, "z2s: S12");
	TEST_EQUAL(U13, S13, "z2s: S13");
	TEST_EQUAL(U21, S21, "z2s: S21");
	TEST_EQUAL(U22, S22, "z2s: S22");
	TEST_EQUAL(U23, S23, "z2s: S23");
	TEST_EQUAL(U31, S31, "z2s: S31");
	TEST_EQUAL(U32, S32, "z2s: S32");
	TEST_EQUAL(U33, S33, "z2s: S33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_z2yn(*z, *u, 3);
	TEST_EQUAL(U11, Y11, "z2y: Y11");
	TEST_EQUAL(U12, Y12, "z2y: Y12");
	TEST_EQUAL(U13, Y13, "z2y: Y13");
	TEST_EQUAL(U21, Y21, "z2y: Y21");
	TEST_EQUAL(U22, Y22, "z2y: Y22");
	TEST_EQUAL(U23, Y23, "z2y: Y23");
	TEST_EQUAL(U31, Y31, "z2y: Y31");
	TEST_EQUAL(U32, Y32, "z2y: Y32");
	TEST_EQUAL(U33, Y33, "z2y: Y33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2sn(*y, *u, z0, 3);
	TEST_EQUAL(U11, S11, "y2s: S11");
	TEST_EQUAL(U12, S12, "y2s: S12");
	TEST_EQUAL(U13, S13, "y2s: S13");
	TEST_EQUAL(U21, S21, "y2s: S21");
	TEST_EQUAL(U22, S22, "y2s: S22");
	TEST_EQUAL(U23, S23, "y2s: S23");
	TEST_EQUAL(U31, S31, "y2s: S31");
	TEST_EQUAL(U32, S32, "y2s: S32");
	TEST_EQUAL(U33, S33, "y2s: S33");

	(void)memset((void *)u, 0, sizeof(u));
	vnaconv_y2zn(*y, *u, 3);
	TEST_EQUAL(U11, Z11, "y2z: Z11");
	TEST_EQUAL(U12, Z12, "y2z: Z12");
	TEST_EQUAL(U13, Z13, "y2z: Z13");
	TEST_EQUAL(U21, Z21, "y2z: Z21");
	TEST_EQUAL(U22, Z22, "y2z: Z22");
	TEST_EQUAL(U23, Z23, "y2z: Z23");
	TEST_EQUAL(U31, Z31, "y2z: Z31");
	TEST_EQUAL(U32, Z32, "y2z: Z32");
	TEST_EQUAL(U33, Z33, "y2z: Z33");

	vnaconv_s2zin(*s, zi, z0, 3);
	if (opt_v) {
	    cmatrix_print("zi", zi, 3, 1);
	}

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_s2zimn(*s, v, z0, 3, 3);
	TEST_EQUAL(v[0], zi[0], "s2zimn: zi0");
	TEST_EQUAL(v[1], zi[1], "s2zimn: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_z2zin(*z, v, z0, 3);
	TEST_EQUAL(v[0], zi[0], "z2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "z2zi: zi1");
	TEST_EQUAL(v[2], zi[2], "z2zi: zi1");

	(void)memset((void *)v, 0, sizeof(v));
	vnaconv_y2zin(*y, v, z0, 3);
	TEST_EQUAL(v[0], zi[0], "y2zi: zi0");
	TEST_EQUAL(v[1], zi[1], "y2zi: zi1");
	TEST_EQUAL(v[2], zi[2], "y2zi: zi1");


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
