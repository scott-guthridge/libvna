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

#ifndef _TEST_H
#define _TEST_H

#include <complex.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * test_result_t
 */
typedef enum test_result {
    T_PASS	= 0,
    T_FAIL	= 1,
    T_SKIPPED	= 77,
    T_ERROR	= 99
} test_result_t;

/* progname: test program name */
extern char *progname;

/* test_isequal_eps: maximum allowed normalized error in test_isequal */
extern double test_isequal_eps;

/* report the result of the test to stdout */
extern void test_report(test_result_t result);

/* test_crandn: generate a 2d normally distributed random complex number */
extern double complex test_crandn();

/* test_crandn_nz: like test_crandn, except with magitude >= 0.1 */
double complex test_crandn_nz();

/* test_crandn_nrz: like test_crand_nz, but with angle in 20-160, 200-340 degrees */
extern double complex test_crandn_nrz();

/* test_is_equal: test if two values are equal */
extern bool test_isequal(double complex actual, double complex expected);

/* test_is_equal: test if two values are equal with label */
extern bool test_isequal_label(double complex actual, double complex expected,
	const char *label);

/* test_print_cmatrix: print an m by n serialized complex matrix */
extern void test_print_cmatrix(const char *tag, double complex *a,
	int m, int n);

/* test_init_isequal: init test_isequal_eps based on the machine precision */
extern void test_init_isequal();

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _TEST_H */
