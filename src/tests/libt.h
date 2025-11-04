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

#ifndef _LIBT_H
#define _LIBT_H

#include <complex.h>
#include <stdarg.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * libt_result_t
 */
typedef enum libt_result {
    T_PASS	= 0,
    T_FAIL	= 1,
    T_SKIPPED	= 77,
    T_ERROR	= 99
} libt_result_t;

/* progname: test program name */
extern char *progname;

extern bool opt_a;
extern int opt_v;

/* libt_isequal_eps: maximum allowed normalized error in libt_isequal */
extern double libt_isequal_eps;

/* libt_isequal_init: init libt_isequal_eps based on the machine precision */
extern void libt_isequal_init();

/* libt_isequal_d: test if two doubles are approximately equal */
extern bool libt_isequal_d(double actual, double expected);

/* libt_isequal_d_rpt: test for equality and report miscompare if not */
extern bool libt_isequal_d_rpt(const char *prefix, double actual,
	double expected);

/* libt_isequal_c: test if two complex numbers are approximately equal */
extern bool libt_isequal_c(double complex actual, double complex expected);

extern bool libt_isequal_c_rpt(const char *prefix, double complex actual,
	double complex expected);

/* libt_isequal: test if two values are equal */
extern bool libt_isequal(double complex actual, double complex expected);

/* libt_isequal: test if two values are equal with label */
extern bool libt_isequal_label(double complex actual, double complex expected,
	const char *label);

/* libt_print_cmatrix: print an m by n serialized complex matrix */
extern void libt_print_cmatrix(const char *tag, const double complex *a,
	int m, int n);

/* libt_randu: uniformally distributed numbers between min and max */
extern double libt_randu(double min, double max);

/* libt_randn: return a normally distributed random number */
extern double libt_randn();

/* libt_randn: return a pair of normally distributed random numbers */
extern double libt_randn2(double *);

/* libt_rand_nsmm: return a truncated Rice(nu, sigma) random number */
extern double libt_rand_nsmm(double nu, double sigma, double min, double max);

/* libt_fail: report a test failure and abort if opt_a */
extern void libt_fail(const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 1, 2)));
#else
    ;
#endif

/* libt_error: report an error in the test itself and exit */
extern void libt_error(const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 1, 2)))
	__attribute__ ((__noreturn__));
#else
    ;
#endif

    /* report the result of the test to stdout */
    extern void libt_report(libt_result_t result);


#ifdef __cplusplus
    } /* extern "C" */
#endif

#endif /* _LIBT_H */
