/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2026 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "libt.h"
#include "vnadata.h"

#define SQRT2_2		0.707106781186547524400844362104

/*
 * libt_isequal_eps: maximum allowed normalized error in libt_isequal
 */
double libt_isequal_eps = 0.0;

/*
 * libt_isequal_init: initialize libt_isequal_eps based on machine precision
 */
void libt_isequal_init()
{
    double eps = 0.5;

    /*
     * Find the smallest number that when added to 1.0 compares
     * greater than one.
     */
    while (1.0 + 0.5 * eps > 1.0) {
	eps *= 0.5;
    }

    /*
     * Set eps to half the available precision minus one decimal place.
     */
    libt_isequal_eps = 10.0 * sqrt(eps);
}

/*
 * libt_isequal_d: test if two doubles are approximately equal
 *   @actual: actual value
 *   @expected: expected value
 */
bool libt_isequal_d(double actual, double expected)
{
    double error = fabs(actual - expected);
    double scale = fabs(expected);

    /*
     * If magnitude of the expected value is more than 1, normalize
     * the error to the expected value.
     */
    if (scale > 1.0) {
	error /= scale;
    }
    return error <= libt_isequal_eps;
}

/*
 * libt_isequal_c: test if two complex numbers are approximately equal
 *   @actual: actual value
 *   @expected: expected value
 */
bool libt_isequal_c(double complex actual, double complex expected)
{
    double error = cabs(actual - expected);
    double scale = cabs(expected);

    /*
     * If magnitude of the expected value is more than 1, normalize
     * the error to the expected value.
     */
    if (scale > 1.0) {
	error /= scale;
    }
    return error <= libt_isequal_eps;
}

/*
 * libt_isequal_d_rpt: test if two doubles are approximately equal
 *   @actual: actual value
 *   @expected: expected value
 */
bool libt_isequal_d_rpt(const char *prefix, double actual, double expected)
{
    double error = fabs(actual - expected);
    double scale = fabs(expected);

    /*
     * If magnitude of the expected value is more than 1, normalize
     * the error to the expected value.
     */
    if (scale > 1.0) {
	error /= scale;
    }
    if (error > libt_isequal_eps) {
	if (prefix != NULL) {
	    (void)printf("%s: ", prefix);
	}
	(void)printf("data miscompare: %f%+fj != %f%+fj (%f)",
		creal(actual), cimag(actual),
		creal(expected), cimag(expected),
		error);
	return false;
    }
    return true;
}

/*
 * libt_isequal_c_rpt: test if actual and expected are approximately equal
 */
bool libt_isequal_c_rpt(const char *prefix, double complex actual,
	double complex expected)
{
    double error = cabs(actual - expected);
    double scale = cabs(expected);

    /*
     * If magnitude of the expected value is more than 1, normalize
     * the error to the expected value.
     */
    if (scale > 1.0) {
	error /= scale;
    }
    if (error > libt_isequal_eps) {
	if (prefix != NULL) {
	    (void)printf("%s: ", prefix);
	}
	(void)printf("data miscompare: %f%+fj != %f%+fj (%f)",
		creal(actual), cimag(actual),
		creal(expected), cimag(expected),
		error);
	return false;
    }
    return true;
}

/*
 * libt_isequal_label: test if actual and expected are approximately equal
 */
bool libt_isequal_label(double complex actual, double complex expected,
	const char *label)
{
    if (!libt_isequal_c_rpt(label, actual, expected)) {
	(void)printf("\n");
	return false;
    }
    return true;
}

/*
 * libt_isequal: test if actual and expected are approximately equal
 */
bool libt_isequal(double complex actual, double complex expected)
{
    return libt_isequal_label(actual, expected, NULL);
}

/*
 * libt_randu: uniformally distributed numbers between min and max
 */
double libt_randu(double min, double max)
{
    return min + (max - min) * random() / (double)RANDOM_MAX;
}

/*
 * libt_randn2: return a pair of normally distributed random numbers
 */
double libt_randn2(double *second)
{
    double u1 = (random() + 1.0) / (RANDOM_MAX + 1.0); /* Box Muller method */
    double u2 = (double)random() / (RANDOM_MAX + 1.0);
    double r = sqrt(-2.0 * log(u1));
    double a = 2 * M_PI * u2;

    if (second != NULL) {
	*second = r * sin(a);
    }
    return r * cos(a);
}

/*
 * libt_randn: generate a normally distributed random number
 */
double libt_randn()
{
    return libt_randn2(NULL);
}

/*
 * libt_print_cmatrix: print an m by n serialized complex matrix
 */
void libt_print_cmatrix(const char *tag, const double complex *a, int m, int n)
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
 * libt_print_vnadata: print a vnadata structure
 */
void libt_print_vnadata(const char *tag, const vnadata_t *vdp)
{
    const int frequencies = vnadata_get_frequencies(vdp);
    const int rows = vnadata_get_rows(vdp);
    const int columns = vnadata_get_columns(vdp);
    const int ports = MAX(rows, columns);

    (void)printf("%s (%s) type %s:\n", tag, vnadata_get_name(vdp),
	    vnadata_get_type_name(vnadata_get_type(vdp)));
    if (!vnadata_has_fz0(vdp)) {
	printf("  z0:\n");
	for (int port = 0; port < ports; ++port) {
	    double complex z0 = vnadata_get_z0(vdp, port);

	    (void)printf("    %f%+fj\n", creal(z0), cimag(z0));
	}
	(void)printf("\n");
    } else {
	printf("  fz0:\n");
	for (int findex = 0; findex < frequencies; ++findex) {
	    (void)printf("    %2d:", findex);
	    for (int port = 0; port < ports; ++port) {
		double complex z0 = vnadata_get_fz0(vdp, findex, port);

		(void)printf(" %f%+fj", creal(z0), cimag(z0));
	    }
	    (void)printf("\n");
	}
	(void)printf("\n");
    }
    (void)printf("  data:\n");
    for (int findex = 0; findex < frequencies; ++findex) {
	(void)printf("    findex %d (%f)\n", findex,
		vnadata_get_frequency(vdp, findex));
	for (int row = 0; row < rows; ++row) {
	    (void)printf("   ");
	    for (int column = 0; column < columns; ++column) {
		double complex value = vnadata_get_cell(vdp, findex,
			row, column);
		(void)printf(" %f%+fj", creal(value), cimag(value));
	    }
	    (void)printf("\n");
	}
	(void)printf("\n");
    }
    (void)printf("\n");
}

/*
 * libt_error: report an error in the test itself and exit
 */
void libt_error(const char *format, ...)
{
    va_list ap;

    va_start(ap, format);
    (void)fprintf(stderr, "%s: ", progname);
    (void)vfprintf(stderr, format, ap);
    va_end(ap);

    exit(99);
}

/*
 * libt_fail: report a test failure and abort on opt_a
 */
void libt_fail(const char *format, ...)
{
    va_list ap;

    va_start(ap, format);
    (void)vprintf(format, ap);
    va_end(ap);
    if (opt_a) {
	abort();
    }
}

/* report the result of the test to stdout */
void libt_report(libt_result_t result)
{
    const char *result_string;

    switch (result) {
    case T_PASS:
	result_string = "PASS";
	break;
    case T_SKIPPED:
	result_string = "SKIPPED";
	break;
    case T_ERROR:
	result_string = "ERROR";
	break;
    default:
	result_string = "FAIL";
	break;
    }
    (void)printf("%s %s\n", progname, result_string);
}
