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

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "libt.h"

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
     * Set eps to half the available precision.
     */
    libt_isequal_eps = sqrt(eps);
}

/*
 * libt_isequal: test if actual and expected are approximately equal
 */
bool libt_isequal(double complex actual, double complex expected)
{
    double d = cabs(expected);

    if (d < 1.0) {
	d = 1.0;
    }
    if (cabs(actual - expected) / d > libt_isequal_eps) {
	printf("data miscompare: %f%+fj != %f%+fj (%f)\n",
		creal(actual), cimag(actual), creal(expected), cimag(expected),
		cabs(actual - expected) / d);
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
    double d = cabs(expected);

    if (d < 1.0) {
	d = 1.0;
    }
    if (cabs(actual - expected) / d > libt_isequal_eps) {
	printf("%s: %f%+fj != %f%+fj (%f)\n", label,
		creal(actual), cimag(actual), creal(expected), cimag(expected),
		cabs(actual - expected) / d);
	return false;
    }
    return true;
}

/*
 * libt_crandn: generate a random complex number where real and imaginary parts
 *	are normally distributed with zero mean and unit standard deviation
 */
double complex libt_crandn()
{
    double u1 = (random() + 1.0) / RAND_MAX;	/* Box Muller method */
    double u2 = (double)random() / RAND_MAX;
    double r = sqrt(-2.0 * log(u1));
    double a = 2 * M_PI * u2;

    return r * (cos(a) + I * sin(a));
}

/*
 * libt_crandn_nz: like libt_crandn, except with magitude >= 0.1
 */
double complex libt_crandn_nz()
{
    double u1 = (random() + 1.0) / RAND_MAX;
    double u2 = (double)random() / RAND_MAX;
    double r = 0.1 + 0.9 * sqrt(-2.0 * log(u1));
    double a = 2 * M_PI * u2;

    return r * (cos(a) + I * sin(a));
}

/*
 * libt_crandn_nrz: like libt_crandn_nz, angle in 20-160, 200-340 degrees
 */
double complex libt_crandn_nrz()
{
    double u1 = (random() + 1.0) / RAND_MAX;
    double u2 = (double)random() / RAND_MAX;
    double r = 0.1 + 0.9 * sqrt(-2.0 * log(u1));
    double d = (2.0 * u2 - 1.0) * 140.0;
    double a;

    /*
     * Coming into this expression, d is +/- 140.  If non-negative,
     * shift it to 20 ... 160.  If negative, shift it to -160 ... -20
     */
    if (d >= 0.0) {
	d += 20.0;
    } else {
	d -= 20.0;
    }
    a = M_PI / 180.0 * d;

    return r * (cos(a) + I * sin(a));
}

/*
 * libt_print_cmatrix: print an m by n serialized complex matrix
 */
void libt_print_cmatrix(const char *tag, double complex *a, int m, int n)
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
