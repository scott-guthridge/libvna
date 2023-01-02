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
#include <stdio.h>
#include <stdlib.h>
#include "vnacommon_internal.h"

/*
 * Make sure consecutive x-terms differ by at least this amount to
 * avoid an ill conditioned matrix.
 */
#define MIN_DX  0.0001

/*
 * Spline coefficient indices.
 */
enum {
    B = 0,
    C = 1,
    D = 2
};

/*
 * _vnacommon_spline_calc: find natural cubic spline coefficients
 *   @n:        number of spline segments
 *   @x_vector  n+1 element vector of x values
 *   @y_vector  n+1 element vector of y values
 *   @c_vector: n element vector of b,c,d tuples
 *
 * The resulting coefficients go into the following interpolation
 * polynomical:
 *   y(x) = y_vector[i] + c_vector[i][B] dx + c_vector[i][C] dx^2 +
 *                        c_vector[i][D] dx^3
 *   where:
 *     dx = x - x_vector[i],
 *     x_vector[i] <= x <= x_vector[i+1]
 */
int _vnacommon_spline_calc(int n, const double *x_vector,
	const double *y_vector, double (*c_vector)[3])
{
    int i;
    int rv = 0;
    double *hp = NULL;  /* widths of each segment */
    double *mp = NULL;  /* slopes of each segment */
    double *up = NULL;  /* major diagonal of matrix */
    double *vp = NULL;  /* right-hand terms of matrix */
    double *sp = NULL;  /* second derivative at x[i] */

    /*
     * Special-case a single element.
     */
    if (n < 2) {
	return 0;
    }

    /*
     * Allocate temporary vectors.
     */
    if ((hp = (double *)malloc(n * sizeof(double))) == NULL) {
	rv = -1;
	goto out;
    }
    if ((mp = (double *)malloc(n * sizeof(double))) == NULL) {
	rv = -1;
	goto out;
    }
    if ((up = (double *)malloc(n * sizeof(double))) == NULL) {
	rv = -1;
	goto out;
    }
    if ((vp = (double *)malloc(n * sizeof(double))) == NULL) {
	rv = -1;
	goto out;
    }
    if ((sp = (double *)malloc((n + 1) * sizeof(double))) == NULL) {
	rv = -1;
	goto out;
    }

    /*
     * Find the segment widths and slopes, and test that the x[i]
     * terms are increasing.
     */
    for (i = 0; i < n; ++i) {
	hp[i] =  x_vector[i+1] - x_vector[i];
	mp[i] = (y_vector[i+1] - y_vector[i]) / hp[i];
	if (hp[i] < MIN_DX) {
	    /* error reported by caller */
	    errno = EINVAL;
	    return -1;
	}
    }

    /*
     * Do gaussian elimination, taking advantage of the fact that
     * the matrix is sparse (tri-diagonal).
     */
    up[0] = 2.0 * (hp[0] + hp[1]);
    vp[0] = 6.0 * (mp[1] - mp[0]);
    for (i = 1; i < n - 1; ++i) {
	up[i] = 2.0 * (hp[i] + hp[i+1]) - hp[i] * hp[i]   / up[i-1];
	vp[i] = 6.0 * (mp[i+1] - mp[i]) - hp[i] * vp[i-1] / up[i-1];
    }

    /*
     * Do back-substitution to solve for the second derivative
     * terms.  The natural cubic spline is defined to have zero
     * second derivatives at the two endpoints.
     */
    sp[n] = 0.0;
    for (i = n-1; i >= 1; --i) {
	sp[i] = (vp[i-1] - hp[i] * sp[i+1]) / up[i-1];
    }
    sp[0] = 0.0;

    /*
     * Compute the coefficients.
     */
    for (i = 0; i < n; ++i) {
	c_vector[i][B] = (y_vector[i+1] - y_vector[i]) / hp[i]
	    - hp[i] / 3.0 * sp[i] - hp[i] / 6.0 * sp[i+1];
	c_vector[i][C] = sp[i] / 2.0;
	c_vector[i][D] = (sp[i+1] - sp[i]) / (6.0 * hp[i]);
    }
    rv = 0;

out:
    /*
     * Free the temporary vectors.
     */
    if (sp != NULL)
	free((void *)sp);
    if (vp != NULL)
	free((void *)vp);
    if (up != NULL)
	free((void *)up);
    if (mp != NULL)
	free((void *)mp);
    if (hp != NULL)
	free((void *)hp);

    return rv;
}

/*
 * _vnacommon_spline_eval: evaluate the spline at x
 *   @n:        number of spline segments
 *   @x_vector  n+1 element vector of x values
 *   @y_vector  n+1 element vector of y values
 *   @c_vector: n element vector of b,c,d tuples
 *   @x:        independent parameter
 */
double _vnacommon_spline_eval(int n, const double *x_vector,
	const double *y_vector, const double (*c_vector)[3], double x)
{
    int i;
    int low, high;
    double dx;

    /*
     * Validate n.
     */
    if (n < 1) {
	errno = EINVAL;
	return HUGE_VAL;
    }

    /*
     * Special-case one element.
     */
    if (n == 1) {
	return y_vector[0];
    }

    /*
     * Special case linear extrapolation on the left.
     */
    if (x < x_vector[0]) {
	double m; /* first derivative at x[0] */

	m  = c_vector[0][B];
	return m * (x - x_vector[0]) + y_vector[0];
    }

    /*
     * Special case linear extrapolation on the right.
     */
    if (x >= x_vector[n]) {
	double m;	/* first derivative at x[n] */

	dx = x_vector[n] - x_vector[n-1];
	m  = c_vector[n-1][B] + dx * (2.0 * c_vector[n-1][C] +
		dx * 3.0 * c_vector[n-1][D]);
	return m * (x - x_vector[n]) + y_vector[n];
    }

    /*
     * Binary search for the interval containing x and evaluate.
     */
    low  = 0;
    high = n-1;
    for (;;) {
	i = (low + high) / 2;
	if (low >= high) {
	    break;
	}
	if (x <  x_vector[i]) {
	    high = i - 1;
	    continue;
	}
	if (x >= x_vector[i+1]) {
	    low = i + 1;
	    continue;
	}
	assert(x >= x_vector[i]);
	assert(x <= x_vector[i+1]);
	break;
    }
    dx = x - x_vector[i];
    return y_vector[i] + dx * (c_vector[i][B] +
			       dx * (c_vector[i][C] +
				     dx * c_vector[i][D]));
}
