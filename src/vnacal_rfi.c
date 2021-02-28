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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vnacal_internal.h"


#define EPS	1.0e-25

/*
 * Rational Function Interpolation
 *     Interpolate within an m-wide window of points using a ratio of
 *     polynomials, f(x) as the interpolation function.
 *
 *                 n0 + n1 x + n2 x^2 + ...
 *         f(x) = --------------------------
 *                  1 + d1 x + d2 x^2 + ...
 *
 *     For m odd, the order of both numerator and denominator is (m-1)/2.
 *     For m even, the denominator has order m/2 and numerator has order
 *     one less.
 */

/*
 * _vnacal_rfi: apply rational function interpolation
 *   @xp:      vector of x points
 *   @yp:      vector of y points
 *   @n:       length of xp and yp
 *   @m:       order (number of points that determine the interpolation)
 *   @segment: left x index that bounds x (used as hint on entry)
 *   @x:       dependent variable to interpolate
 */
double complex _vnacal_rfi(const double *xp, double complex *yp,
	int n, int m, int *segment, double x)
{
    int base;
    int nearest;
    int cur;
    double complex y;
    double complex c[m], d[m];

    assert(n >= 1);
    assert(m <= n);

    /*
     * Special-case extrapolation on the left.
     */
    if (x <= xp[0]) {
	return yp[0];
    }

    /*
     * Special-case extrapolation on the right.
     */
    if (x >= xp[n - 1]) {
	return yp[n - 1];
    }

    /*
     * Bound *segment to 0 .. n-1.
     */
    if (*segment < 0) {
	*segment = 0;
    } else if (*segment > n - 1) {
	*segment = n - 1;
    }

    /*
     * Using *segment as a hint, find the segment that bounds x.
     */
    if (x < xp[*segment]) {
	do {
	    --*segment;
	} while (x < xp[*segment]);
    } else if (x > xp[*segment + 1]) {
	do {
	    ++*segment;
	} while (x > xp[*segment + 1]);
    }

    /*
     * If x is equal to one of the bounds, return the associated y.
     * Otherise, find the xp index nearest x.
     */
    {
	double dx1, dx2;

	dx1 = fabs(x - xp[*segment]);
	if (dx1 <= EPS) {
	    return yp[*segment];
	}
	dx2 = fabs(x - xp[*segment + 1]);
	if (dx2 <= EPS) {
	    return yp[*segment + 1];
	}
	if (dx1 <= dx2 || m < 2) {
	    nearest = *segment;
	} else {
	    nearest = *segment + 1;
	}
    }

    /*
     * Find the base index of the m-wide window best centered around x.
     */
    if (m & 1) {
	base = nearest - ((m - 1) / 2);
    } else {
	base = *segment - ((m / 2) - 1);
    }
    if (base < 0) {
	base = 0;
    } else if (base + m > n) {
	base = n - m;
    }

    /*
     * Compute the rational function interpolation of x using the
     * Burlirch-Stoer algorithm.  This code is based on that given in
     * "Numerical Recipes in C", Press, Teukolsky, Vetterling, Flannery -
     * Cambridge Univ. Press - 1992 - pp. 112-113.
     */
    cur = nearest - base;
    assert(base >= 0 && base <= n - m);
    assert(cur >= 0 && cur < m);
    for (int i = 0; i < m; ++i) {
	c[i] = yp[base + i];
	d[i] = yp[base + i] + EPS;
    }
    y = yp[base + cur--];
    for (int i = 0; i < m - 1; ++i) {
	int j;

	for (j = 0; j < m - i - 1; ++j) {
	    double complex c_d = c[j + 1] - d[j];
	    double complex dx1 = x - xp[base + j];
	    double complex dx2 = x - xp[base + i + j + 1];
	    double complex den = dx1 * d[j] - dx2 * c[j + 1];

	    c[j] = c_d * dx1 * d[j]     / den;
	    d[j] = c_d * dx2 * c[j + 1] / den;
	}
	if (2 * (cur + 1) < m - i) {
	    assert(cur + 1 >= 0 && cur + 1 < m - i);
	    y += c[cur + 1];
	} else {
	    assert(cur >= 0 && cur < m - i);
	    y += d[cur--];
	}
    }
    return y;
}
