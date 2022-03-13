/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacommon_internal.h"


/*
 * _vnacommon_lu: replace A with its LU decomposition
 *  @a:         serialized input and output matrix
 *  @row_index: vector of n ints; on return maps original row to current row
 *  @n:         dimensions of A (which must be square)
 *
 *   L is placed below the major diagonal.  Its own major diagonal,
 *   which is not stored in the result matrix, is all ones.  U is placed
 *   in the upper diagonal and above.
 *
 * Returns the determinant of the matrix.
 */
double complex _vnacommon_lu(complex double *a, int *row_index, int n)
{
    double complex d = 1.0;	/* determinant */
    double row_scale[n];	/* 1 / largest magitude in each row */

#define A(i, j)		((a)[(i) * n + (j)])

    /*
     * Find row_scale.  Initialize row_index.
     */
    for (int i = 0; i < n; ++i) {
	double max = 0.0;

	for (int j = 0; j < n; ++j) {
	    double temp = cabs(A(i, j));

	    if (temp > max) {
		max = temp;
	    }
	}
	row_scale[i] = max;
	row_index[i] = i;
    }

    /*
     * Compute LU decomposition using Crout's method working column
     * by column.
     */
    for (int j = 0; j < n; ++j) {
	int    best_index = j;
	double best_value = 0.0;

	/*
	 * Compute U terms above the major diagonal.
	 */
	for (int i = 0; i < j; ++i) {
	    double complex s = 0.0;

	    s = A(i, j);
	    for (int k = 0; k < i; ++k) {
		s -= A(i, k) * A(k, j);
	    }
	    A(i, j) = s;
	}

	/*
	 * Compute the diagonal U term and L terms below the diagonal.
	 */
	for (int i = j; i < n; ++i) {
	    double complex s = A(i, j);
	    double temp;

	    for (int k = 0; k < j; ++k) {
		s -= A(i, k) * A(k, j);
	    }
	    A(i, j) = s;

	    /*
	     * Find the row with best value for the pivot position.
	     */
	    if ((temp = row_scale[i] * cabs(s)) > best_value) {
		best_index = i;
		best_value = temp;
	    }
	}

	/*
	 * Move the row with best pivot value into the pivot position.
	 */
	if (best_index != j) {
	    int itemp;

	    /*
	     * Swap rows.
	     */
	    for (int k = 0; k < n; ++k) {
		double complex ctemp;

		ctemp = A(best_index, k);
		A(best_index, k) = A(j, k);
		A(j, k) = ctemp;
	    }

	    /*
	     * Swap row indices.
	     */
	    itemp = row_index[best_index];
	    row_index[best_index] = row_index[j];
	    row_index[j] = itemp;

	    /*
	     * Move the row scale value down.
	     */
	    row_scale[best_index] = row_scale[j];

	    /*
	     * Negate the determinant.
	     */
	    d *= -1.0;
	}
	d *= A(j, j);

	/*
	 * Divide L terms by the pivot.
	 */
	if (j != n - 1) {
	    double complex scale = 1.0 / A(j, j);

	    for (int i = j + 1; i < n; ++i)
		A(i, j) *= scale;
	}
    }
    return d;
}
#undef A
