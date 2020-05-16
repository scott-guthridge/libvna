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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * interpolate: interpolate error terms for frequency
 *   @etsp: pointer to vnacal_etermset_t
 *   @etp: pointer to vnacal_error_terms_t
 *   @segment: pointer to current segment index
 *   @index: error term index
 *   @f: frequency
 */
static double complex interpolate(vnacal_etermset_t *etsp,
	vnacal_error_terms_t *etp, int *segment, int index, double f)
{
    assert(etsp->ets_frequencies >= 1);
    return _vnacal_rfi(etsp->ets_frequency_vector,
	    etp->et_data_vectors[index],
	    etsp->ets_frequencies,
	    MIN(etsp->ets_frequencies, VNACAL_MAX_M), segment, f);
}

/*
 * vnacal_apply: apply the calibration to measured values (simple form)
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 *   @frequencies: number of frequencies in matrix and s_parameters
 *   @frequency_vector: vector of increasing frequency points
 *   @matrix: vrows x vcolumns matrix of measured voltage ratios
 *   @s_parameters: caller-allocated vnadata_t structure
 */
int vnacal_apply(vnacal_t *vcp, int set, int frequencies,
	const double *frequency_vector, const double complex *const *matrix,
	vnadata_t *s_parameters)
{
    int rows, columns;
    double fmin, fmax;
    int segment = 0;
    vnacal_etermset_t *etsp;

    /*
     * Sanity check the arguments.
     */
    assert(vcp != NULL);
    if (set < 0 || set >= vcp->vc_sets) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply: invalid set index (%d)", set);
	return -1;
    }
    rows    = vnacal_get_rows(vcp, set);
    columns = vnacal_get_columns(vcp, set);
    etsp    = vcp->vc_set_vector[set];
    if (frequencies < 0) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply: %d: invalid frequency count", frequencies);
	return -1;
    }
    for (int i = 0; i < frequencies - 1; ++i) {
	if (frequency_vector[i] >= frequency_vector[i + 1]) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_apply: non-increasing frequencies");
	    return -1;
	}
    }
    fmin = _vnacal_etermset_get_fmin_bound(etsp);
    if (frequency_vector[0] < fmin) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply: frequency out of bounds %.3e < %.3e",
		frequency_vector[0], etsp->ets_frequency_vector[0]);
	return -1;
    }
    fmax = _vnacal_etermset_get_fmax_bound(etsp);
    if (frequency_vector[frequencies - 1] > fmax) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply: frequency out of bounds %.3e > %.3e",
		frequency_vector[frequencies - 1],
		etsp->ets_frequency_vector[etsp->ets_frequencies - 1]);
	return -1;
    }

    /*
     * Set up the output matrix.
     */
    if (vnadata_init(s_parameters, frequencies, rows, columns, VPT_S) == -1) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "vnadata_init: %s", strerror(errno));
	return -1;
    }

    /*
     * For each frequency index...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = etsp->ets_frequency_vector[findex];
	double complex a[columns * columns];
	double complex b[rows * columns];
	double complex s[rows * columns];
#define A(i, j)	(a[(i) * columns + (j)])
#define B(i, j)	(b[(i) * columns + (j)])
#define S(i, j)	(s[(i) * columns + (j)])
#define M(i, j) (matrix[(i) * columns + (j)][findex])

	/*
	 * Initialize A to the (columns x columns) identity matrix.
	 */
	for (int i = 0; i < columns; ++i) {
	    for (int j = 0; j < columns; ++j) {
		A(i, j) = (i == j) ? 1.0 : 0.0;
	    }
	}
	/*
	 * Compute the S parameters as follows.  First find B:
	 *
	 *   B = (M - E0) ./ E1
	 *
	 * where M is the matrix of measured voltage ratios, E0 is
	 * the matrix of directivity (diagonal terms) and leakage
	 * (off-diagonal terms) error terms, E1 is the matrix of
	 * reflection tracking (diagonal terms) and transmission tracking
	 * (off-diagonal terms) error terms, and the ./ operator is
	 * element by element division.  Then once we have B, find S by
	 * solving the matrix equation:
	 *
	 *   S (I + E2 .* B) = B
	 *
	 * for S, where I is the (columns x columns) identify matrix,
	 * E2 is the matrix of the port mismatch error terms, and the
	 * .* operator is element-by-element multiplication.  If the M
	 * and S matrices are rectangular, than the E2 .* B expression
	 * won't be the same size as I -- we fix that by truncating E2
	 * .* B if it's larger than I, or filling in missing cells with
	 * zeros if it's smaller.
	 */
	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		vnacal_error_terms_t *etp =
		    &etsp->ets_error_term_matrix[(row) * columns + (column)];
		double complex e0 = interpolate(etsp, etp, &segment, 0, f);
		double complex e1 = interpolate(etsp, etp, &segment, 1, f);
		double complex e2 = interpolate(etsp, etp, &segment, 2, f);
		double complex x = (M(row, column) - e0) / e1;

		if (row < columns) {
		    A(row, column) += e2 * x;
		}
		B(row, column) = x;
	    }
	}
	_vnacommon_mrdivide(s, b, a, rows, columns);
	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		(void)vnadata_set_cell(s_parameters, findex, row, column,
			S(row, column));
	    }
	}
#undef M
#undef S
#undef B
#undef A
    }
    return 0;
}
