/*
 * Vector Network Analyzer Calibration Library
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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _vnacal_interpolate: interpolate error terms for frequency
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
 * vnacal_input_alloc: allocate an application of libvnacal
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 *   @rows: number of rows in the S-parameter matrix
 *   @columns: number of columns in the S-parameter matrix
 *   @frequencies: number of measured frequencies
 *   @frequency_vector: vector of measured frequencies
 *
 * The frequencies given in the input don't have to be the same as
 * those in the calibration; however, they may not extend outside of
 * the calibration range.  In other words, the library will interpolate,
 * but it will not extrapolate.
 */
vnacal_input_t *vnacal_input_alloc(vnacal_t *vcp, int set,
	int rows, int columns, int frequencies)
{
    vnacal_input_t *vip;

    if (vcp == NULL) {
	errno = EINVAL;
	return NULL;
    }
    if (set < 0 || set >= vcp->vc_sets) {
	_vnacal_error(vcp, "vnacal_input_alloc: "
		"invalid set index (%d)", set);
	errno = EINVAL;
	return NULL;
    }
    if (rows < 0 || columns < 0) {
	_vnacal_error(vcp, "vnacal_input_alloc: "
		"invalid dimension (%d x %d)", rows, columns);
	errno = EINVAL;
	return NULL;
    }
    if (frequencies < 0) {
	_vnacal_error(vcp, "vnacal_input_alloc: "
		"invalid frequency count (%d)", frequencies);
	errno = EINVAL;
	return NULL;
    }
    if ((vip = malloc(sizeof(vnacal_input_t))) == NULL) {
	_vnacal_error(vcp, "malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)vip, 0, sizeof(*vip));
    vip->vi_vcp		= vcp;
    vip->vi_set		= set;
    vip->vi_rows              = rows;
    vip->vi_columns           = columns;
    vip->vi_frequencies       = frequencies;
    vip->vi_frequencies_valid = false;
    if ((vip->vi_frequency_vector = calloc(frequencies,
		    sizeof(double))) == NULL) {
	_vnacal_error(vcp, "calloc: %s", strerror(errno));
	vnacal_input_free(vip);
	return NULL;
    }
    if ((vip->vi_matrix = calloc(rows * columns,
		    sizeof(double complex *))) == NULL) {
	_vnacal_error(vcp, "calloc: %s", strerror(errno));
	vnacal_input_free(vip);
	return NULL;
    }
    for (int cell = 0; cell < rows * columns; ++cell) {
	if ((vip->vi_matrix[cell] = calloc(frequencies,
			sizeof(double complex))) == NULL) {
	    _vnacal_error(vcp, "calloc: %s", strerror(errno));
	    vnacal_input_free(vip);
	    return NULL;
	}
    }
    if ((vip->vi_counts = calloc(rows * columns, sizeof(int))) == NULL) {
	_vnacal_error(vcp, "calloc: %s", strerror(errno));
	vnacal_input_free(vip);
	return NULL;
    }
    if ((vip->vi_map = calloc(rows * columns, sizeof(int))) == NULL) {
	_vnacal_error(vcp, "calloc: %s", strerror(errno));
	vnacal_input_free(vip);
	return NULL;
    }
    for (int i = 0; i < rows * columns; ++i) {
	vip->vi_map[i] = -1;
    }
    return vip;
}

/*
 * vnacal_input_set_frequency_vector: set the frequency vector
 *   @vip: pointer to opaque object returned from vnacal_input_alloc
 *   @frequency_vector: vector of measured frequencies
 */
int vnacal_input_set_frequency_vector(vnacal_input_t *vip,
	const double *frequency_vector)
{
    vnacal_t *vcp;
    vnacal_etermset_t *etsp;
    double fmin, fmax;

    if (vip == NULL) {
	errno = EINVAL;
	return -1;
    }
    vcp = vip->vi_vcp;
    etsp = vcp->vc_set_vector[vip->vi_set];
    if (frequency_vector == NULL) {
	_vnacal_error(vcp, "vnacal_input_set_frequency_vector: "
		"invalid NULL frequency vector");
	errno = EINVAL;
	return -1;
    }
    for (int i = 0; i < vip->vi_frequencies - 1; ++i) {
	if (frequency_vector[i] >= frequency_vector[i + 1]) {
	    errno = EINVAL;
	    return -1;
	}
    }
    fmin = _vnacal_etermset_get_fmin_bound(etsp);
    fmax = _vnacal_etermset_get_fmax_bound(etsp);
    if (frequency_vector[0] < fmin) {
	_vnacal_error(vcp, "vnacal_input_set_frequency_vector: "
		"frequency out of bounds %.3e < %.3e",
		frequency_vector[0], etsp->ets_frequency_vector[0]);
	errno = EINVAL;
	return -1;
    }
    if (frequency_vector[vip->vi_frequencies - 1] > fmax) {
	_vnacal_error(vcp, "vnacal_input_set_frequency_vector: "
		"frequency out of bounds %.3e > %.3e",
		frequency_vector[vip->vi_frequencies - 1],
		etsp->ets_frequency_vector[etsp->ets_frequencies - 1]);
	errno = EINVAL;
	return -1;
    }
    (void)memcpy((void *)vip->vi_frequency_vector, (void *)frequency_vector,
	    vip->vi_frequencies * sizeof(double));
    vip->vi_frequencies_valid = true;
    return 0;
}

/*
 * vnacal_input_add_vector: add a vector of measurements to the input
 *   @vip: pointer to opaque object returned from vnacal_input_alloc
 *   @row: measured DUT port (zero-based)
 *   @column: driven DUT port (zero-based)
 *   @vector: vector of voltage measurements, one per frequency
 *
 * Repeated calls to this function on the same row and column average
 * the vectors given.
 */
int vnacal_input_add_vector(vnacal_input_t *vip,
	int row, int column, const double complex *vector)
{
    vnacal_t *vcp;
    vnacal_etermset_t *etsp;
    int map;
    int cell;
    double complex *destination_vector;

    if (vip == NULL) {
	errno = EINVAL;
	return -1;
    }
    vcp = vip->vi_vcp;
    etsp = vcp->vc_set_vector[vip->vi_set];
    cell = row * vip->vi_columns + column;
    if (row < 0 || row >= vip->vi_rows) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid row: %d", row);
	errno = EINVAL;
	return -1;
    }
    if (column < 0 || column >= vip->vi_columns) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid column: %d", column);
	errno = EINVAL;
	return -1;
    }
    if (vector == NULL) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid NULL vector");
	errno = EINVAL;
	return -1;
    }
    if (etsp->ets_rows * etsp->ets_columns == 2) {
	map = (row == column) ? 0 : 1;
    } else if (row < etsp->ets_rows && column < etsp->ets_columns) {
	map = etsp->ets_columns * row + column;
    } else {
	_vnacal_error(vcp, "vnacal_input_add_vector: ambiguous DUT "
		"to VNA port map: use vnacal_input_add_vector_map "
		"instead");
	errno = EINVAL;
	return -1;
    }
    if (vip->vi_map[cell] != -1 && vip->vi_map[cell] != map) {
	_vnacal_error(vcp, "vnacal_input_add_vector: inconsistent "
		"DUT to VNA port mapping %d,%d -> %d,%d (previously %d,%d)",
		row, column, map / etsp->ets_columns, map % etsp->ets_columns,
		vip->vi_map[cell] / etsp->ets_columns,
		vip->vi_map[cell] % etsp->ets_columns);
	errno = EINVAL;
	return -1;
    }
    vip->vi_map[cell] = map;
    destination_vector = vip->vi_matrix[cell];
    for (int findex = 0; findex < vip->vi_frequencies; ++findex) {
	destination_vector[findex] += vector[findex];
    }
    ++vip->vi_counts[cell];
    return 0;
}

/*
 * vnacal_input_add_mapped_vector: add a vector of measurements to the input
 *   @vip: pointer to opaque object returned from vnacal_input_alloc
 *   @vrow: VNA detector port (zero-based)
 *   @vcolumn: VNA driving port (zero-based)
 *   @drow: measured DUT port (zero-based)
 *   @dcolumn: driven DUT port (zero-based)
 *   @vector: vector of voltage measurements, one per frequency
 *
 * Repeated calls to this function on the same row and column average
 * the vectors given.
 */
int vnacal_input_add_mapped_vector(vnacal_input_t *vip,
	int vrow, int vcolumn, int drow, int dcolumn,
	const double complex *vector)
{
    vnacal_t *vcp;
    vnacal_etermset_t *etsp;
    int map;
    int cell;
    double complex *destination_vector;

    if (vip == NULL) {
	errno = EINVAL;
	return -1;
    }
    vcp = vip->vi_vcp;
    etsp = vcp->vc_set_vector[vip->vi_set];
    cell = drow * vip->vi_columns + dcolumn;
    if (vrow < 0 || vrow >= etsp->ets_rows) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid vrow: %d", vrow);
	errno = EINVAL;
	return -1;
    }
    if (vcolumn < 0 || vcolumn >= etsp->ets_columns) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid vcolumn: %d", vcolumn);
	errno = EINVAL;
	return -1;
    }
    if (drow < 0 || drow >= vip->vi_rows) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid drow: %d", drow);
	errno = EINVAL;
	return -1;
    }
    if (dcolumn < 0 || dcolumn >= vip->vi_columns) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid dcolumn: %d", dcolumn);
	errno = EINVAL;
	return -1;
    }
    if (vector == NULL) {
	_vnacal_error(vcp, "vnacal_input_add_vector: "
	   "invalid NULL vector");
	errno = EINVAL;
	return -1;
    }
    map = etsp->ets_columns * vrow + vcolumn;
    if (vip->vi_map[cell] != -1 && vip->vi_map[cell] != map) {
	_vnacal_error(vcp, "vnacal_input_add_vector: inconsistent "
		"DUT to VNA port mapping %d,%d -> %d,%d (previously %d,%d)",
		drow, dcolumn, map / etsp->ets_columns, map % etsp->ets_columns,
		vip->vi_map[cell] / etsp->ets_columns,
		vip->vi_map[cell] % etsp->ets_columns);
	errno = EINVAL;
	return -1;
    }
    vip->vi_map[cell] = map;
    destination_vector = vip->vi_matrix[cell];
    for (int findex = 0; findex < vip->vi_frequencies; ++findex) {
	destination_vector[findex] += vector[findex];
    }
    ++vip->vi_counts[cell];
    return 0;
}

/*
 * GET_VALUE: get the given average value from the vnacal_input_t
 */
#define GET_VALUE(vip, s_cell, findex) \
	((vip)->vi_counts[s_cell] == 0 ? 0.0 : \
	 ((vip)->vi_matrix[s_cell][findex] / \
	 (double)(vip)->vi_counts[s_cell]))

/*
 * vnacal_input_get_value: get the given uncalibrated value
 *   @vip: pointer to opaque object returned from vnacal_input_alloc
 *   @row: measured DUT port (zero-based)
 *   @column: driven DUT port (zero-based)
 *   @findex: frequency index
 */
double complex vnacal_input_get_value(vnacal_input_t *vip,
	int row, int column, int findex)
{
    vnacal_t *vcp;
    int cell;

    /*
     * Validate arguments
     */
    if (vip == NULL) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    vcp = vip->vi_vcp;
    if (row < 0 || row >= vip->vi_rows) {
	_vnacal_error(vcp, "vnacal_input_get_value: "
	   "invalid row: %d", row);
	errno = EINVAL;
	return HUGE_VAL;
    }
    if (column < 0 || column >= vip->vi_columns) {
	_vnacal_error(vcp, "vnacal_input_get_value: "
	   "invalid column: %d", column);
	errno = EINVAL;
	return HUGE_VAL;
    }
    if (findex < 0 || findex >= vip->vi_frequencies) {
	_vnacal_error(vcp, "vnacal_input_get_value: "
	   "invalid findex: %d", findex);
	errno = EINVAL;
	return HUGE_VAL;
    }

    /*
     * Calculate the value.
     */
    cell = row * vip->vi_columns + column;
    return GET_VALUE(vip, cell, findex);
}

/*
 * vnacal_input_apply: apply the calibration and generate S-parameters
 *   @vip: pointer to opaque object returned from vnacal_input_alloc
 *   @s_parameters: caller-allocated vnadata_t object
 */
int vnacal_input_apply(const vnacal_input_t *vip,
	vnadata_t *s_parameters)
{
    vnacal_t *vcp;

    /*
     * Validate arguments
     */
    if (vip == NULL) {
	errno = EINVAL;
	goto error;
    }
    vcp = vip->vi_vcp;
    if (s_parameters == NULL) {
	_vnacal_error(vcp, "vnacal_input_apply: "
	   "invalid NULL s_parameters");
	errno = EINVAL;
	goto error;
    }
    if (!vip->vi_frequencies_valid) {
	_vnacal_error(vcp, "vnacal_input_apply: "
	   "no frequency vector given");
	errno = EINVAL;
	goto error;
    }

    /*
     * Calculate the S parameters
     */
    {
	int rows = vip->vi_rows;
	int columns = vip->vi_columns;
	int frequencies = vip->vi_frequencies;
	int segment = 0;
	double *frequency_vector = vip->vi_frequency_vector;
	vnacal_etermset_t *etsp = vcp->vc_set_vector[vip->vi_set];
	double complex a[columns * columns];
	double complex b[rows * columns];
	double complex s[rows * columns];
#define A(i, j)	(a[(i) * columns + (j)])
#define B(i, j)	(b[(i) * columns + (j)])
#define S(i, j)	(s[(i) * columns + (j)])

	/*
	 * Set up the output matrix.
	 */
	if (vnadata_init(s_parameters, frequencies, rows, columns,
		    VPT_S) == -1) {
	    _vnacal_error(vcp, "_vnacal_input_apply: %s",
		    strerror(errno));
	    goto error;
	}
	(void)vnadata_set_frequency_vector(s_parameters, frequency_vector);
	(void)vnadata_set_all_z0(s_parameters, etsp->ets_z0);

	/*
	 * Compute actual S-parameters.
	 */
	for (int findex = 0; findex < frequencies; ++findex) {
	    /*
	     * Initialize A to the identity matrix.
	     */
	    for (int i = 0; i < columns; ++i) {
		for (int j = 0; j < columns; ++j) {
		    A(i, j) = (i == j) ? 1.0 : 0.0;
		}
	    }

	    /*
	     * Form a (columns x columns) matrix A and (rows x columns)
	     * matrix B what we'll use to solve for the S parameter matrix.
	     *
	     * B represents the voltages emanating from each DUT port.
	     * We set A = I + E*B, where I is the (columns x columns)
	     * identity matrix, E is the matrix of e11(or e22) error
	     * terms, expanded out with zeros to be of dimension (columns
	     * x columns), B is the B matrix above, but also expanded
	     * out with zeros to match the dimensions, and operator '*'
	     * is element-by-element multiplication.  The matrix division
	     * S = B / A then produces S.
	     */
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    int s_cell = row * columns + column;/* cell of S matrix */
		    int c_cell;				/* cell of cal matrix */
		    vnacal_error_terms_t *etp;
		    double complex u, v;

		    /*
		     * Get the calibration matrix cell, row and column.  If no
		     * mapping has been filled in (because no input vector was
		     * supplied in this cell), default to a sane value:
		     * diagonal elements to 0 and off-diagonal to 1.
		     */
		    if ((c_cell = vip->vi_map[s_cell]) == -1) {
			c_cell = (row == column) ? 0 : 1;
		    }

		    /*
		     * Find the error terms and calculate a and b matrices.
		     */
		    etp = &etsp->ets_error_term_matrix[c_cell];
		    {
			double complex e00_e30 = interpolate(etsp, etp,
				&segment, 0, frequency_vector[findex]);
			double complex e10e01_e10e32 = interpolate(etsp, etp,
				&segment, 1, frequency_vector[findex]);
			double complex e11_e22 = interpolate(etsp, etp,
				&segment, 2, frequency_vector[findex]);
			double complex measured;

			measured = GET_VALUE(vip, s_cell, findex);
			v = (measured - e00_e30) / e10e01_e10e32;
			u = e11_e22 * v;

		    }
		    if (row < columns) {
			A(row, column) += u;
		    }
		    B(row, column) = v;
		}
	    }

	    /*
	     * Calculate the actual S parameters and copy to s_parameters.
	     */
	    _vnacommon_mrdivide(s, b, a, rows, columns);
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    if (vnadata_set_cell(s_parameters, findex, row, column,
			    S(row, column)) != 0) {
			_vnacal_error(vcp,
				"_vnacal_input_get_s_parameters: "
				"vnadata_set_cell: %s",
			    strerror(errno));
			goto error;
		    }
		}
	    }
	}
    }
    return 0;

error:
    return -1;
}
#undef A
#undef B
#undef F

/*
 * vnacal_input_free: free the input object
 *   @vip: pointer to opaque object returned from vnacal_input_alloc
 */
void vnacal_input_free(vnacal_input_t *vip)
{
    if (vip != NULL) {
	int rows = vip->vi_rows;
	int columns = vip->vi_columns;

	free((void *)vip->vi_map);
	free((void *)vip->vi_counts);
	if (vip->vi_matrix != NULL) {
	    for (int cell = 0; cell < rows * columns; ++cell) {
		free((void *)vip->vi_matrix[cell]);
	    }
	    free((void *)vip->vi_matrix);
	}
	free((void *)vip->vi_frequency_vector);
	free((void *)vip);
    }
}
