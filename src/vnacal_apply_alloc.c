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
 * vnacal_apply_alloc: allocate an application of libvnacal
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 *   @drows: number of rows in the DUT S-parameter matrix
 *   @dcolumns: number of columns in the DUT S-parameter matrix
 *   @frequencies: number of measured frequencies
 *
 * The frequencies given in the input don't have to be the same as
 * those in the calibration; however, they may not extend outside of
 * the calibration range.  In other words, the library will interpolate,
 * but it will not extrapolate.
 */
vnacal_apply_t *vnacal_apply_alloc(vnacal_t *vcp, int set,
	int drows, int dcolumns, int dfrequencies)
{
    vnacal_apply_t *vap;
    int dcells = drows * dcolumns;

    /*
     * Sanity check the arguments.
     */
    if (vcp == NULL) {
	errno = EINVAL;
	return NULL;
    }
    if (set < 0 || set >= vcp->vc_sets) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply_alloc: invalid set index (%d)", set);
	return NULL;
    }
    if (drows < 0 || dcolumns < 0) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply_alloc: invalid dimension (%d x %d)",
		drows, dcolumns);
	return NULL;
    }
    if (dfrequencies < 0) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply_alloc: invalid frequency count (%d)",
		dfrequencies);
	return NULL;
    }

    /*
     * Allocate the vnacal_apply_t structure.
     */
    if ((vap = malloc(sizeof(vnacal_apply_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)vap, 0, sizeof(*vap));
    vap->va_vcp			= vcp;
    vap->va_set			= set;
    vap->va_vrows		= vnacal_get_rows(vcp, set);
    vap->va_vcolumns		= vnacal_get_columns(vcp, set);
    vap->va_drows		= drows;
    vap->va_dcolumns		= dcolumns;
    vap->va_equations		= 0;
    vap->va_frequencies_valid	= false;

    /*
     * Allocate the equation bitmap.
     */
    if ((vap->va_bitmap = calloc((dcells + 31) / 32,
		    sizeof(uint32_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	vnacal_apply_free(vap);
	return NULL;
    }

    /*
     * Allocate the data vector.  This object holds the frequency
     * vector and two internal matrices, A and B side-by-side.  A is
     * a coefficient matrix with number of columns equal to the number
     * of S-parameters (drows * dcolumns), and B is a column vector of
     * constant terms.  Each row represents an equation in a linear
     * system of equations: A s = B, where s is a column vector of
     * all the S-parameters taken row by row.  In order to solve for
     * all S-parameters, it's necessary that va_data has no fewer than
     * drows * dcolumns rows (equations).  It may have more than that,
     * however, i.e. the system may be over-determined, in which case we
     * find the solution with least squared error using QR decomposition
     * and back substitution.
     *
     * The number of columns in va_data is drows x dcolumns (row of A)
     * plus 1 (row of B).  The number of rows is not known, but it must
     * be at least drows * dcolumns.  We use that as an initial guess,
     * and resize as necessary if additional equations are given.
     *
     * The va_equations variable gives the number of rows of va_data
     * that are actually in-use.
     */
    if ((vap->va_data = vnadata_alloc_and_init(dfrequencies, drows * dcolumns,
		    drows * dcolumns + 1, VPT_UNDEF)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	vnacal_apply_free(vap);
	return NULL;
    }
    return vap;
}

/*
 * vnacal_apply_set_frequency_vector: set the frequency vector
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @frequency_vector: vector of measured frequencies
 */
int vnacal_apply_set_frequency_vector(vnacal_apply_t *vap,
	const double *frequency_vector)
{
    vnacal_t *vcp;
    int dfrequencies;
    vnacal_etermset_t *etsp;
    double fmin, fmax;

    assert(vap != NULL);
    assert(frequency_vector != NULL);
    dfrequencies = vnadata_get_frequencies(vap->va_data);
    vcp = vap->va_vcp;
    etsp = vcp->vc_set_vector[vap->va_set];
    if (frequency_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply_set_frequency_vector: "
		"invalid NULL frequency vector");
	return -1;
    }
    for (int i = 0; i < dfrequencies - 1; ++i) {
	if (frequency_vector[i] >= frequency_vector[i + 1]) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_apply_set_frequency_vector: "
		    "error: frequencies must be ascending");
	    return -1;
	}
    }
    fmin = _vnacal_etermset_get_fmin_bound(etsp);
    if (frequency_vector[0] < fmin) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply_set_frequency_vector: "
		"frequency out of bounds %.3e < %.3e",
		frequency_vector[0], etsp->ets_frequency_vector[0]);
	return -1;
    }
    fmax = _vnacal_etermset_get_fmax_bound(etsp);
    if (frequency_vector[dfrequencies - 1] > fmax) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_apply_set_frequency_vector: "
		"frequency out of bounds %.3e > %.3e",
		frequency_vector[dfrequencies - 1],
		etsp->ets_frequency_vector[etsp->ets_frequencies - 1]);
	return -1;
    }
    vnadata_set_frequency_vector(vap->va_data, frequency_vector);
    vap->va_frequencies_valid = true;
    return 0;
}

/*
 * check_map: validate bitmap
 */
static int check_map(const char *function, vnacal_apply_t *vap, const int *map)
{
    int vports = MAX(vap->va_vrows, vap->va_vcolumns);
    int dports = MAX(vap->va_drows, vap->va_dcolumns);
    uint32_t bitmap[(dports + 31) / 32];
    vnacal_t *vcp = vap->va_vcp;

    (void)memset((void *)bitmap, 0, sizeof(bitmap));
    for (int vport = 0; vport < vports; ++vport) {
	int dport = map[vport];

	if (dport == -1) {
	    continue;
	}
	if (dport < 0 || dport >= dports) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: map[%d] = %d must be in -1..%d",
		    function, vport, dport, dports - 1);
	    return -1;
	}
	if (bitmap[dport / 32] & (1U << (dport % 32))) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: DUT port %d given multiple times in map",
		    function, dport);
	    return -1;
	}
	bitmap[dport / 32] |= 1U << (dport % 32);
    }
    return 0;
}

/*
 * apply_add_column_common: common code to add a column of measurements
 *   @function: name of user-callable function
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @vcolumn: column of the calibration matrix being added
 *   @stride: spacing between vector elements (can be slice of matrix)
 *   @column_vector: pointer to first element of the vector
 *   @map: map of VNA port to DUT port, MAX(vrows, vcolumns) long or NULL
 */
static int apply_add_column_common(const char *function, vnacal_apply_t *vap,
	int vcolumn, int stride, const double complex *const *column_vector,
	const int *map)
{
    vnacal_t *vcp = vap->va_vcp;
    vnacal_etermset_t *etsp = vcp->vc_set_vector[vap->va_set];
    vnadata_t *vdp = vap->va_data;
    const double *frequency_vector = vnadata_get_frequency_vector(vdp);
    int vrows    = vap->va_vrows;
    int vcolumns = vap->va_vcolumns;
    int dfrequencies = vnadata_get_frequencies(vdp);
    int drows    = vap->va_drows;
    int dcolumns = vap->va_dcolumns;
    int dcolumn;
    int segment = 0;

    /*
     * Check arguments.
     */
    if (vcolumn < 0 || vcolumn >= vcolumns) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: invalid vcolumn %d: must be between 0 and %d",
		function, vcolumn, vcolumns - 1);
	return -1;
    }
    if (map != NULL) {
	if (check_map(function, vap, map) == -1) {
	    return -1;
	}
	dcolumn = map[vcolumn];
    } else {
	dcolumn = vcolumn;
    }

    /*
     * If we're being offered measurements for a DUT column that's not
     * connected or isn't wanted, ignore them.
     */
    if (dcolumn < 0 || dcolumn >= dcolumns) {
	return 0;
    }

    /*
     * For each measurement...
     */
    for (int vrow = 0; vrow < vrows; ++vrow) {
	int drow = (map != NULL) ? map[vrow] : vrow;

	/*
	 * If the measured row isn't connected or is outside the
	 * bounds of the DUT matrix, skip it.
	 */
	if (drow < 0 || drow >= drows) {
	    continue;
	}

	/*
	 * Otherwise, we need to add an equation for this measurement.
	 * Resize va_data as required.
	 */
	if (vap->va_equations >= vnadata_get_rows(vdp)) {
	    int new_rows = MAX(4,
		    vap->va_equations + vap->va_equations / 2);

	    if (vnadata_resize(vdp, dfrequencies, new_rows,
			drows * dcolumns + 1, VPT_UNDEF) == -1) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"malloc: %s", strerror(errno));
		return -1;
	    }
	}

	/*
	 * For each frequency...
	 */
	for (int findex = 0; findex < dfrequencies; ++findex) {
	    double f = frequency_vector[findex];
	    double complex *data_row =
		&vnadata_get_matrix(vdp, findex)[
		vap->va_equations * (drows * dcolumns + 1)];
#define A(drow, dcolumn)	(data_row[dcolumns * (drow) + (dcolumn)])
#define B()			(data_row[drows * dcolumns])
#define V(vrow, findex)		(column_vector[stride * (vrow)][findex])

	    /*
	     * Start by adding a coefficient of 1 for the S-parameter
	     * primarily associated with this equation.
	     */
	    A(drow, dcolumn) = 1.0;

	    /*
	     * Next, for each connected port, add the contribution of
	     * the port mismatch term.  At the same time, fill in the
	     * elements of B vector.
	     */
	    for (int k = 0; k < vrows; ++k) {
		int d = (map != NULL) ? map[k] : k;

		if (k == vrow || (d >= 0 && d < dcolumns)) {
		    vnacal_error_terms_t *etp =
			&etsp->ets_error_term_matrix[k * vcolumns + vcolumn];
		    double complex e0 = interpolate(etsp, etp, &segment, 0, f);
		    double complex e1 = interpolate(etsp, etp, &segment, 1, f);
		    double complex x;

		    x = (V(k, findex) - e0) / e1;
		    if (d >= 0 && d < dcolumns) {
			double complex e2;

			e2 = interpolate(etsp, etp, &segment, 2, f);
			A(drow, d) += e2 * x;
		    }
		    if (k == vrow) {
			B() = x;
		    }
		}
	    }

	    /*
	     * Record that we have an equation for the associated
	     * S parameter.  This is used in vnacal_apply_get_data to
	     * catch under-determined systems, i.e. where the supplied
	     * measurements aren't sufficient to solve for the S
	     * parameters.
	     */
	    vap->va_bitmap[(drow * dcolumns + dcolumn) / 32] |=
		1U << (drow * dcolumns + dcolumn) % 32;
#undef V
#undef B
#undef A
	}
	++vap->va_equations;
    }
    return 0;
}

/*
 * vnacal_apply_add_column: add a column of measurements with map
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @vcolumn: column of the calibration matrix being added
 *   @column_vector: pointer to vector, vrows long, of pointers
 *      to vectors of voltage measurements, dfrequencies long
 *   @map: map of VNA port to DUT port, MAX(vrows, vcolumns) long
 *
 * When the DUT has more ports than the VNA, this function gives the
 * mapping for this set of measurements.  The map must be the greater
 * dimension of the calibration matrix.  DUT port indices start at zero.
 * The special value -1 indicates that the VNA port and any unused DUT
 * ports are connected to terminators.
 */
int vnacal_apply_add_column(vnacal_apply_t *vap, int vcolumn,
	const double complex *const *column_vector, const int *map)
{
    assert(vap != NULL);
    assert(column_vector != NULL);
    return apply_add_column_common("vnacal_apply_add_column", vap,
	    vcolumn, /*stride*/1, column_vector, map);
}

/*
 * vnacal_apply_add_matrix: add matrix of vector with port map
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @matrix: vrows x vcolumns matrix of measured voltage ratios
 *   @map: map of VNA port to DUT port, MAX(vrows, vcolumns) long
 */
int vnacal_apply_add_matrix(vnacal_apply_t *vap,
	const double complex *const *matrix, const int *map)
{
    assert(vap != NULL);
    assert(matrix != NULL);
    for (int vcolumn = 0; vcolumn < vap->va_vcolumns; ++vcolumn) {
	int rc;

	rc = apply_add_column_common("vnacal_apply_add_matrix", vap,
			vcolumn, vap->va_vcolumns, &matrix[vcolumn], map);
	if (rc == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * vnacal_apply_get_data: create the s-parameter matrix
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @s_parameters: caller-allocated vnadata_t structure
 */
int vnacal_apply_get_data(const vnacal_apply_t *vap, vnadata_t *s_parameters)
{
    vnacal_t *vcp;
    int drows, dcolumns, dfrequencies;

    assert(vap != NULL);
    assert(s_parameters != NULL);
    vcp = vap->va_vcp;
    drows = vap->va_drows;
    dcolumns = vap->va_dcolumns;
    dfrequencies = vnadata_get_frequencies(vap->va_data);

    /*
     * Test that each S parameter has an associated equation.
     */
    for (int drow = 0; drow < drows; ++drow) {
	for (int dcolumn = 0; dcolumn < dcolumns; ++dcolumn) {
	    int dcell = drow * dcolumns + dcolumn;

	    if (!(vap->va_bitmap[dcell / 32] & (1U << (dcell % 32)))) {
		_vnacal_error(vcp, VNAERR_USAGE, "vnacal_apply_get_data: "
			"no equation for S%d%s%d", drow + 1,
			(drows >= 10 || dcolumns >= 10) ? "," : "",
			dcolumn + 1);
		return -1;
	    }
	}
    }

    /*
     * Init the output matrix.
     */
    if (vnadata_init(s_parameters, dfrequencies, vap->va_drows,
		vap->va_dcolumns, VPT_S) == -1) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"vnadata_init: %s", strerror(errno));
	return -1;
    }
    (void)vnadata_set_frequency_vector(s_parameters,
	    vnadata_get_frequency_vector(vap->va_data));
    (void)vnadata_set_all_z0(s_parameters,
	    vcp->vc_set_vector[vap->va_set]->ets_z0);

    /*
     * Calculate S parameters.
     */
    for (int findex = 0; findex < dfrequencies; ++findex) {
	int equations = vap->va_equations;
	complex double a[equations][drows * dcolumns];
	complex double b[equations];

	for (int equation = 0; equation < equations; ++equation) {
	    for (int cell = 0; cell < drows * dcolumns; ++cell) {
		a[equation][cell] = vnadata_get_cell(vap->va_data,
			findex, equation, cell);
	    }
	    b[equation] = vnadata_get_cell(vap->va_data,
		    findex, equation, drows * dcolumns);
	}
	_vnacommon_qrsolve(vnadata_get_matrix(s_parameters, findex),
		*a, b, equations, drows * dcolumns, 1);
    }
    return 0;
}

/*
 * vnacal_apply_free: free the input object
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 */
void vnacal_apply_free(vnacal_apply_t *vap)
{
    if (vap != NULL) {
	vnadata_free(vap->va_data);
	free((void *)vap->va_bitmap);
	free((void *)vap);
    }
}
