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

#define VNADATA_NO_BOUNDS_CHECK

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"
#include "vnadata.h"

/*
 * _vnacal_make_data_parameter_matrix: make a parameter matrix from data
 *   @function: called function name
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vdp: network parameter data for a calibration standard
 *   @parameter_matrix: caller-allocated matrix to receive result
 *   @parameter_matrix_size: allocation in bytes of the result matrix
 *
 * Fills parameter_matrix with parameter indices representing the network
 * parameter data.  Depending on the dimensions, the resulting buffer
 * is suitable for passing to the vnacal_new_add_* functions with
 * vnacal_new_add_mapped_matrix_m being the most general.  Rows and
 * columns must match the dimensions of the data and currently must be
 * square.  These parameters are unnecessary, but included for defensive
 * programming against overrunning the caller's buffer.  Caller can delete
 * the returned parameters by a call to vnacal_delete_parameter_matrix.
 */
static int _vnacal_make_data_parameter_matrix(const char *function,
	vnacal_t *vcp, const vnadata_t *vdp,
	int *parameter_matrix, size_t parameter_matrix_size)
{
    vnadata_t *vdp_copy = NULL;
    vnacal_standard_t *stdp = NULL;
    vnacal_data_standard_t *vdsp = NULL;
    int rows, columns, ports;
    int frequencies;
    bool has_fz0;
    double *frequency_vector;
    double complex **data_matrix;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (vdp == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: vdp cannot be NULL", function);
	return -1;
    }
    rows = vnadata_get_rows(vdp);
    columns = vnadata_get_columns(vdp);
    if (rows * columns * sizeof(int) > parameter_matrix_size) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: insufficient result matrix allocation", function);
	return -1;
    }

    /*
     * Init the parameter matrix for the error clean-up code at the end.
     */
    for (int cell = 0; cell < rows * columns; ++cell) {
	parameter_matrix[cell] = -1;
    }

    /*
     * If the network parameter data is not S parameters, convert.
     */
    if (vnadata_get_type(vdp) != VPT_S) {
	vdp_copy = vnadata_alloc(vcp->vc_error_fn, vcp->vc_error_arg);
	if (vdp_copy == NULL)
	    goto error;
	if (vnadata_convert(vdp, vdp_copy, VPT_S) == -1)
	    goto error;
	vdp = vdp_copy;
    }

    /*
     * Parameter matrices can in general be rectangular to support
     * partially known standards, e.g. maybe we know only a single
     * column of the S parameter matrix of some standard when making the
     * calibration in T8 or TE10 parameters.  The data for such a standard
     * could be saved in a vnadata_t structure as VPT_UNDEF, and we could
     * accept that here.  The main problem, however, is that rectangular
     * S-parameter matrices are not reference impedance renormalizable, so
     * the z0 values would have to match at evaluation time.  We may add
     * that functionality later, but for now, support only square data.
     */
    assert(rows == columns);
    ports = rows;
    frequencies = vnadata_get_frequencies(vdp);

    /*
     * Allocate and init the vnacal_standard_t structure.
     */
    if ((stdp = malloc(sizeof(vnacal_standard_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	goto error;
    }
    (void)memset((void *)stdp, 0, sizeof(*stdp));
    stdp->std_type = VNACAL_DATA;
    if ((stdp->std_name = strdup(vnadata_get_name(vdp))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "strdup: %s", strerror(errno));
	goto error;
    }
    stdp->std_ports = ports;
    stdp->std_refcount = 0;
    stdp->std_vcp = vcp;
    vdsp = &stdp->std_data_standard;
    vdsp->vds_frequencies = frequencies;
    if ((frequency_vector = calloc(frequencies, sizeof(double))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto error;
    }
    (void)memcpy((void *)frequency_vector,
	    (void *)vnadata_get_frequency_vector(vdp),
	    frequencies * sizeof(double));
    vdsp->vds_frequency_vector = frequency_vector;
    vdsp->vds_has_fz0 = has_fz0 = vnadata_has_fz0(vdp);
    vdsp->vds_segment = 0;
    if (has_fz0) {
	double complex **z0_vector_vector;

	z0_vector_vector = calloc(ports, sizeof(double complex *));
	if (z0_vector_vector == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    goto error;
	}
	vdsp->u.vds_z0_vector_vector = z0_vector_vector;
	for (int port = 0; port < ports; ++port) {
	    double complex *vector;

	    vector = calloc(frequencies, sizeof(double complex));
	    if (vector == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s",
			strerror(errno));
		goto error;
	    }
	    z0_vector_vector[port] = vector;
	    for (int findex = 0; findex < frequencies; ++findex) {
		vector[findex] = vnadata_get_fz0(vdp, findex, port);
	    }
	}
    } else {
	double complex *z0_vector;

	if ((z0_vector = calloc(ports, sizeof(double complex))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	    goto error;
	}
	vdsp->u.vds_z0_vector = z0_vector;
	(void)memcpy((void *)z0_vector, (void *)vnadata_get_z0_vector(vdp),
		ports * sizeof(double complex));
    }
    data_matrix = calloc(rows * columns, sizeof(double complex *));
    if (data_matrix == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto error;
    }
    vdsp->vds_data = data_matrix;
    for (int row = 0; row < rows; ++row) {
	for (int column = 0; column < columns; ++column) {
	    const int cell = row * columns + column;
	    double complex *vector;

	    vector = calloc(frequencies, sizeof(double complex));
	    if (vector == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s",
			strerror(errno));
		goto error;
	    }
	    data_matrix[cell] = vector;
	    (void)vnadata_get_to_vector(vdp, row, column, vector);
	}
    }
    data_matrix = NULL;  /* now owned by standard */

    /*
     * Fill the resulting parameter matrix.
     */
    for (int row = 0; row < rows; ++row) {
	for (int column = 0; column < columns; ++column) {
	    vnacal_parameter_t *vpmrp;

	    vpmrp = _vnacal_alloc_parameter(__func__, vcp);
	    if (vpmrp == NULL) {
		goto error;
	    }
	    vpmrp->vpmr_type = VNACAL_DATA;
	    vpmrp->vpmr_stdp = stdp;
	    vpmrp->vpmr_row = row;
	    vpmrp->vpmr_column = column;
	    parameter_matrix[columns * row + column] = vpmrp->vpmr_index;
	}
    }
    stdp->std_refcount = rows * columns;
    vnadata_free(vdp_copy);
    return ports;

error:
    for (int cell = 0; cell < rows * columns; ++cell) {
	int idx = parameter_matrix[cell];
	vnacal_parameter_t *vpmrp;

	if (idx == -1) {
	    break;
	}
	vpmrp = vcp->vc_parameter_collection.vprmc_vector[idx];
	vpmrp->vpmr_type = VNACAL_NEW; /* prevent deletion of standard */
	_vnacal_release_parameter(vpmrp);
    }
    _vnacal_free_standard(stdp);
    vnadata_free(vdp_copy);
    return -1;
}

/*
 * vnacal_make_data_parameter: make a parameter from 1x1 network parameter data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vdp: network parameter data for a calibration standard
 *
 * Dimensions of vdp must be 1x1. Data must be convertable to
 * S-parameters.  Automatically handles parameter
 * conversion, interpolation and renormalization as needed.
 */
int vnacal_make_data_parameter(vnacal_t *vcp, const vnadata_t *vdp)
{
    int parameter;

    if (_vnacal_make_data_parameter_matrix(__func__,
		vcp, vdp, &parameter, sizeof(parameter)) == -1) {
	return -1;
    }
    return parameter;
}

/*
 * vnacal_make_data_parameter_matrix: make a parameter matrix from data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vdp: network parameter data for a calibration standard
 *   @parameter_matrix: caller-allocated matrix to receive result
 *   @parameter_matrix_size: size in bytes of the result matrix
 *
 * Fill parameter_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_* functions.  Automatically handles parameter
 * conversion, interpolation and renormalization.  Data must be
 * convertable to S-parameters.  The parameter_matrix_size parameter is
 * the allocation in bytes of the result matrix, used to protect against
 * buffer overrun.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
int vnacal_make_data_parameter_matrix(vnacal_t *vcp,
	const vnadata_t *vdp, int *parameter_matrix,
	size_t parameter_matrix_size)
{
    return _vnacal_make_data_parameter_matrix(__func__,
	    vcp, vdp, parameter_matrix, parameter_matrix_size);
}
