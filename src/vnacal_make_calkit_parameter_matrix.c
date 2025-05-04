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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/*
 * vnacal_make_calkit_parameter_matrix: make parameter matrix for kit standard
 *   @function: name of user-called function
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vcdp: data describing the calkit standard
 *   @parameter_matrix: caller-supplied result matrix
 *   @parameter_matrix_size: size in bytes of the result matrix
 *
 * Fill parameter_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_single_reflect* or vnadata_new_add_line* functions.
 * The parameter_matrix_size parameter is the allocation in bytes of the
 * result matrix, used to protect against buffer overrun.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
static int _vnacal_make_calkit_parameter_matrix(const char *function,
	vnacal_t *vcp, const vnacal_calkit_data_t *vcdp,
	int *parameter_matrix, size_t parameter_matrix_size)
{
    vnacal_standard_t *stdp = NULL;
    const char *name;
    int ports = 0;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (vcdp == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: vcdp cannot be NULL", function);
	return -1;
    }
    if ((name = _vnacal_get_calkit_name(vcdp, &ports)) == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: vnacal_calkit_data_t structure is not valid", function);
	return -1;
    }
    if (ports * ports * sizeof(int) > parameter_matrix_size) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: insufficient result matrix allocation", function);
	return -1;
    }
    for (int cell = 0; cell < ports * ports; ++cell) {
	parameter_matrix[cell] = -1;
    }

    /*
     * Allocate and init the standard.
     */
    if ((stdp = malloc(sizeof(vnacal_standard_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	goto error;
    }
    (void)memset((void *)stdp, 0, sizeof(*stdp));
    stdp->std_type = VNACAL_CALKIT;
    if ((stdp->std_name = strdup(name)) == NULL) {
	goto error;
    }
    stdp->std_ports = ports;
    stdp->std_refcount = 0;
    stdp->std_vcp = vcp;
    stdp->std_calkit_data = *vcdp;

    /*
     * Add the parameter structures.
     */
    for (int row = 0; row < ports; ++row) {
	for (int column = 0; column < ports; ++column) {
	    const int cell = ports * row + column;
	    vnacal_parameter_t *vpmrp;

	    vpmrp = _vnacal_alloc_parameter(function, vcp);
	    if (vpmrp == NULL) {
		goto error;
	    }
	    vpmrp->vpmr_type = VNACAL_CALKIT;
	    vpmrp->vpmr_stdp = stdp;
	    vpmrp->vpmr_row = row;
	    vpmrp->vpmr_column = column;
	    parameter_matrix[cell] = vpmrp->vpmr_index;
	}
    }
    stdp->std_refcount = ports * ports;
    return ports;

error:
    for (int row = 0; row < ports; ++row) {
	for (int column = 0; column < ports; ++column) {
	    const int cell = ports * row + column;
	    int idx = parameter_matrix[cell];
	    vnacal_parameter_t *vpmrp;

	    if (idx == -1) {
		break;
	    }
	    vpmrp = vcp->vc_parameter_collection.vprmc_vector[idx];
	    vpmrp->vpmr_type = VNACAL_NEW; /* prevent deletion of standard */
	    _vnacal_release_parameter(vpmrp);
	}
    }
    _vnacal_free_standard(stdp);
    return -1;
}

/*
 * vnacal_make_calkit_parameter: make a parameter for a one-port kit standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vcdp: data describing the calkit standard
 */
int vnacal_make_calkit_parameter(vnacal_t *vcp,
	const vnacal_calkit_data_t *vcdp)
{
    int parameter;

    if (_vnacal_make_calkit_parameter_matrix(__func__, vcp, vcdp,
	    &parameter, sizeof(parameter)) == -1) {
	return -1;
    }
    return parameter;
}

/*
 * vnacal_make_calkit_parameter_matrix: make parameter matrix for kit standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vcdp: data describing the calkit standard
 *   @parameter_matrix: caller-supplied result matrix
 *   @parameter_matrix_size: size in bytes of the result matrix
 *
 * Fill parameter_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_single_reflect* or vnadata_new_add_line* functions.
 * The parameter_matrix_size parameter is the allocation in bytes of the
 * result matrix, used to protect against buffer overrun.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
int vnacal_make_calkit_parameter_matrix(vnacal_t *vcp,
	const vnacal_calkit_data_t *vcdp, int *parameter_matrix,
	size_t parameter_matrix_size)
{
    return _vnacal_make_calkit_parameter_matrix(__func__, vcp, vcdp,
	    parameter_matrix, parameter_matrix_size);
}
