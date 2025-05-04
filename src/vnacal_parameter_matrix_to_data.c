/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2025 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * vnacal_parameter_matrix_to_data: convert parameter matrix to parameter data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter_matrix: index of parameter
 *   @rows: number of rows in parameter matrix
 *   @columns: number of columns in parameter matrix
 *   @vdp: takes frequency vector, z0 and type as input; returns data
 */
int vnacal_parameter_matrix_to_data(vnacal_t *vcp,
	const int *parameter_matrix, int rows, int columns, vnadata_t *vdp)
{
    vnadata_parameter_type_t requested_type;
    int frequencies;
    vnacal_parameter_t **vpmrp_matrix = NULL;
    vnacal_parameter_matrix_map_t *vpmmp = NULL;
    int rc = -1;

    /*
     * Validate parameters.
     */
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	goto out;
    }
    if (parameter_matrix == NULL || rows < 0 || columns < 0) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid parameter matrix",
		__func__);
	goto out;
    }
    if (vdp == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: vdp cannot be null",
		__func__);
	goto out;
    }
    requested_type = vnadata_get_type(vdp);
    frequencies = vnadata_get_frequencies(vdp);
    if (vnadata_get_rows(vdp) != rows || vnadata_get_columns(vdp) != columns) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: data matrix dimensions must be %d x %d to match "
		"parameter matrix",
		__func__, rows, columns);
	goto out;
    }
    if (rows != columns && requested_type != VPT_UNDEF) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: parameter type must be VPT_UNDEF when parameter "
		"matrix is rectangular",
		__func__);
	goto out;
    }

    /*
     * If the requested type is not VPT_UNDEF or VPT_S, set it to VPT_S.
     */
    if (requested_type != VPT_UNDEF && requested_type != VPT_S) {
	if (vnadata_set_type(vdp, VPT_S) == -1) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "vnadata_set_type: %s",
		    strerror(errno));
	    goto out;
	}
    }

    /*
     * Form the parameter pointer matrix from the integer parameter matrix.
     */
    vpmrp_matrix = calloc(rows * columns, sizeof(vnacal_parameter_t *));
    if (vpmrp_matrix == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto out;
    }
    for (int cell = 0; cell < rows * columns; ++cell) {
	vnacal_parameter_t *vpmrp;

	vpmrp = _vnacal_get_parameter(vcp, parameter_matrix[cell]);
	if (vpmrp == NULL) {
	    goto out;
	}
	vpmrp_matrix[cell] = vpmrp;
    }

    /*
     * Analyze the parameter matrix.
     */
    if ((vpmmp = _vnacal_analyze_parameter_matrix(__func__, vcp,
		    vpmrp_matrix, rows, columns, /*initial=*/false)) == NULL) {
	goto out;
    }

    /*
     * Evaluate the parameter matrix at each frequency and store
     * the result into vdp.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	const double complex *z0_vector;
	double complex *data_matrix;
	double f;

	z0_vector = vnadata_get_fz0_vector(vdp, findex);
	data_matrix = vnadata_get_matrix(vdp, findex);
	f = vnadata_get_frequency(vdp, findex);
	rc = _vnacal_eval_parameter_matrix_i(__func__,
		vpmmp, f, z0_vector, data_matrix);
	if (rc == -1) {
	    goto out;
	}
    }

    /*
     * If needed, convert the type.
     */
    if (requested_type != VPT_UNDEF && requested_type != VPT_S) {
	if (vnadata_convert(vdp, vdp, requested_type) == -1) {
	    abort();
	}
    }
    rc = 0;

out:
    _vnacal_free_parameter_matrix_map(vpmmp);
    free((void *)vpmrp_matrix);
    return rc;
}
