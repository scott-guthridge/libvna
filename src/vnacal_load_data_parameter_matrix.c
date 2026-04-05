/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2026 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include "vnaconv.h"
#include "vnadata.h"


/*
 * vnacal_load_data_parameter_matrix: load a parameter matrix from file
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @filename: filename of Touchstone of NPD file to load
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
int vnacal_load_data_parameter_matrix(vnacal_t *vcp,
	const char *filename, int *parameter_matrix,
	size_t parameter_matrix_size)
{
    vnadata_t *vdp = NULL;
    int rv = -1;

    if ((vdp = vnadata_alloc(vcp->vc_error_fn, vcp->vc_error_arg)) == NULL) {
	goto out;
    }
    if (vnadata_load(vdp, filename) == -1) {
	goto out;
    }
    if ((rv = vnacal_make_data_parameter_matrix(vcp, vdp,
		parameter_matrix, parameter_matrix_size)) == -1) {
	goto out;
    }
    rv = 0;

out:
    vnadata_free(vdp);
    return rv;
}
