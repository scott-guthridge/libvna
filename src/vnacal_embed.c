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
#include "vnadata.h"


/*
 * vnacal_embed: embed DUT into test fixture using network parameter data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vdp_in: network parameter data to embed
 *   @vdp_out: takes frequency vector, z0 and type as input; returns data
 *   @fixture_matrix: parameter matrix of fixture
 *   @fixture_ports: number of fixture ports (must be even)
 *
 *   For in-place conversion, vdp_out can be the same as vdp_in.
 */
int vnacal_embed(vnacal_t *vcp, vnadata_t *vdp_in, vnadata_t *vdp_out,
	const int *fixture_matrix, int fixture_ports)
{
    const int rows = vnadata_get_rows(vdp_in);
    const int columns = vnadata_get_columns(vdp_in);
    const int ports = MAX(rows, columns);
    int target_parameter_matrix[rows][columns];
    int result_parameter_matrix[rows][columns];
    vnadata_parameter_type_t parameter_type = vnadata_get_type(vdp_out);
    int rv = -1;

    /*
     * Init matrices before any "goto out" for delete parameter matrix.
     */
    _vnacal_init_parameter_matrix(&target_parameter_matrix[0][0],
	    rows, columns);
    _vnacal_init_parameter_matrix(&result_parameter_matrix[0][0],
	    rows, columns);

    assert(rows == columns);  /* enforced by vnadata for VPT_S */
    if (fixture_ports != 2 * ports) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: fixture_ports must be exactly two times DUT ports",
		__func__);
	goto out;
    }
    if (vnadata_get_rows(vdp_out) != ports ||
	    vnadata_get_columns(vdp_out) != ports) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: vdp_out must have same dimensions as vdp_in", __func__);
	goto out;
    }
    if (_vnacal_make_data_parameter_matrix(__func__, vcp, vdp_in,
		&target_parameter_matrix[0][0],
		sizeof(target_parameter_matrix)) == -1) {
	goto out;
    }
    if ((rv = _vnacal_embed_parameter_matrix(__func__, &_vnacal_embed_ops,
	    vcp, &target_parameter_matrix[0][0], ports, fixture_matrix,
	    &result_parameter_matrix[0][0],
	    sizeof(result_parameter_matrix))) == -1) {
	goto out;
    }
    if (vnadata_set_type(vdp_out, VPT_S) == -1) {
	goto out;
    }
    if (vnacal_parameter_matrix_to_data(vcp, &result_parameter_matrix[0][0],
		rows, columns, vdp_out) == -1) {
	goto out;
    }
    if (vnadata_convert(vdp_out, vdp_out, parameter_type) == -1) {
	goto out;
    }
    rv = 0;

out:
    (void)vnacal_delete_parameter_matrix(vcp,
	    &result_parameter_matrix[0][0], rows, columns);
    (void)vnacal_delete_parameter_matrix(vcp,
	    &target_parameter_matrix[0][0], rows, columns);
    return rv;
}
