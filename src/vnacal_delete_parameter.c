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
#include "vnacal_internal.h"


/*
 * vnacal_delete_parameter: delete the parameter with given index
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: parameter to delete
 */
int vnacal_delete_parameter(vnacal_t *vcp, int parameter)
{
    vnacal_parameter_t *vpmrp;

    if (parameter < VNACAL_PREDEFINED_PARAMETERS) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_delete_parameter: %s",
		strerror(errno));
	return -1;
    }
    vpmrp = _vnacal_get_parameter(vcp, parameter);
    if (vpmrp == NULL || vpmrp->vpmr_deleted) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_delete_parameter: "
		"%d: nonexistent parameter", parameter);
	return -1;
    }
    vpmrp->vpmr_deleted = true;
    _vnacal_release_parameter(vpmrp);
    return 0;
}
