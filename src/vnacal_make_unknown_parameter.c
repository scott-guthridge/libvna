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
 * vnacal_make_unknown_parameter: create an unknown parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @initial_guess: index of scalar or vector parameter serving as guess
 */
int vnacal_make_unknown_parameter(vnacal_t *vcp, int initial_guess)
{
    vnacal_parameter_t *vstdp_other;
    vnacal_parameter_t *vpmrp;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vstdp_other = _vnacal_get_parameter(vcp, initial_guess);
    if (vstdp_other == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_unknown_parameter: "
		"initial_guess must refer to a valid scalar or vector "
		"parameter");
	return -1;
    }
    vpmrp = _vnacal_alloc_parameter("vnacal_make_unknown_parameter", vcp);
    if (vpmrp == NULL) {
	return -1;
    }
    _vnacal_hold_parameter(vstdp_other);
    vpmrp->vpmr_type = VNACAL_UNKNOWN;
    vpmrp->vpmr_other = vstdp_other;
    return vpmrp->vpmr_index;
}
