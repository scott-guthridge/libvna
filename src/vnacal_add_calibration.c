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
 * vnacal_add_calibration: add a newly solved calibration to the vnacal_t
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @name: name of the new calibration
 *   @vnp: pointer to vnacal_new_t structure
 *
 * If name exists, this function replaces the previous calibration.
 */
int vnacal_add_calibration(vnacal_t *vcp, const char *name, vnacal_new_t *vnp)
{
    int ci;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_add_calibration: invalid vnp");
	return -1;
    }
    if (vnp->vn_vcp != vcp) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_add_calibration: "
		"new calibration not associated with this vnacal_t structure");
	return -1;
    }
    if (vnp->vn_calibration == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_add_calibration: "
		"need to call vnacal_new_solve first");
	return -1;
    }
    if ((ci = _vnacal_add_calibration_common(__func__, vcp,
		    vnp->vn_calibration, name)) == -1) {
	return -1;
    }
    vnp->vn_calibration = NULL;

    return ci;
}
