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
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_delete_calibration: delete a calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
int vnacal_delete_calibration(vnacal_t *vcp, int ci)
{
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (ci >= 0 && ci < vcp->vc_calibration_allocation) {
	vnacal_calibration_t *calp = vcp->vc_calibration_vector[ci];

	if (calp != NULL) {
	    _vnacal_calibration_free(vcp->vc_calibration_vector[ci]);
	    vcp->vc_calibration_vector[ci] = NULL;
	    return 0;
	}
    }
    errno = ENOENT;
    return -1;
}
