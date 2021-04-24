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

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_find_calibration: find a calibration by name
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @name: name of calibration to find
 */
int vnacal_find_calibration(const vnacal_t *vcp, const char *name)
{
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    for (int ci = 0; ci < vcp->vc_calibration_allocation; ++ci) {
	vnacal_calibration_t *calp = vcp->vc_calibration_vector[ci];

	if (calp != NULL && strcmp(calp->cal_name, name) == 0) {
	    return ci;
	}
    }
    errno = ENOENT;
    return -1;
}
