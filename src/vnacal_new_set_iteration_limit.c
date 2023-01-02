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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
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
#include "vnacal_new_internal.h"



/*
 * vnacal_new_set_iteration_limit: set iteration limit for iterative solutions
 *   @vnp: pointer to vnacal_new_t structure
 *   @iterations: maximum number of iterations before failing
 */
int vnacal_new_set_iteration_limit(vnacal_new_t *vnp, int iterations)
{
    vnacal_t *vcp;

    /*
     * Validate arguments.
     */
    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vcp = vnp->vn_vcp;
    if (iterations < 1) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_iteration_limit: "
		"iterations must be at least 1");
	return -1;
    }
    vnp->vn_iteration_limit = iterations;
    return 0;
}
