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
 * _vnacal_alloc: allocate a calibration structure
 *   @function: name of user-called function
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
vnacal_t *_vnacal_alloc(const char *function,
	vnaerr_error_fn_t *error_fn, void *error_arg)
{
    vnacal_t *vcp;

    /*
     * Allocate the vnacal_t.
     */
    if ((vcp = (vnacal_t *)malloc(sizeof(vnacal_t))) == NULL) {
	if (error_fn != NULL) {
	    int saved_errno = errno;
	    char message[80];

	    (void)snprintf(message, sizeof(message), "malloc: %s",
		    strerror(errno));
	    message[sizeof(message)-1] = '\000';
	    errno = saved_errno;
	    (*error_fn)(message, error_arg, VNAERR_SYSTEM);
	    errno = saved_errno;
	}
	return NULL;
    }
    (void)memset((void *)vcp, 0, sizeof(vnacal_t));
    vcp->vc_magic = VC_MAGIC;
    vcp->vc_error_fn = error_fn;
    vcp->vc_error_arg = error_arg;
    if (_vnacal_setup_parameter_collection(function, vcp) == -1) {
	vnacal_free(vcp);
	return NULL;
    }
    vcp->vc_fprecision = VNACAL_DEFAULT_DATA_PRECISION;
    vcp->vc_dprecision = VNACAL_DEFAULT_FREQUENCY_PRECISION;
    vcp->vc_new_head.l_forw = &vcp->vc_new_head;
    vcp->vc_new_head.l_back = &vcp->vc_new_head;

    return vcp;
}


/*
 * vnacal_create: create a new calibration structure
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
vnacal_t *vnacal_create(vnaerr_error_fn_t *error_fn, void *error_arg)
{
    return _vnacal_alloc("vnacal_create", error_fn, error_arg);
}
