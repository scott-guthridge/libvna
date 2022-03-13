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
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_free: free a vnacal_t structure
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
void vnacal_free(vnacal_t *vcp)
{
    if (vcp != NULL && vcp->vc_magic == VC_MAGIC) {
	while (vcp->vc_new_head.l_forw != &vcp->vc_new_head) {
	    vnacal_new_t *vnp = (vnacal_new_t *)((char *)(vcp->vc_new_head.
			l_forw) - offsetof(vnacal_new_t, vn_next));

	    vnacal_new_free(vnp);
	}
	(void)vnaproperty_delete(&vcp->vc_properties, ".");
	assert(vcp->vc_properties == NULL);
	_vnacal_teardown_parameter_collection(vcp);
	vcp->vc_magic = -1;
	free((void *)vcp->vc_filename);
	free((void *)vcp);
    }
}
