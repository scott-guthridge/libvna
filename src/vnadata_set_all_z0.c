/*
 * Electrical Network Parameter Conversion Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * vnadata_set_all_z0: set all z0's to the same value
 *   @vdp:  vnadata object pointer
 *   @z0:   new value
 */
int vnadata_set_all_z0(vnadata_t *vdp, double complex z0)
{
    vnadata_internal_t *vdip;
    int ports;

    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    ports = MAX(vdp->vd_rows, vdp->vd_columns);
    if (vdip->vdi_flags & VF_PER_F_Z0) {
	if (_vnadata_convert_to_z0(vdip) == -1) {
	    return -1;
	}
    }
    for (int port = 0; port < ports; ++port) {
	vdip->vdi_z0_vector[port] = z0;
    }
    return 0;
}
