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

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * vnadata_get_z0: return the z0 value for the given port
 *   @vdp:  a pointer to the vnadata_t structure
 *   @port: port number (zero-based)
 */
double complex vnadata_get_z0(const vnadata_t *vdp, int port)
{
    vnadata_internal_t *vdip;
    int ports;

    if (vdp == NULL) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    ports = MAX(vdp->vd_rows, vdp->vd_columns);
    if (port < 0 || port > ports) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"vnadata_get_z0: port index: %d: out of bounds", port);
	return HUGE_VAL;
    }
    if (vdip->vdi_flags & VF_PER_F_Z0) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"vnadata_get_z0: per-frequency z0 values are in-use: "
		"vnadata_get_fz0 instead");
	return HUGE_VAL;
    }
    return vdip->vdi_z0_vector[port];
}
