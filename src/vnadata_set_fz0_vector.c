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
 * vnadata_set_fz0_vector: set the z0 vector for the given frequency
 *   @vdp:       a pointer to the vnadata_t structure
 *   @findex:    frequency index
 *   @z0_vector: new values (length is max of rows and columns)
 */
int vnadata_set_fz0_vector(vnadata_t *vdp, int findex,
	const double complex *z0_vector)
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
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"vnadata_set_fz0_vector: invalid frequency index: %d", findex);
	return -1;
    }
    ports = MAX(vdp->vd_rows, vdp->vd_columns);
    if (!(vdip->vdi_flags & VF_PER_F_Z0)) {
	if (_vnadata_convert_to_fz0(vdip) == -1) {
	    return -1;
	}
    }
    (void)memcpy((void *)vdip->vdi_z0_vector_vector[findex],
	    (void *)z0_vector, ports * sizeof(double complex));
    return 0;
}
