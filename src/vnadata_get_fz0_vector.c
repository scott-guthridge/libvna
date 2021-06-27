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

#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * vnadata_get_fz0_vector: return the z0 vector for the given frequency
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 */
const double complex *vnadata_get_fz0_vector(const vnadata_t *vdp, int findex)
{
    vnadata_internal_t *vdip;

    if (vdp == NULL) {
	errno = EINVAL;
	return NULL;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    if (findex < 0 || findex > vdp->vd_frequencies) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"vnadata_get_fz0_vector: frequency index: %d: out of bounds",
		findex);
	return NULL;
    }
    if (vdip->vdi_flags & VF_PER_F_Z0) {
	return vdip->vdi_z0_vector_vector[findex];
    }
    return vdip->vdi_z0_vector;
}
