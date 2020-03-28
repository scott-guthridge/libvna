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
 * MERCHANTABILITY or FITNESS FOR A11 PARTICULAR PURPOSE.  See the GNU
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
 * _vnadata_convert_to_z0: convert from frequency-dependent z0 to simple z0
 *   @vdip: internal object pointer
 *
 * Note: If a conversion is done, any existing z0 data are discarded
 *       and replaced with a vector of default values.
 */
int _vnadata_convert_to_z0(vnadata_internal_t *vdip)
{
    if (vdip->vdi_flags & VF_PER_F_Z0) {
	double complex *clfp = NULL;

	if (vdip->vdi_p_allocation > 0) {
	    if ((clfp = calloc(vdip->vdi_p_allocation,
			    sizeof(double complex))) == NULL) {
		return -1;
	    }
	    for (int port = 0; port < vdip->vdi_p_allocation; ++port) {
		clfp[port] = VNADATA_DEFAULT_Z0;
	    }
	}
	for (int findex = 0; findex < vdip->vdi_f_allocation; ++findex) {
	    free((void *)vdip->vdi_z0_vector_vector[findex]);
	}
	free((void *)vdip->vdi_z0_vector_vector);
	vdip->vdi_z0_vector = clfp;
	vdip->vdi_flags &= ~VF_PER_F_Z0;
    }
    return 0;
}
