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
 * _vnadata_convert_to_z0: convert from simple z0 to frequency-dependent z0
 *   @vdip: internal object pointer
 *
 * Note: If a conversion is done, the simple vector is copied to all
 *       frequency rows, preserving the values.
 */
int _vnadata_convert_to_fz0(vnadata_internal_t *vdip)
{
    if (!(vdip->vdi_flags & VF_PER_F_Z0)) {
	double complex **clfpp = NULL;

	if (vdip->vdi_f_allocation > 0) {
	    clfpp = calloc(vdip->vdi_f_allocation, sizeof(double complex *));
	    if (clfpp == NULL) {
		return -1;
	    }
	    if (vdip->vdi_p_allocation > 0) {
		for (int findex = 0; findex < vdip->vdi_f_allocation;
			++findex) {
		    if ((clfpp[findex] = calloc(vdip->vdi_p_allocation,
				    sizeof(double complex))) == NULL) {
			while (--findex >= 0) {
			    free((void *)clfpp[findex]);
			}
			free((void *)clfpp);
			return -1;
		    }
		    for (int port = 0; port < vdip->vdi_p_allocation; ++port) {
			clfpp[findex][port] = vdip->vdi_z0_vector[port];
		    }
		}
	    }
	}
	free((void *)vdip->vdi_z0_vector);
	vdip->vdi_z0_vector_vector = clfpp;
	vdip->vdi_flags |= VF_PER_F_Z0;
    }
    return 0;
}
