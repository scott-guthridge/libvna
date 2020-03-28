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

#include <assert.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * vnadata_add_frequency: add a new frequency entry
 *   @vdp: object pointer
 *   @frequency: new frequency value
 */
int vnadata_add_frequency(vnadata_t *vdp, double frequency)
{
    vnadata_internal_t *vdip;

    /*
     * Check parameters
     */
    if (vdp == NULL || frequency < 0.0) {
	errno = EINVAL;
	return -1;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }

    /*
     * Extend the frequency allocation as needed.
     */
    if (vdp->vd_frequencies + 1 > vdip->vdi_f_allocation) {
	int old_allocation = vdip->vdi_f_allocation;
	int new_allocation = MAX(50, old_allocation + old_allocation / 2);

	if (_vnadata_extend_f(vdip, new_allocation) == -1) {
	    return -1;
	}
    }

    /*
     * Add the new frequency.
     */
    vdp->vd_frequency_vector[vdp->vd_frequencies++] = frequency;

    return 0;
}
