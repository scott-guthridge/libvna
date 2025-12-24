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

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * vnadata_set_name: set a name for this device
 *   @vdp: a pointer to the vnadata_t structure
 *   @name: name for the device
 *
 *   Sets a name for the device.  The name is limited to 31 characters
 *   plus a terminating NUL.  If the given name is longer than this, it
 *   will be truncated.
 */
int vnadata_set_name(vnadata_t *vdp, const char *name)
{
    vnadata_internal_t *vdip;

    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)strncpy(vdip->vdi_name, name, VNADATA_MAX_NAME);
    vdip->vdi_name[VNADATA_MAX_NAME] = '\000';
    vdip->vdi_flags |= VF_NAME_SET;
    return 0;
}
