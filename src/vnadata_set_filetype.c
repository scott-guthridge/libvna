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
 * vnadata_set_filetype: set the file type
 *   @vdp: a pointer to the vnadata_t structure
 *   @type: file type
 *
 *   The default type is VNADATA_FILETYPE_AUTO where the library tries
 *   to intuit the type from the filename.
 */
int vnadata_set_filetype(vnadata_t *vdp, vnadata_filetype_t filetype)
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
    switch (filetype) {
    case VNADATA_FILETYPE_AUTO:
    case VNADATA_FILETYPE_TOUCHSTONE1:
    case VNADATA_FILETYPE_TOUCHSTONE2:
    case VNADATA_FILETYPE_NPD:
	break;

    default:
	_vnadata_error(vdip, VNAERR_USAGE,
		"vnadata_set_filetype: invalid filetype");
	return -1;
    }
    vdip->vdi_filetype = filetype;

    return 0;
}
