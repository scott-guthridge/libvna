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
 * _vnadata_set_simple_format: set the format string (1 parameter)
 *   @vdip:   internal parameter matrix
 *   @type: parameter type
 *   @format: coordinate system
 */
int _vnadata_set_simple_format(vnadata_internal_t *vdip,
	vnadata_parameter_type_t type, vnadata_format_t format)
{
    vnadata_format_descriptor_t *vfdp_new = NULL;
    int rc;

    /*
     * Allocate the new format vector, length 1.
     */
    if ((vfdp_new = malloc(sizeof(vnadata_format_descriptor_t))) == NULL) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	goto out;
    }
    (void)memset((void *)vfdp_new, 0, sizeof(*vfdp_new));
    vfdp_new->vfd_parameter = type;
    vfdp_new->vfd_format = format;

    /*
     * Install the new vector.
     */
    free((void *)vdip->vdi_format_vector);
    vdip->vdi_format_vector = vfdp_new;
    vdip->vdi_format_count = 1;

    /*
     * Update the format string.
     */
    if (_vnadata_update_format_string(vdip) == -1) {
	goto out;
    }
    rc = 0;

out:
    return rc;
}
