/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * MAX_FORMAT: maximum length of a format specifier
 */
#define MAX_FORMAT	5

/*
 * _vnadata_update_format_string: regenerate vdi_format_string
 *   @vdip:   internal parameter matrix
 */
int _vnadata_update_format_string(vnadata_internal_t *vdip)
{
    char *new_string = NULL;
    char *cur;

    free((void *)vdip->vdi_format_string);
    vdip->vdi_format_string = NULL;
    if (vdip->vdi_format_count != 0) {
	if ((new_string = malloc(vdip->vdi_format_count *
			(MAX_FORMAT + 1))) == NULL) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "malloc: %s", strerror(errno));
	    return -1;
	}
	cur = new_string;
	for (int i = 0;;) {
	    const char *name;

	    name = _vnadata_format_to_name(&vdip->vdi_format_vector[i]);
	    (void)strcpy(cur, name);
	    cur += strlen(cur);
	    if (++i >= vdip->vdi_format_count) {
		break;
	    }
	    *cur++ = ',';
	}
	*cur = '\000';
	vdip->vdi_format_string = new_string;
    }
    return 0;
}
