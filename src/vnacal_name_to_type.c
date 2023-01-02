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

#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_name_to_type: convert error-term type name to enum
 *   @name: name of error term type (case insensitive)
 */
vnacal_type_t vnacal_name_to_type(const char *name)
{
    if (name == NULL) {
	return VNACAL_NOTYPE;
    }
    switch (name[0]) {
    case 'E':
    case 'e':
	if (strcasecmp(name, "E12") == 0) {
	    return VNACAL_E12;
	}
	break;

    case 'T':
    case 't':
	switch (name[1]) {
	case '1':
	    if (strcasecmp(name, "T16") == 0) {
		return VNACAL_T16;
	    }
	    break;

	case '8':
	    if (strcasecmp(name, "T8") == 0) {
		return VNACAL_T8;
	    }
	    break;

	case 'E':
	case 'e':
	    if (strcasecmp(name, "TE10") == 0) {
		return VNACAL_TE10;
	    }
	    break;

	default:
	    break;
	}
	break;

    case 'U':
    case 'u':
	switch (name[1]) {
	case '1':
	    if (strcasecmp(name, "U16") == 0) {
		return VNACAL_U16;
	    }
	    break;

	case '8':
	    if (strcasecmp(name, "U8") == 0) {
		return VNACAL_U8;
	    }
	    break;

	case 'E':
	case 'e':
	    if (strcasecmp(name, "UE10") == 0) {
		return VNACAL_UE10;
	    }
	    if (strcasecmp(name, "UE14") == 0) {
		return VNACAL_UE14;
	    }
	    break;

	default:
	    break;
	}
	break;

    }
    return VNACAL_NOTYPE;
}
