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
 * vnadata_get_type_name: convert parameter type to name
 *   @type: parameter type
 */
const char *vnadata_get_type_name(vnadata_parameter_type_t type)
{
    switch (type) {
    case VPT_UNDEF:
	return "undefined";
    case VPT_S:
	return "S";
    case VPT_T:
	return "T";
    case _VPT_U:
	return "U";
    case VPT_Z:
	return "Z";
    case VPT_Y:
	return "Y";
    case VPT_H:
	return "H";
    case VPT_G:
	return "G";
    case VPT_A:
	return "A";
    case VPT_B:
	return "B";
    case VPT_ZIN:
	return "Zin";
    default:
	break;
    }
    return NULL;
}
