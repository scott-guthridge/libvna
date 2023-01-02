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
 * vnacal_type_to_name: convert error term type to name
 *   @type: type of error terms
 */
const char *vnacal_type_to_name(vnacal_type_t type)
{
    switch (type) {
    case VNACAL_T8:
	return "T8";

    case VNACAL_U8:
	return "U8";

    case VNACAL_TE10:
	return "TE10";

    case VNACAL_UE10:
	return "UE10";

    case VNACAL_T16:
	return "T16";

    case VNACAL_U16:
	return "U16";

    case VNACAL_UE14:
	return "UE14";

    case _VNACAL_E12_UE14:
	return "E12_UE14";

    case VNACAL_E12:
	return "E12";

    default:
	break;
    }
    return NULL;
}
