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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vnadata_internal.h"


/*
 * _vnadata_parse_filename: try to determine filetype from filename
 *   @filename: input or output filename
 */
vnadata_filetype_t _vnadata_parse_filename(const char *filename, int *ports)
{
    const char *suffix;

    /*
     * Set *ports to the default of unknown (-1).
     */
    if (ports != NULL) {
	*ports = -1;
    }

    /*
     * If there's no file extension, return VNADATA_FILETYPE_AUTO to
     * indicate unknown.
     */
    if ((suffix = strrchr(filename, '.')) == NULL) {
	return VNADATA_FILETYPE_AUTO;
    }
    ++suffix;

    /*
     * If the file has a suffix of .ts, assume touchstone 2.
     */
    if (strcasecmp(suffix, "ts") == 0) {
	return VNADATA_FILETYPE_TOUCHSTONE2;
    }

    /*
     * If the file has a suffix matching s[0-9]+p, assume touchstone 1.
     */
    if (suffix[0] == 's') {
	const char *cp = suffix + 1;

	if (isdigit(*cp)) {
	    do {
		++cp;
	    } while (isdigit(*cp));
	    if (cp[0] == 'p' && cp[1] == '\000') {
		if (ports != NULL) {
		    *ports = atoi(&suffix[1]);
		}
		return VNADATA_FILETYPE_TOUCHSTONE1;
	    }
	}
    }

    /*
     * If the file ends in .npd, return NPD type.
     */
    if (strcasecmp(suffix, "npd") == 0) {
	return VNADATA_FILETYPE_NPD;
    }

    /*
     * None of the above.  Return VNADATA_FILETYPE_AUTO to indicate
     * unknown.
     */
    return VNADATA_FILETYPE_AUTO;
}
