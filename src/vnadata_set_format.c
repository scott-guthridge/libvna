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
 * parse_format: helper function for vnadata_set_format
 *   @vfdp: format descriptor
 *   @format: format string
 */
static int parse_format(vnadata_format_descriptor_t *vfdp, const char *format)
{
    const char *cur = format;

    /*
     * Set defaults.
     */
    vfdp->vfd_parameter = VPT_UNDEF;
    vfdp->vfd_format = VNADATA_FORMAT_REAL_IMAG;

    /*
     * Switch on the first letter of the parameter specifier.
     */
    switch (cur[0]) {
    case '\000':
	goto bad_format;

    case 'a':
	++cur;
	vfdp->vfd_parameter = VPT_A;
	goto parse_coordinates;

    case 'b':
	++cur;
	vfdp->vfd_parameter = VPT_B;
	goto parse_coordinates;

    case 'd':
	goto parse_coordinates;		/* db */

    case 'g':
	++cur;
	vfdp->vfd_parameter = VPT_G;
	goto parse_coordinates;

    case 'h':
	++cur;
	vfdp->vfd_parameter = VPT_H;
	goto parse_coordinates;

    case 'i':
	if (cur[1] == 'l') {
	    cur += 2;
	    vfdp->vfd_parameter = VPT_S;
	    vfdp->vfd_format = VNADATA_FORMAT_IL;
	    break;
	}
	break;

    case 'm':
	goto parse_coordinates;		/* ma */

    case 'p':
	if (strncmp(cur, "prc", 3) == 0) {
	    cur += 3;
	    vfdp->vfd_parameter = VPT_ZIN;
	    vfdp->vfd_format = VNADATA_FORMAT_PRC;
	    break;
	}
	if (strncmp(cur, "prl", 3) == 0) {
	    cur += 3;
	    vfdp->vfd_parameter = VPT_ZIN;
	    vfdp->vfd_format = VNADATA_FORMAT_PRL;
	    break;
	}
	break;

    case 'r':
	if (cur[1] == 'l') {
	    cur += 2;
	    vfdp->vfd_parameter = VPT_S;
	    vfdp->vfd_format = VNADATA_FORMAT_RL;
	    break;
	}
	goto parse_coordinates;		/* ri */

    case 's':
	if (strncmp(cur, "src", 3) == 0) {
	    cur += 3;
	    vfdp->vfd_parameter = VPT_ZIN;
	    vfdp->vfd_format = VNADATA_FORMAT_SRC;
	    break;
	}
	if (strncmp(cur, "srl", 3) == 0) {
	    cur += 3;
	    vfdp->vfd_parameter = VPT_ZIN;
	    vfdp->vfd_format = VNADATA_FORMAT_SRL;
	    break;
	}
	++cur;
	vfdp->vfd_parameter = VPT_S;
	goto parse_coordinates;

    case 't':
	++cur;
	vfdp->vfd_parameter = VPT_T;
	goto parse_coordinates;

    case 'u':
	++cur;
	vfdp->vfd_parameter = VPT_U;
	goto parse_coordinates;

    case 'v':
	if (strncmp(cur, "vswr", 4) == 0) {
	    cur += 4;
	    vfdp->vfd_parameter = VPT_S;
	    vfdp->vfd_format = VNADATA_FORMAT_VSWR;
	    break;
	}
	break;

    case 'y':
	++cur;
	vfdp->vfd_parameter = VPT_Y;
	goto parse_coordinates;

    case 'z':
	if (strncmp(cur, "zin", 3) == 0) {
	    cur += 3;
	    vfdp->vfd_parameter = VPT_ZIN;
	    goto parse_coordinates;
	}
	++cur;
	vfdp->vfd_parameter = VPT_Z;
	goto parse_coordinates;

    parse_coordinates:
	switch (*cur) {
	case 'd':
	    if (strncmp(cur, "db", 2) == 0) {
		if (vfdp->vfd_parameter == VPT_ZIN) {	/* no ZindB */
		    goto bad_format;
		}
		cur += 2;
		vfdp->vfd_format = VNADATA_FORMAT_DB_ANGLE;
	    }
	    break;

	case 'm':
	    if (strncmp(cur, "ma", 2) == 0) {
		cur += 2;
		vfdp->vfd_format = VNADATA_FORMAT_MAG_ANGLE;
	    }
	    break;

	case 'r':
	    if (strncmp(cur, "ri", 2) == 0) {
		cur += 2;
		vfdp->vfd_format = VNADATA_FORMAT_REAL_IMAG;
	    }
	    break;

	default:
	    break;
	}
	break;

    default:
	break;
    }
    if (*cur == '\000') {
	return 0;
    }
bad_format:
    return -1;
}

/*
 * vnadata_set_format: set the format string
 *   @vdp: a pointer to the vnadata_t structure
 *   @format: a comma-separated case-insensitive list of the following:
 *     [{S,Z,Y,T,H,G,A,B}][{ri,ma,dB}]
 *     {il,rl}
 *     zin[{ri,ma}]
 *     {prc,prl,src,srl}
 *     vswr
 */
int vnadata_set_format(vnadata_t *vdp, const char *format)
{
    vnadata_internal_t *vdip;
    vnadata_format_descriptor_t *vfdp_new = NULL;
    size_t length;
    char *format_copy = NULL;
    char *cur;
    int nfields = 0;
    int rc = -1;

    /*
     * Validate pointer.
     */
    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }

    /*
     * If format is NULL, clear the format.
     */
    if (format == NULL) {
	goto update;
    }
    length = strlen(format);

    /*
     * Make a copy of the format string with spaces removed and
     * upper-case converted to lower.  Count comma-separated
     * fields and convert all commas to NUL's.
     */
    if ((format_copy = malloc(length + 1)) == NULL) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return -1;
    }
    nfields = 1;
    cur = format_copy;
    for (const char *cp = format; *cp != '\000'; ++cp) {
	if (*cp > 0x7e) {
	    _vnadata_error(vdip, VNAERR_USAGE, "vnadata_set_format: "
		    "invalid char '\\%02x' in format", *cp);
	    goto out;
	}
	if (isspace(*cp)) {
	    continue;
	}
	if (*cp == ',') {
	    ++nfields;
	    *cur++ = '\000';
	    continue;
	}
	*cur++ = isupper(*cp) ? tolower(*cp) : *cp;
    }
    *cur = '\000';

    /*
     * Allocate a new format vector.
     */
    if ((vfdp_new = calloc(nfields,
		    sizeof(vnadata_format_descriptor_t))) == NULL) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	goto out;
    }

    /*
     * Parse the formats and fill in each element of the vector.
     */
    cur = format_copy;
    for (int i = 0;;) {
	if (parse_format(&vfdp_new[i], cur) == -1) {
	    _vnadata_error(vdip, VNAERR_USAGE,
		    "invalid format specifier: \"%s\"", cur);
	    goto out;
	}
	if (++i >= nfields) {
	    break;
	}
	cur += strlen(cur) + 1;
    }

    /*
     * Replace the current format vector.
     */
update:
    free((void *)vdip->vdi_format_vector);
    vdip->vdi_format_vector = vfdp_new;
    vfdp_new = NULL;
    vdip->vdi_format_count = nfields;

    /*
     * Update the format string.
     */
    if (_vnadata_update_format_string(vdip) == -1) {
	goto out;
    }
    rc = 0;

out:
    free((void *)vfdp_new);
    free((void *)format_copy);
    return rc;
}
