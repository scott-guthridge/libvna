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
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"

/*
 * _vnadata_set_name_from_filename: provide a default name from filename
 *   @vdip:   internal parameter matrix
 *   @filename: filename from which to derive the vnadata name
 */
void _vnadata_set_name_from_filename(vnadata_internal_t *vdip,
	const char *filename)
{
    /*
     * Set only if a name hasn't been explicitly set.
     */
    if (!(vdip->vdi_flags & VF_NAME_SET)) {
	const char *start;
	const char *cp;
	size_t n = VNADATA_MAX_NAME;

	/*
	 * Strip any path prefix.
	 */
	if ((start = strrchr(filename, '/')) != NULL) {
	    ++start;
	} else {
	    start = filename;
	}
	/*
	 * Strip extension.
	 */
	if ((cp = strrchr(start, '.')) != NULL) {
	    if (cp - start < n) {
		n = cp - start;
	    }
	}
	(void)strncpy(vdip->vdi_name, start, n);
	vdip->vdi_name[n] = '\000';
	vdip->vdi_flags |= VF_FILENAME_SEEN;
    }
}

/*
 * vnadata_load_common: load network parameters from a file
 *   @vdip:   internal parameter matrix
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 */
static int vnadata_load_common(vnadata_internal_t *vdip,
	FILE *fp, const char *filename)
{
    vnadata_t *vdp = &vdip->vdi_vd;
    int filename_ports = -1;

    /*
     * Determine the filetype and call the appropriate parser.  If the
     * filetype can be determined from the filename, use that type.
     * Otherwise, if the vnadata_t structure already has a filetype,
     * keep existing type.  If all else fails, default to the native
     * space-separated field NPD type.
     */
    {
	vnadata_filetype_t filetype;

	filetype = _vnadata_parse_filename(filename, &filename_ports);
	if (filetype != VNADATA_FILETYPE_AUTO) {
	    vdip->vdi_filetype = filetype;
	} else if (vdip->vdi_filetype == VNADATA_FILETYPE_AUTO) {
	    vdip->vdi_filetype = VNADATA_FILETYPE_NPD;
	}
    }
    switch (vdip->vdi_filetype) {
    case VNADATA_FILETYPE_TOUCHSTONE1:
    case VNADATA_FILETYPE_TOUCHSTONE2:
	/*
	 * Load touchstone formats.  Warn if filename_ports doesn't
	 * match actual number of ports.
	 */
	if (_vnadata_load_touchstone(vdip, fp, filename) == -1) {
	    return -1;
	}
	if (filename_ports != -1 && filename_ports != vdp->vd_columns) {
	    _vnadata_error(vdip, VNAERR_WARNING,
		    "%s: warning: filename suggests %d port(s) but found %d",
		    filename, filename_ports, vdp->vd_columns);
	}
	break;

    case VNADATA_FILETYPE_NPD:
	/*
	 * Load native NPD format.
	 */
	if (_vnadata_load_npd(vdip, fp, filename) == -1) {
	    return -1;
	}
	break;

    default:
	abort();
	/*NOTREACHED*/
    }
    _vnadata_set_name_from_filename(vdip, filename);
    return 0;
}

/*
 * vnadata_load: load network parameters from filename
 *   @vdp: a pointer to the vnadata_t structure
 *   @filename: file to load
 */
int vnadata_load(vnadata_t *vdp, const char *filename)
{
    vnadata_internal_t *vdip;
    FILE *fp;
    int rv;

    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    vdip = VDP_TO_VDIP(vdp);
    if (vdip->vdi_magic != VDI_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if ((fp = fopen(filename, "r")) == NULL) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"fopen: %s: %s", filename, strerror(errno));
	return -1;
    }
    rv = vnadata_load_common(vdip, fp, filename);
    (void)fclose(fp);
    if (rv == -1) {
	return -1;
    }
    return 0;
}

/*
 * vnadata_load: load network parameters from a file pointer
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 */
int vnadata_fload(vnadata_t *vdp, FILE *fp, const char *filename)
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
    if (vnadata_load_common(vdip, fp, filename) == -1) {
	return -1;
    }
    return 0;
}
