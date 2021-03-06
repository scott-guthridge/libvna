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

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "vnafile_internal.h"


/*
 * _vnafile_load_common: load network parameters from a file
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @vdp: output data (reshaped as needed)
 */
static int _vnafile_load_common(vnafile_t *vfp,
	FILE *fp, const char *filename, vnadata_t *vdp)
{
    int filename_ports = -1;

    if (vfp->vf_type != VNAFILE_NATIVE) {
	vnafile_type_t type;

	type = _vnafile_find_type(filename, &filename_ports);
	if (vfp->vf_type == VNAFILE_AUTO) {
	    vfp->vf_type = type;
	}
    }
    switch (vfp->vf_type) {
    case VNAFILE_TOUCHSTONE1:
    case VNAFILE_TOUCHSTONE2:
	/*
	 * Load touchstone formats.  Warn if filename_ports doesn't
	 * match actual ports, except allow .s2p as wild.
	 */
	if (_vnafile_load_touchstone(vfp, fp, filename, vdp) == -1) {
	    return -1;
	}
	if (filename_ports != -1 && filename_ports != 2 &&
		filename_ports != vdp->vd_columns) {
	    _vnafile_error(vfp, VNAERR_WARNING,
		    "%s: warning: filename suggests %d port(s) but %d found",
		    filename, filename_ports, vdp->vd_columns);
	}
	break;

    case VNAFILE_NATIVE:
	/*
	 * Load native format.
	 */
	if (_vnafile_load_native(vfp, fp, filename, vdp) == -1) {
	    return -1;
	}
	break;

    default:
	abort();
	/*NOTREACHED*/
    }
    return 0;
}

/*
 * vnafile_load: load network parameters from filename
 *   @filename: file to load
 *   @type: file format (use VNAFILE_AUTO to determine from filename)
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional opaque argument passed through to the error function
 *   @vdp: output data (reshaped as needed)
 */
vnafile_t *vnafile_load(const char *filename, vnafile_type_t type,
	vnaerr_error_fn_t *error_fn, void *error_arg, vnadata_t *vdp)
{
    vnafile_t *vfp;
    FILE *fp;
    int rv;

    if ((vfp = vnafile_alloc(error_fn, error_arg)) == NULL) {
	return NULL;
    }
    if ((fp = fopen(filename, "r")) == NULL) {
	_vnafile_error(vfp, VNAERR_SYSTEM,
		"fopen: %s: %s", filename, strerror(errno));
	vnafile_free(vfp);
	return NULL;
    }
    rv = _vnafile_load_common(vfp, fp, filename, vdp);
    (void)fclose(fp);
    if (rv == -1) {
	vnafile_free(vfp);
	return NULL;
    }
    return vfp;
}

/*
 * vnafile_load: load network parameters from a file pointer
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional opaque argument passed through to the error function
 *   @vdp: output data (reshaped as needed)
 */
vnafile_t *vnafile_fload(FILE *fp, const char *filename, vnafile_type_t type,
	vnaerr_error_fn_t *error_fn, void *error_arg, vnadata_t *vdp)
{
    vnafile_t *vfp;

    if ((vfp = vnafile_alloc(error_fn, error_arg)) == NULL) {
	return NULL;
    }
    if (_vnafile_load_common(vfp, fp, filename, vdp) == -1) {
	vnafile_free(vfp);
	return NULL;
    }
    return vfp;
}
