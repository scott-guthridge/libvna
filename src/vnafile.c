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

#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include "vnafile_internal.h"


/*
 * vnafile_error: report an error
 *   @vfp: pointer to vnafile_t
 *   @format: printf format string
 */
void _vnafile_error(const vnafile_t *vfp, const char *format, ...)
{
    va_list ap;
    char *msg = NULL;
    int rv;

    if (vfp->vf_error_fn != NULL) {
	va_start(ap, format);
	rv = vasprintf(&msg, format, ap);
	va_end(ap);
	if (rv == -1) {
	    (*vfp->vf_error_fn)(strerror(errno), vfp->vf_error_arg);
	    goto out;
	}
	(*vfp->vf_error_fn)(msg, vfp->vf_error_arg);
    }

out:
    free((void *)msg);
}

/*
 * parse_format: helper function for vnafile_set_format
 *   @vffp: format descriptor
 *   @format: format string
 */
static int parse_format(vnafile_format_t *vffp, const char *format)
{
    const char *cur = format;

    /*
     * Set defaults.
     */
    vffp->vff_parameter = VPT_UNDEF;
    vffp->vff_format = VNAFILE_FORMAT_REAL_IMAG;

    /*
     * Switch on the first letter of the parameter specifier.
     */
    switch (cur[0]) {
    case '\000':
	goto bad_format;

    case 'a':
	++cur;
	vffp->vff_parameter = VPT_A;
	goto parse_coordinates;

    case 'b':
	++cur;
	vffp->vff_parameter = VPT_B;
	goto parse_coordinates;

    case 'd':
	goto parse_coordinates;		/* db */

    case 'g':
	++cur;
	vffp->vff_parameter = VPT_G;
	goto parse_coordinates;

    case 'h':
	++cur;
	vffp->vff_parameter = VPT_H;
	goto parse_coordinates;

    case 'i':
	if (cur[1] == 'l') {
	    cur += 2;
	    vffp->vff_parameter = VPT_S;
	    vffp->vff_format = VNAFILE_FORMAT_IL;
	    break;
	}
	break;

    case 'm':
	goto parse_coordinates;		/* ma */

    case 'p':
	if (strncmp(cur, "prc", 3) == 0) {
	    cur += 3;
	    vffp->vff_parameter = VPT_ZIN;
	    vffp->vff_format = VNAFILE_FORMAT_PRC;
	    break;
	}
	if (strncmp(cur, "prl", 3) == 0) {
	    cur += 3;
	    vffp->vff_parameter = VPT_ZIN;
	    vffp->vff_format = VNAFILE_FORMAT_PRL;
	    break;
	}
	break;

    case 'r':
	if (cur[1] == 'l') {
	    cur += 2;
	    vffp->vff_parameter = VPT_S;
	    vffp->vff_format = VNAFILE_FORMAT_RL;
	    break;
	}
	goto parse_coordinates;		/* ri */

    case 's':
	if (strncmp(cur, "src", 3) == 0) {
	    cur += 3;
	    vffp->vff_parameter = VPT_ZIN;
	    vffp->vff_format = VNAFILE_FORMAT_SRC;
	    break;
	}
	if (strncmp(cur, "srl", 3) == 0) {
	    cur += 3;
	    vffp->vff_parameter = VPT_ZIN;
	    vffp->vff_format = VNAFILE_FORMAT_SRL;
	    break;
	}
	++cur;
	vffp->vff_parameter = VPT_S;
	goto parse_coordinates;

    case 't':
	++cur;
	vffp->vff_parameter = VPT_T;
	goto parse_coordinates;

    case 'v':
	if (strncmp(cur, "vswr", 4) == 0) {
	    cur += 4;
	    vffp->vff_parameter = VPT_S;
	    vffp->vff_format = VNAFILE_FORMAT_VSWR;
	    break;
	}
	break;

    case 'y':
	++cur;
	vffp->vff_parameter = VPT_Y;
	goto parse_coordinates;

    case 'z':
	if (strncmp(cur, "zin", 3) == 0) {
	    cur += 3;
	    vffp->vff_parameter = VPT_ZIN;
	    goto parse_coordinates;
	}
	++cur;
	vffp->vff_parameter = VPT_Z;
	goto parse_coordinates;

    parse_coordinates:
	switch (*cur) {
	case 'd':
	    if (strncmp(cur, "db", 2) == 0) {
		cur += 2;
		vffp->vff_format = VNAFILE_FORMAT_DB_ANGLE;
	    }
	    break;

	case 'm':
	    if (strncmp(cur, "ma", 2) == 0) {
		cur += 2;
		vffp->vff_format = VNAFILE_FORMAT_MAG_ANGLE;
	    }
	    break;

	case 'r':
	    if (strncmp(cur, "ri", 2) == 0) {
		cur += 2;
		vffp->vff_format = VNAFILE_FORMAT_REAL_IMAG;
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
 * MAX_FORMAT: maximum length of a format specifier
 */
#define MAX_FORMAT	4

/*
 * _vnafile_format_to_name: return the given format descriptor as a string
 *   @vffp: format descriptor
 */
const char *_vnafile_format_to_name(const vnafile_format_t *vffp)
{
    switch (vffp->vff_format) {
    case VNAFILE_FORMAT_PRC:
	return "PRC";

    case VNAFILE_FORMAT_PRL:
	return "PRL";

    case VNAFILE_FORMAT_SRC:
	return "SRC";

    case VNAFILE_FORMAT_SRL:
	return "SRL";

    case VNAFILE_FORMAT_IL:
	return "IL";

    case VNAFILE_FORMAT_RL:
	return "RL";

    case VNAFILE_FORMAT_VSWR:
	return "VSWR";

    default:
	switch (vffp->vff_parameter) {
	case VPT_UNDEF:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "ri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "ma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "dB";
	    default:
		break;
	    }
	    break;

	case VPT_S:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Sri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Sma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "SdB";
	    default:
		break;
	    }
	    break;

	case VPT_T:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Tri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Tma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "TdB";
	    default:
		break;
	    }
	    break;

	case VPT_Z:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Zri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Zma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "ZdB";
	    default:
		break;
	    }
	    break;

	case VPT_Y:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Yri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Yma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "YdB";
	    default:
		break;
	    }
	    break;

	case VPT_H:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Hri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Hma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "HdB";
	    default:
		break;
	    }
	    break;

	case VPT_G:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Gri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Gma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "GdB";
	    default:
		break;
	    }
	    break;

	case VPT_A:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Ari";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Ama";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "AdB";
	    default:
		break;
	    }
	    break;

	case VPT_B:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Bri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Bma";
	    case VNAFILE_FORMAT_DB_ANGLE:
		return "BdB";
	    default:
		break;
	    }
	    break;

	case VPT_ZIN:
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_REAL_IMAG:
		return "Zinri";
	    case VNAFILE_FORMAT_MAG_ANGLE:
		return "Zinma";
	    default:
		break;
	    }
	    break;

	default:
	    break;
	}
    }
    abort();
    /*NOTREACHED*/
}

/*
 * _vnafile_update_format_string: regenerate vfp_format_string
 *   @vfp: pointer to the structure returned from vnafile_alloc
 */
int _vnafile_update_format_string(vnafile_t *vfp)
{
    char *new_string = NULL;
    char *cur;

    free((void *)vfp->vf_format_string);
    vfp->vf_format_string = NULL;
    if ((new_string = malloc(vfp->vf_format_count *
		    (MAX_FORMAT + 1))) == NULL) {
	_vnafile_error(vfp, "malloc: %s", strerror(errno));
	return -1;
    }
    cur = new_string;
    for (int i = 0;;) {
	const char *name = _vnafile_format_to_name(&vfp->vf_format_vector[i]);

	(void)strcpy(cur, name);
	cur += strlen(cur);
	if (++i >= vfp->vf_format_count) {
	    break;
	}
	*cur++ = ',';
    }
    *cur = '\000';
    vfp->vf_format_string = new_string;
    return 0;
}

/*
 * vnafile_alloc: allocate the vnafile_t parameter structure
 *   @error_fn:  optional error reporting function
 *   @error_arg: optional opaque argument passed through to the error function
 */
vnafile_t *vnafile_alloc(vnafile_error_fn_t *error_fn, void *error_arg)
{
    vnafile_t *vfp;

    vfp = (vnafile_t *)malloc(sizeof(vnafile_t));
    if (vfp == NULL) {
	if (error_fn != NULL) {
	    char buf[80];

	    (void)snprintf(buf, sizeof(buf), "malloc: %s", strerror(errno));
            buf[sizeof(buf)-1] = '\000';
            (*error_fn)(error_arg, buf);
	}
	return NULL;
    }
    (void)memset((void *)vfp, 0, sizeof(vnafile_t));
    vfp->vf_magic      = VF_MAGIC;
    vfp->vf_error_fn   = error_fn;
    vfp->vf_error_arg  = error_arg;
    vfp->vf_type       = VNAFILE_AUTO;
    vfp->vf_fprecision = 7;
    vfp->vf_dprecision = 6;
    if ((vfp->vf_format_vector = malloc(sizeof(vnafile_format_t))) == NULL) {
	_vnafile_error(vfp, "malloc: %s", strerror(errno));
	vnafile_free(vfp);
	return NULL;
    }
    (void)memset((void *)vfp->vf_format_vector, 0, sizeof(vnafile_format_t));
    vfp->vf_format_vector[0].vff_parameter = VPT_UNDEF;
    vfp->vf_format_vector[0].vff_format = VNAFILE_FORMAT_REAL_IMAG;
    vfp->vf_format_count = 1;
    if (_vnafile_update_format_string(vfp) == -1) {
	vnafile_free(vfp);
	return NULL;
    }
    return vfp;
}

/*
 * vnafile_get_file_type: return the file type
 *   @vfp: pointer to the structure returned from vnafile_alloc
 */
vnafile_type_t vnafile_get_file_type(const vnafile_t *vfp)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return (vnafile_type_t)-1;
    }
    return vfp->vf_type;
}

/*
 * vnafile_set_file_type: set the file type
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @type: file type
 *
 *   The default type is VNAFILE_AUTO where the library tries to intuit
 *   the type from the filename.
 */
int vnafile_set_file_type(vnafile_t *vfp, vnafile_type_t type)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    switch (type) {
    case VNAFILE_AUTO:
    case VNAFILE_NATIVE:
    case VNAFILE_TOUCHSTONE1:
    case VNAFILE_TOUCHSTONE2:
	break;

    default:
	_vnafile_error(vfp, "vnafile_set_file_type: invalid type");
	errno = EINVAL;
	return -1;
    }
    vfp->vf_type = type;
    return 0;
}

/*
 * vnafile_get_format: current the format string
 *   @vfp: pointer to the structure returned from vnafile_alloc
 */
const char *vnafile_get_format(const vnafile_t *vfp)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    return vfp->vf_format_string;
}

/*
 * vnafile_set_format: set the format string
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @format: a comma-separated case-insensitive list of the following:
 *     [{S,Z,Y,T,H,G,A,B}][{ri,ma,dB}]
 *     {il,rl}
 *     zin[{ri,ma}]
 *     {prc,prl,src,srl}
 *     vswr
 *
 *   If format is NULL or empty string, default to "ri".
 */
int vnafile_set_format(vnafile_t *vfp, const char *format)
{
    vnafile_format_t *vffp_new = NULL;
    size_t length;
    char *format_copy = NULL;
    char *cur;
    int nfields = 1;
    int rc = -1;

    /*
     * Validate pointer.
     */
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }

    /*
     * If format is NULL, default it to "ri".
     */
    if (format == NULL) {
	format = "ri";
    }
    length = strlen(format);

    /*
     * Make a copy of the format string with spaces removed and
     * upper-case converted to lower.  Count comma-separated
     * fields and convert all commas to NUL's.
     */
    if ((format_copy = malloc(length + 1)) == NULL) {
	_vnafile_error(vfp, "malloc: %s", strerror(errno));
	return -1;
    }
    cur = format_copy;
    for (const char *cp = format; *cp != '\000'; ++cp) {
	if (*cp > 0x7e) {
	    _vnafile_error(vfp, "vnafile_set_format: "
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
    if ((vffp_new = calloc(nfields, sizeof(vnafile_format_t))) == NULL) {
	_vnafile_error(vfp, "malloc: %s", strerror(errno));
	goto out;
    }

    /*
     * Parse the formats and fill in each element of the vector.
     */
    cur = format_copy;
    for (int i = 0;;) {
	if (parse_format(&vffp_new[i], cur) == -1) {
	    _vnafile_error(vfp, "invalid format specifier: \"%s\"", cur);
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
    free((void *)vfp->vf_format_vector);
    vfp->vf_format_vector = vffp_new;
    vffp_new = NULL;
    vfp->vf_format_count = nfields;

    /*
     * Update the format string.
     */
    if (_vnafile_update_format_string(vfp) == -1) {
	goto out;
    }
    rc = 0;

out:
    free((void *)vffp_new);
    free((void *)format_copy);
    return rc;
}

/*
 * _vnafile_set_simple_format: set the format string (1 parameter)
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @parameter: parameter type
 *   @format: coordinate system
 */
int _vnafile_set_simple_format(vnafile_t *vfp,
	vnadata_parameter_type_t parameter, vnafile_format_type_t format)
{
    vnafile_format_t *vffp_new = NULL;
    int rc;

    /*
     * Allocate the new format vector, length 1.
     */
    if ((vffp_new = malloc(sizeof(vnafile_format_t))) == NULL) {
	_vnafile_error(vfp, "malloc: %s", strerror(errno));
	goto out;
    }
    (void)memset((void *)vffp_new, 0, sizeof(*vffp_new));
    vffp_new->vff_parameter = parameter;
    vffp_new->vff_format = format;

    /*
     * Install the new vector.
     */
    free((void *)vfp->vf_format_vector);
    vfp->vf_format_vector = vffp_new;
    vfp->vf_format_count = 1;

    /*
     * Update the format string.
     */
    if (_vnafile_update_format_string(vfp) == -1) {
	goto out;
    }
    rc = 0;

out:
    return rc;
}

/*
 * vnafile_get_fprecision: get the frequency value precision
 *   @vfp: pointer to the structure returned from vnafile_alloc
 */
int vnafile_get_fprecision(const vnafile_t *vfp)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    return vfp->vf_fprecision;
}

/*
 * vnafile_set_fprecision: set the frequency value precision
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @precision: precision in decimal places (1..n) or VNAFILE_MAX_PRECISION
 */
int vnafile_set_fprecision(vnafile_t *vfp, int precision)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (precision < 1) {
	_vnafile_error(vfp, "vnafile_set_fprecision: invalid precision: %s",
		precision);
	errno = EINVAL;
	return -1;
    }
    vfp->vf_fprecision = precision;
    return 0;
}

/*
 * vnafile_get_dprecision: set the data value precision
 *   @vfp: pointer to the structure returned from vnafile_alloc
 */
int vnafile_get_dprecision(const vnafile_t *vfp)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    return vfp->vf_dprecision;
}

/*
 * vnafile_set_dprecision: set the data value precision for vnafile_save
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @precision: precision in decimal places (1..n) or VNAFILE_MAX_PRECISION
 */
int vnafile_set_dprecision(vnafile_t *vfp, int precision)
{
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (precision < 1) {
	_vnafile_error(vfp, "vnafile_set_fprecision: invalid precision: %s",
		precision);
	errno = EINVAL;
	return -1;
    }
    vfp->vf_dprecision = precision;
    return 0;
}

/*
 * vnafile_free: free the structure obtained from vnafile_alloc
 *   @vfp: pointer to structure to free
 */
void vnafile_free(vnafile_t *vfp)
{
    if (vfp != NULL && vfp->vf_magic == VF_MAGIC) {
	vfp->vf_magic = -1;
	free((void *)vfp->vf_format_vector);
	free((void *)vfp->vf_format_string);
	free((void *)vfp);
    }
}
