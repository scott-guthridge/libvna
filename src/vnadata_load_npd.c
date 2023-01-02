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

#define VNADATA_NO_BOUNDS_CHECK
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
 * npd_record_type_t: type of record returned from scan_line
 */
typedef enum npd_record_type {
    T_KVERSION,
    T_KROWS,		/* for compatibility */
    T_KCOLUMNS,		/* for compatibility */
    T_KPORTS,
    T_KFREQUENCIES,
    T_KPARAMETERS,
    T_KFPRECISION,
    T_KDPRECISION,
    T_KZ0,
    T_DATA,
    T_EOF
} npd_record_type_t;

/*
 * npd_scan_state_t: scanner state
 */
typedef struct npd_scan_state {
    vnadata_internal_t	       *nss_vdip;
    FILE		       *nss_fp;
    const char		       *nss_filename;
    bool			nss_start_of_line;
    int				nss_line;
    int				nss_char;
    npd_record_type_t		nss_record_type;
    size_t			nss_text_size;
    size_t			nss_text_allocation;
    size_t			nss_field_count;
    size_t			nss_field_allocation;
    char		       *nss_text;
    int		               *nss_fields;
} npd_scan_state_t;

/*
 * GET_CHAR: read the next character
 *   @nssp: scanner state
 */
#define GET_CHAR(nssp)	\
	((nssp)->nss_char = getc((nssp)->nss_fp))

/*
 * FIELD: get the field at the specified index
 *   @nssp: scanner state
 *   @index: field index
 */
#define FIELD(nssp, index) \
    (&(nssp)->nss_text[(nssp)->nss_fields[index]])

/*
 * add_char: add a character to nss_text
 *   @nssp: scanner state
 *   @c:    character to add
 */
static int add_char(npd_scan_state_t *nssp, char c)
{
    vnadata_internal_t *vdip = nssp->nss_vdip;

    if (nssp->nss_text_size >= nssp->nss_text_allocation) {
	size_t new_allocation = MAX(81, 2 * nssp->nss_text_allocation);
	char *cp;

	if ((cp = realloc(nssp->nss_text, new_allocation)) == NULL) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "realloc: %s", strerror(errno));
	    return -1;
	}
	nssp->nss_text = cp;
	nssp->nss_text_allocation = new_allocation;
    }
    nssp->nss_text[nssp->nss_text_size++] = c;
    return 0;
}

/*
 * start_field: start a new field
 *   @nssp: scanner state
 */
static int start_field(npd_scan_state_t *nssp)
{
    vnadata_internal_t *vdip = nssp->nss_vdip;

    if (nssp->nss_field_count >= nssp->nss_field_allocation) {
	size_t new_allocation = MAX(9, 2 * nssp->nss_field_allocation);
	int *ip;

	ip = realloc(nssp->nss_fields, new_allocation * sizeof(int));
	if (ip == NULL) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "realloc: %s", strerror(errno));
	    return -1;
	}
	nssp->nss_fields = ip;
	nssp->nss_field_allocation = new_allocation;
    }
    nssp->nss_fields[nssp->nss_field_count] = nssp->nss_text_size;

    return 0;
}

/*
 * end_field: end the current field
 *   @nssp: scanner state
 */
static void end_field(npd_scan_state_t *nssp)
{
    add_char(nssp, '\000');
    ++nssp->nss_field_count;
}

/*
 * scan_line: scan an input line
 *   @nssp: scanner state
 */
static int scan_line(npd_scan_state_t *nssp)
{
    vnadata_internal_t *vdip = nssp->nss_vdip;

    /*
     * Start a new line.
     */
    nssp->nss_text_size = 0;
    nssp->nss_field_count = 0;
    for (;;) {
	/*
	 * Handle delayed advancement to the next line.  We do this
	 * so that nss_line remains accurate in _vnadata_load_npd.
	 */
	if (nssp->nss_start_of_line) {
	    assert(nssp->nss_char == '\n');
	    GET_CHAR(nssp);
	    ++nssp->nss_line;
	    nssp->nss_start_of_line = false;
	}

	/*
	 * At EOF, return whatever we have (possibly zero fields).
	 */
	if (nssp->nss_char == EOF) {
	    break;
	}

	/*
	 * Handle newline.
	 */
	if (nssp->nss_char == '\n') {
	    nssp->nss_start_of_line = true;	/* leave the newline char */
	    if (nssp->nss_field_count == 0) {	/* skip blank lines */
		continue;
	    }
	    break;
	}

	/*
	 * Skip whitespace.
	 */
	if (isascii(nssp->nss_char) && isspace(nssp->nss_char) &&
		nssp->nss_char != '\n') {
	    do {
		GET_CHAR(nssp);
	    } while (isascii(nssp->nss_char) && isspace(nssp->nss_char) &&
		    nssp->nss_char != '\n');
	    continue;
	}

	/*
	 * Skip comments (unless they're really keywords)
	 */
	if (nssp->nss_char == '#') {
	    GET_CHAR(nssp);
	    if (nssp->nss_char != ':') {
		while (nssp->nss_char != '\n' && nssp->nss_char != EOF) {
		    GET_CHAR(nssp);
		}
		continue;
	    }
	    GET_CHAR(nssp);
	    if (!isalpha(nssp->nss_char)) {
		while (nssp->nss_char != '\n' && nssp->nss_char != EOF) {
		    GET_CHAR(nssp);
		}
		continue;
	    }
	    if (start_field(nssp) == -1) {
		return -1;
	    }
	    if (add_char(nssp, '#') == -1) {
		return -1;
	    }
	    if (add_char(nssp, ':') == -1) {
		return -1;
	    }
	    goto in_field;
	}

	/*
	 * Read a field, using separate cases for the separator of space
	 * which we interpret as any 7-bit ASCII whitespace except newline,
	 * and other characters.
	 */
	if (start_field(nssp) == -1) {
	    return -1;
	}
    in_field:
	while (!isascii(nssp->nss_char) || !isspace(nssp->nss_char)) {
	    if (nssp->nss_char == EOF) {
		break;
	    }
	    if (add_char(nssp, nssp->nss_char) == -1) {
		return -1;
	    }
	    GET_CHAR(nssp);
	}
	end_field(nssp);
    }

    /*
     * Determine the record type.
     */
    if (nssp->nss_field_count == 0) {
	nssp->nss_record_type = T_EOF;
	return 0;
    }
    if (FIELD(nssp, 0)[0] == '#') {
	assert(FIELD(nssp, 0)[1] == ':');
	switch (FIELD(nssp, 0)[2]) {
	case 'c':
	    if (strcmp(&FIELD(nssp, 0)[2], "columns") == 0) {
		nssp->nss_record_type = T_KCOLUMNS;
		return 0;
	    }
	    break;

	case 'd':
	    if (strcmp(&FIELD(nssp, 0)[2], "dprecision") == 0) {
		nssp->nss_record_type = T_KDPRECISION;
		return 0;
	    }
	    break;

	case 'f':
	    if (strcmp(&FIELD(nssp, 0)[2], "frequencies") == 0) {
		nssp->nss_record_type = T_KFREQUENCIES;
		return 0;
	    }
	    if (strcmp(&FIELD(nssp, 0)[2], "fprecision") == 0) {
		nssp->nss_record_type = T_KFPRECISION;
		return 0;
	    }
	    break;

	case 'p':
	    if (strcmp(&FIELD(nssp, 0)[2], "parameters") == 0) {
		/* special-case: join the parameter fields by comma */
		for (size_t s = 3; s < nssp->nss_field_count; ++s) {
		    FIELD(nssp, s)[-1] = ',';
		}
		nssp->nss_field_count = 2;
		nssp->nss_record_type = T_KPARAMETERS;
		return 0;
	    }
	    if (strcmp(&FIELD(nssp, 0)[2], "ports") == 0) {
		nssp->nss_record_type = T_KPORTS;
		return 0;
	    }
	    break;

	case 'r':
	    if (strcmp(&FIELD(nssp, 0)[2], "rows") == 0) {
		nssp->nss_record_type = T_KROWS;
		return 0;
	    }
	    break;

	case 'v':
	    if (strcmp(&FIELD(nssp, 0)[2], "version") == 0) {
		nssp->nss_record_type = T_KVERSION;
		return 0;
	    }
	    break;

	case 'z':
	    if (strcmp(&FIELD(nssp, 0)[2], "z0") == 0) {
		nssp->nss_record_type = T_KZ0;
		return 0;
	    }
	    break;

	default:
	    break;
	}
	_vnadata_error(vdip, VNAERR_SYNTAX,
		"%s (line %d) error: unrecognized keyword: %s",
		nssp->nss_filename, nssp->nss_line, FIELD(nssp, 0));
	return -1;
    }
    nssp->nss_record_type = T_DATA;
    return 0;
}

/*
 * convert_int: convert an integer
 *   @field: string to convert
 *   @value: address to receive result
 */
static bool convert_int(const char *field, int *value)
{
    char *end;

    *value = strtol(field, &end, 0);
    if (end == field) {
	return false;
    }
    while (isspace(*end))
	++end;
    return *end == '\000';
}

/*
 * convert_double: convert a double
 *   @field: string to convert
 *   @value: address to receive result
 */
static bool convert_double(const char *field, double *value)
{
    char *end;

    *value = strtod(field, &end);
    if (end == field) {
	return false;
    }
    while (isspace(*end))
	++end;
    return *end == '\000';
}

/*
 * expect_nnint_arg: expect a single non-negative argument
 *   @nssp:  scanner state
 *   @value: address to receive value
 */
static int expect_nnint_arg(npd_scan_state_t *nssp, int *value)
{
    vnadata_internal_t *vdip = nssp->nss_vdip;
    int temp;

    if (nssp->nss_field_count != 2) {
	_vnadata_error(vdip, VNAERR_SYNTAX,
		"%s (line %d) error: one argument expected after %s",
		nssp->nss_filename, nssp->nss_line,
		&FIELD(nssp, 0)[2]);
	return -1;
    }
    if (!convert_int(FIELD(nssp, 1), &temp) || temp < 0) {
	_vnadata_error(vdip, VNAERR_SYNTAX,
		"%s (line %d) error: non-negative integer expected after %s",
		nssp->nss_filename, nssp->nss_line,
		FIELD(nssp, 0));
	return -1;
    }
    *value = temp;
    return 0;
}

/*
 * _vnadata_load_npd: load matrix data in libvna NPD format
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @type: convert parameters to this type (any if VPT_UNDEF)
 */
int _vnadata_load_npd(vnadata_internal_t *vdip, FILE *fp, const char *filename)
{
    vnadata_t *vdp = &vdip->vdi_vd;
    npd_scan_state_t nss;
    int rows = -1;
    int columns = -1;
    int ports = -1;
    int frequencies = -1;
    int rv = -1;
    int n_fields = 1;
    int best_quality = 0;
    int best_field = -1;
    int best_drows = -1;
    int best_dcolumns = -1;
    vnadata_parameter_type_t best_type = VPT_UNDEF;
    int parameter_line = -1;
    bool fz0 = false;
    double complex *z0_vector = NULL;
    const vnadata_format_descriptor_t *best_vfdp = NULL;

    (void)memset((void *)&nss, 0, sizeof(nss));
    nss.nss_vdip		= vdip;
    nss.nss_fp			= fp;
    nss.nss_filename		= filename;
    nss.nss_start_of_line	= true;
    nss.nss_line		= 0;
    nss.nss_char		= '\n';
    nss.nss_text_size		= 0;
    nss.nss_text_allocation	= 0;
    nss.nss_field_count		= 0;
    nss.nss_field_allocation	= 0;
    nss.nss_text		= NULL;
    nss.nss_fields		= NULL;
    if (scan_line(&nss) == -1) {
	goto out;
    }
    for (;;) {
	switch (nss.nss_record_type) {
	case T_KVERSION:
	    if (nss.nss_field_count < 2) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"argument expected after %s",
			nss.nss_filename, nss.nss_line,
			FIELD(&nss, 0));
		goto out;
	    }
	    if (strcmp(FIELD(&nss, 1), "1.0") != 0) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"unsupported version %s",
			nss.nss_filename, nss.nss_line, FIELD(&nss, 0));
		goto out;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KPORTS:
	    if (ports != -1) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"redundant ports line", nss.nss_filename, nss.nss_line);
		goto out;
	    }
	    if (expect_nnint_arg(&nss, &ports) == -1) {
		goto out;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KROWS:
	    if (expect_nnint_arg(&nss, &rows) == -1) {
		goto out;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KCOLUMNS:
	    if (expect_nnint_arg(&nss, &columns) == -1) {
		goto out;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KFREQUENCIES:
	    if (expect_nnint_arg(&nss, &frequencies) == -1) {
		goto out;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KPARAMETERS:
	    if (nss.nss_field_count != 2) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"at least one argument expected after %s",
			nss.nss_filename, nss.nss_line,
			FIELD(&nss, 0));
		goto out;
	    }
	    if (vnadata_set_format(vdp, FIELD(&nss, 1)) == -1) {
		goto out;
	    }
	    parameter_line = nss.nss_line;
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KFPRECISION:
	    {
		int temp;

		if (expect_nnint_arg(&nss, &temp) == -1) {
		    goto out;
		}
		if (temp > VNADATA_MAX_PRECISION) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s may not exceed %d",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, 0),
			    VNADATA_MAX_PRECISION);
		    goto out;
		}
		vdip->vdi_fprecision = temp;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_KDPRECISION:
	    {
		int temp;

		if (expect_nnint_arg(&nss, &temp) == -1) {
		    goto out;
		}
		if (temp > VNADATA_MAX_PRECISION) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s may not exceed %d",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, 0),
			    VNADATA_MAX_PRECISION);
		    goto out;
		}
		vdip->vdi_dprecision = temp;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;


	case T_KZ0:
	    if (ports < 0) {	/* compatibility with earlier code */
		if (rows >= 0 && columns >= 0) {
		    if (rows != columns) {
			_vnadata_error(vdip, VNAERR_SYNTAX,
				"%s (line %d) error: rows and columns must be "
				"equal", nss.nss_filename, nss.nss_line);
			goto out;
		    }
		    ports = columns;
		}
	    }
	    if (ports < 0) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"ports must come before #:z0",
			nss.nss_filename, nss.nss_line);
		goto out;
	    }
	    if (nss.nss_field_count == 2 &&
		    strcasecmp(FIELD(&nss, 1), "PER-FREQUENCY") == 0) {
		fz0 = true;
		if (scan_line(&nss) == -1) {
		    goto out;
		}
		continue;
	    }
	    if (nss.nss_field_count != 1 + 2 * ports) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected %d fields after z0",
			nss.nss_filename, nss.nss_line,
			2 * ports);
		goto out;
	    }
	    if (z0_vector == NULL) {
		if ((z0_vector = calloc(ports,
				sizeof(double complex))) == NULL) {
		    _vnadata_error(vdip, VNAERR_SYSTEM,
			    "calloc: %s", strerror(errno));
		    goto out;
		}
	    }
	    for (int port = 0; port < ports; ++port) {
		double re = 0.0, im = 0.0;
		char *cp;

		if (!convert_double(FIELD(&nss, 1 + 2 * port), &re)) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s: expected a numeric argument",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, 1 + 2 * port));
		    goto out;
		}
		if ((cp = strrchr(FIELD(&nss, 2 + 2 * port), 'j')) != NULL) {
		    if (cp[1] == '\000') {
			*cp = '\000';
		    }
		}
		if (!convert_double(FIELD(&nss, 2 + 2 * port), &im)) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s: expected a numeric argument",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, 2 + 2 * port));
		    goto out;
		}
		z0_vector[port] = re + I * im;
	    }
	    if (scan_line(&nss) == -1) {
		goto out;
	    }
	    continue;

	case T_DATA:
	case T_EOF:
	    break;
	}
	break;
    }
    if (ports < 0) {	/* compatibility with earlier code */
	if (rows >= 0 && columns >= 0) {
	    if (rows != columns) {
		_vnadata_error(vdip, VNAERR_SYNTAX,
			"%s (line %d) error: rows and columns must be "
			"equal", nss.nss_filename, nss.nss_line);
		goto out;
	    }
	    ports = columns;
	}
    }
    if (ports < 0) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"required keyword #:ports missing",
		nss.nss_filename, nss.nss_line);
	goto out;
    }
    if (frequencies < 0) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"required keyword #:frequencies missing",
		nss.nss_filename, nss.nss_line);
	goto out;
    }
    if (parameter_line == -1) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"required keyword #:parameters missing",
		nss.nss_filename, nss.nss_line);
	goto out;
    }

    /*
     * If the system impedances are frequency-dependent, add in
     * the Z0 fields.
     */
    if (fz0) {
	n_fields += 2 * ports;
    }

    /*
     * Find the best parameter for the requested type.
     */
    for (int i = 0; i < vdip->vdi_format_count; ++i) {
	const vnadata_format_descriptor_t *vfdp = &vdip->vdi_format_vector[i];
	vnadata_parameter_type_t type;
	int drows = ports, dcolumns = ports;
	int fields = 2 * drows * dcolumns;
	int quality = 0;

	/*
	 * Validate the parameter type against the matrix dimensions
	 * and determine the parameter type.  Determine the dimensions
	 * of the data matrix and the number of fields.
	 */
	type = vfdp->vfd_parameter;
	switch (vfdp->vfd_parameter) {
	case VPT_UNDEF:
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "%s parameter with no type",
		    nss.nss_filename, parameter_line,
		    vnadata_get_type_name(vfdp->vfd_parameter));
	    goto out;

	case VPT_S:
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_IL:
	    case VNADATA_FORMAT_RL:
	    case VNADATA_FORMAT_VSWR:
		fields = ports;
		break;

	    default:
		break;
	    }
	    break;

	case VPT_Z:
	case VPT_Y:
	    break;

	case VPT_T:
	case VPT_U:
	case VPT_H:
	case VPT_G:
	case VPT_A:
	case VPT_B:
	    if (ports != 2) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"%s parameters require a 2x2 matrix",
			nss.nss_filename, nss.nss_line,
			_vnadata_format_to_name(vfdp));
		goto out;
	    }
	    break;

	case VPT_ZIN:
	    drows = 1;
	    dcolumns = ports;
	    fields = 2 * ports;
	    break;

	default:
	    abort();
	    /*NOTREACHED*/
	}

	/*
	 * Find the best parameter type and format.  Matrix types
	 * are always better than Zin.  Then we fine tune based on
	 * the amount of math we have to do convert the parameter.
	 */
	if (vfdp->vfd_parameter != VPT_ZIN) {
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		quality = 6;	/* no conversion */
		break;
	    case VNADATA_FORMAT_MAG_ANGLE:
		quality = 5;	/* complex trig */
		break;
	    case VNADATA_FORMAT_DB_ANGLE:
		quality = 4;	/* exponentiation and complex trig */
		break;
	    default:
		quality = 0;
		break;
	    }
	} else {
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_REAL_IMAG:
		quality = 3;	/* no conversion */
		break;
	    case VNADATA_FORMAT_PRC:
	    case VNADATA_FORMAT_PRL:
	    case VNADATA_FORMAT_SRC:
	    case VNADATA_FORMAT_SRL:
		quality = 2;	/* multiplication and division */
		break;
	    case VNADATA_FORMAT_MAG_ANGLE:
		quality = 1;	/* trigonometry */
		break;
	    default:
		quality = 0;
		break;
	    }
	}
	if (quality > best_quality) {
	    best_quality = quality;
	    best_vfdp = vfdp;
	    best_type = type;
	    best_drows = drows;
	    best_dcolumns = dcolumns;
	    best_field = n_fields;
	}
	n_fields += fields;
    }
    if (best_vfdp == NULL) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"file contains no parameter we can load",
		nss.nss_filename, nss.nss_line);
	goto out;
    }

    /*
     * Set-up the output matrix.
     */
    if (vnadata_init(vdp, best_type, best_drows, best_dcolumns,
		frequencies) == -1) {
	goto out;
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"vnadata_init: %s", strerror(errno));
    }
    if (z0_vector != NULL) {
	if (vnadata_set_z0_vector(vdp, z0_vector) == -1) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "vnadata_set_z0_vector: %s", strerror(errno));
	    goto out;
	}
    } else if (fz0) {
	if ((z0_vector = calloc(ports, sizeof(double complex))) == NULL) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    goto out;
	}
    }

    /*
     * For each frequency, process a data line.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f;

	if (nss.nss_record_type != T_DATA) {
	    if (nss.nss_record_type == T_EOF) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected %d data lines; found only %d",
			nss.nss_filename, nss.nss_line,
			frequencies, findex + 1);
		goto out;
	    }
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected a data line: found %s",
		    nss.nss_filename, nss.nss_line,
		    FIELD(&nss, 0));
	    goto out;
	}
	if (nss.nss_field_count != n_fields) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected %d fields; found %d",
		    nss.nss_filename, nss.nss_line,
		    n_fields, (int)nss.nss_field_count);
	    goto out;
	}
	if (!convert_double(FIELD(&nss, 0), &f)) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "%s: number expected",
		    nss.nss_filename, nss.nss_line,
		    FIELD(&nss, 0));
	    goto out;
	}
	if (vnadata_set_frequency(vdp, findex, f) == -1) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "vnadata_set_frequency: %s", strerror(errno));
	    goto out;
	}
	if (fz0) {
	    for (int port = 0; port < ports; ++port) {
		double re, im;

		if (!convert_double(FIELD(&nss, 1 + 2 * port), &re)) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s: number expected",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, 1 + 2 * port));
		    goto out;
		}
		if (!convert_double(FIELD(&nss, 2 + 2 * port), &im)) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s: number expected",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, 2 + 2 * port));
		    goto out;
		}
		z0_vector[port] = re + I * im;
	    }
	    if (vnadata_set_fz0_vector(vdp, findex, z0_vector) == -1) {
		_vnadata_error(vdip, VNAERR_SYSTEM,
			"vnadata_set_fz0_vector: %s", strerror(errno));
		goto out;
	    }
	}
	for (int row = 0; row < best_drows; ++row) {
	    for (int column = 0; column < best_dcolumns; ++column) {
		int cell = row * best_dcolumns + column;
		double v1, v2;
		double complex value;

		if (!convert_double(FIELD(&nss,
				best_field + 2 * cell), &v1)) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s: number expected",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, best_field + cell));
		    goto out;
		}
		if (!convert_double(FIELD(&nss,
				best_field + 2 * cell + 1), &v2)) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "%s: number expected",
			    nss.nss_filename, nss.nss_line,
			    FIELD(&nss, best_field + cell + 1));
		    goto out;
		}
		switch (best_vfdp->vfd_format) {
		case VNADATA_FORMAT_DB_ANGLE:
		    value = pow(10.0, v1 / 20.0) * cexp(I * M_PI / 180.0 * v2);
		    break;

		case VNADATA_FORMAT_MAG_ANGLE:
		    value = v1 * cexp(I * M_PI / 180.0 * v2);
		    break;

		case VNADATA_FORMAT_REAL_IMAG:
		    value = v1 + I * v2;
		    break;

		case VNADATA_FORMAT_PRC:
		    value = 1.0 / (1.0 / v1 + 2.0 * M_PI * I * f * v2);
		    break;

		case VNADATA_FORMAT_PRL:
		    value = v1 / (1.0 - I * v1 / (2.0 * M_PI * f * v2));
		    break;

		case VNADATA_FORMAT_SRC:
		    value = v1 - I / (2.0 * M_PI * f * v2);
		    break;

		case VNADATA_FORMAT_SRL:
		    value = v1 + 2.0 * M_PI * I * f * v2;
		    break;

		default:
		    abort();
		    /*NOTREACHED*/
		}
		(void)vnadata_set_cell(vdp, findex, row, column, value);
	    }
	}
	if (scan_line(&nss) == -1) {
	    goto out;
	}
    }
    if (nss.nss_record_type != T_EOF) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"extra lines at end of input",
		nss.nss_filename, nss.nss_line);
	goto out;
    }
    rv = 0;

out:
    free((void *)z0_vector);
    free((void *)nss.nss_fields);
    free((void *)nss.nss_text);
    return rv;
}
