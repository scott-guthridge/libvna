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
#include <complex.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "vnacal_internal.h"
#include "vnaproperty_internal.h"


/*
 * Version Codes
 */
typedef enum vncal_version {
    V_UNSUPPORTED = -1,
    V0_2,
    V1_0
} vnacal_version_t;

/*
 * version_t: version table entry
 */
typedef struct version {
    int major;
    int minor;
    vnacal_version_t version;
} version_t;

/*
 * Version Table
 */
const version_t version_table[] = {
    {  1,  0, V1_0 },
    {  0,  2, V0_2 },
    { -1, -1, V_UNSUPPORTED },
};

/*
 * parse_version: parse the version line and return the verson code
 *   @vcp: vnacal structure
 *   @version_line: .vnacal version line
 */
static vnacal_version_t parse_version(vnacal_t *vcp, const char *version_line)
{
    int major_version, minor_version;

    if (sscanf(version_line, "#VNACal %d.%d",
		&major_version, &minor_version) != 2) {
	int old_major, old_minor;

	/*
	 * We maintain compatibility with two older versions of the
	 * calibration file.  Before libvna version 1.0, the #VNACal line
	 * was all upper case and there were two versions supported: 2.x
	 * and 3.x.  The old 2.x supports only E12 terms and stores them
	 * in a different format.  The old 3.x is the same as the new 1.0.
	 */
	if (sscanf(version_line, "#VNACAL %d.%d",
		    &old_major, &old_minor) != 2) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line 1) error: "
		    "expected #VNACal <major>.<minor>", vcp->vc_filename);
	    return V_UNSUPPORTED;
	}
	if (old_major == 2) {
	    major_version = 0;	/* map old 2.x to 0.2 */
	    minor_version = 2;
	} else if (old_major == 3) {
	    major_version = 1;	/* map old 3.x to 1.0 */
	    minor_version = 0;
	} else {
	    _vnacal_error(vcp, VNAERR_VERSION, "%s (line 1) error: "
		    "unsupported pre-release vnacal file version %d.%d",
		    vcp->vc_filename, old_major, old_minor);
	    return V_UNSUPPORTED;
	}
    }
    for (const version_t *vp = version_table; vp->major != -1; ++vp) {
	if (vp->major == major_version && vp->minor == minor_version)
	    return vp->version;
    }
    _vnacal_error(vcp, VNAERR_VERSION, "%s (line 1) error: "
	    "unsupported vnacal file version %d.%d",
	    vcp->vc_filename, major_version, minor_version);
    return -1;
}

/*
 * get_line: get the line number associated with node (+1 for version line)
 *   @node: a node in the property tree
 */
static int get_line(const vnaproperty_t *node)
{
    return _vnaproperty_get_line(node) + 1;
}

/*
 * get_key: get a required key and check the type
 *   @vcp: vnacal structure
 *   @mapping: a mapping node
 *   @key: requested key
 *   @required_type: 'm' (mapping), 'l' (sequence), 's' (scalar), or -1 (null)
 */
static const vnaproperty_t *get_key(vnacal_t *vcp,
	const vnaproperty_t *mapping, const char *key, int required_type)
{
    const vnaproperty_t *vprp;

    if ((vprp = vnaproperty_get_subtree(mapping, "%s", key)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "vnacal_load: %s (line %d): "
		"missing required key %s",
		vcp->vc_filename, get_line(mapping), key);
	return NULL;
    }
    if (vnaproperty_type(vprp, ".") != required_type) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "vnacal_load: %s (line %d): "
		"\"%s\" must have type %s",
		vcp->vc_filename, get_line(vprp), key,
		required_type == 'm' ? "mapping" :
		required_type == 'l' ? "sequence" :
		required_type == 's' ? "scalar" : "null");
	return NULL;
    }
    return vprp;
}

/*
 * top_level_keys: valid keys in top-level mapping
 */
static const char *top_level_keys[] = {
    "calibrations",
    "properties",
    NULL
};

/*
 * calibration_keys: valid keys in each calibration
 */
static const char *calibration_keys[] = {
    "columns",
    "data",
    "frequencies",
    "name",
    "properties",
    "rows",
    "type",
    "z0",
    NULL
};

/*
 * error_terms_keys: valid error term keys
 */
static const char *frequency_keys[] = {
    "el", "em", "er", "f",
    "ti", "tm", "ts", "tx",
    "ui", "um", "us", "ux",
    NULL
};

/*
 * v0_2_top_level_keys: valid top-level keys in version 0.2
 */
static const char *v0_2_top_level_keys[] = {
    "sets",
    NULL
};

/*
 * v0_2_calibration_keys: valid keys in each calibration in version 0.2
 */
static const char *v0_2_calibration_keys[] = {
    "columns",
    "data",
    "frequencies",
    "name",
    "rows",
    "z0",
    NULL
};

/*
 * v0_2_frequency_keys: valid per-frequency keys in version 0.2
 */
static const char *v0_2_frequency_keys[] = {
    "e", "f",
    NULL
};

/*
 * cmp: string compare function for qsort
 *   @arg1: first string
 *   @arg2: second string
 */
static int cmp(const void *arg1, const void *arg2)
{
    return strcmp(*(const char **)arg1, *(const char **)arg2);
}

/*
 * check_mapping: check that all keys in a mapping are expected
 *   @vcp: vnacal structure
 *   @mapping: element expected to be a mapping
 *   @allowed_keys: sorted NULL-terminated vector of valid keys
 */
static int check_mapping(vnacal_t *vcp, const vnaproperty_t *mapping,
	const char *const *allowed_keys)
{
    const char **keys;
    const char *const *ptr1;
    const char *const *ptr2;
    int count;

    assert(vnaproperty_type(mapping, ".") == 'm');
    if ((count = vnaproperty_count(mapping, "{}")) == -1) {
	abort();
    }
    if ((keys = vnaproperty_keys(mapping, "{}")) == NULL) {
	if (errno == ENOMEM) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "malloc: %s", strerror(errno));
	    return -1;
	}
	abort();	/* internal error */
    }
    qsort((void *)keys, count, sizeof(char *), cmp);
    ptr1 = keys;
    ptr2 = allowed_keys;
    for (;;) {
	int rv;

	if (*ptr1 == NULL) {
	    free((void *)keys);
	    return 0;
	}
	if (*ptr2 == NULL || (rv = strcmp(*ptr1, *ptr2)) < 0) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "vnacal_load: %s (line %d): "
		"error: unexpected key: %s", vcp->vc_filename,
		get_line(mapping), *ptr1);
	    free((void *)keys);
	    return -1;
	}
	if (rv == 0) {
	    ++ptr1;
	}
	++ptr2;
    }
}

/*
 * Assign bits to each matrix.
 */
enum {
    TS = 0x0001,
    TI = 0x0002,
    TX = 0x0004,
    TM = 0x0008,

    UM = 0x0010,
    UI = 0x0020,
    UX = 0x0040,
    US = 0x0080,

    EL = 0x0100,
    ER = 0x0200,
    EM = 0x0400,
};

/*
 * matrix_names: names corresponding to each bit position
 */
const char *matrix_names[] = {
    "ts", "ti", "tx", "tm",
    "um", "ui", "ux", "us",
    "el", "er", "em",
};
#define N_MATRIX_NAMES	(sizeof(matrix_names) / (sizeof(const char *)))

/*
 * check_for_stray_matrices: check for extraneous error term matrices
 *   @calp: calibration structure
 *   @vprp_frequency: per-frequency entry of a calibration
 */
static int check_for_stray_matrices(const vnacal_calibration_t *calp,
	const vnaproperty_t *vprp_frequency)
{
    vnacal_t *vcp = calp->cal_vcp;
    uint32_t mask = 0;	/* mask of wanted keys */

    switch (calp->cal_type) {
    case VNACAL_TE10:
	mask |= EL;
	/*FALLTHROUGH*/
    case VNACAL_T8:
    case VNACAL_T16:
	mask |= TS|TI|TX|TM;
	break;

    case VNACAL_UE10:
    case VNACAL_UE14:
	mask |= EL;
	/*FALLTHROUGH*/
    case VNACAL_U8:
    case VNACAL_U16:
	mask |= UM|UI|UX|US;
	break;

    case VNACAL_E12:
	mask |= EL|ER|EM;
	break;

    default:
	abort();
    }
    for (int i; i < N_MATRIX_NAMES; ++i) {
	if (mask & (1U << i)) {
	    continue;
	}
	if (vnaproperty_get_subtree(vprp_frequency, "%s",
		    matrix_names[i]) != NULL) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "key \"%s\" is not expected here",
		    vcp->vc_filename, get_line(vprp_frequency),
		    matrix_names[i]);
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_complex: parse a complex from a vnaproperty node and descriptor
 *   @root: vnaproperty node
 *   @format: property descriptor format (printf-like)
 *   @...: optional arguments depending on format
 */
static double complex parse_complex(const vnaproperty_t *root,
	const char *format, ...)
{
    va_list ap;
    const char *cur;
    char *end;
    double value1 = 0.0, value2 = 0.0;
    int code = 0;

    va_start(ap, format);
    cur = vnaproperty_vget(root, format, ap);
    va_end(ap);
    if (cur == NULL) {
	return HUGE_VAL;
    }
    value1 = strtod(cur, &end);
    if (end != cur) {
	++code;
	cur = end;
	value2 = strtod(cur, &end);
	if (end != cur) {
	    ++code;
	    cur = end;
	}
    }
    while (*cur == ' ' || *cur == '\t' || *cur == '\n') {
	++cur;
    }
    if (*cur == '+') {
	code |= 8;
	++cur;
    } else if (*cur == '-') {
	code |= 16;
	++cur;
    }
    while (*cur == ' ' || *cur == '\t' || *cur == '\n') {
	++cur;
    }
    switch (*cur) {
    case 'I':
    case 'J':
    case 'i':
    case 'j':
	code |= 4;
	++cur;
	break;

    default:
	break;
    }
    while (*cur == ' ' || *cur == '\t' || *cur == '\n') {
	++cur;
    }
    if (*cur != '\000') {
	return HUGE_VAL;
    }
    switch (code) {
    case 1:	/* number */
	return value1;
    case 4:	/* j */
	return I;
    case 5:	/* number j */
	return value1 * I;
    case 6:	/* number number j */
	return value1 + value2 * I;
    case 12:	/* +j */
	return I;
    case 13:	/* number + j */
	return value1 + I;
    case 20:	/* -j */
	return -I;
    case 21:	/* number - j */
	return value1 - I;
    default:
	break;
    }
    return HUGE_VAL;
}

/*
 * parse_type_from_map: parse a required vnacal type from a mapping
 *   @vcp: vnacal structure
 *   @mapping: mapping to parse
 *   @key: required key
 */
static vnacal_type_t parse_type_from_map(vnacal_t *vcp,
	const vnaproperty_t *mapping, const char *key)
{
    const vnaproperty_t *scalar;
    const char *s;
    vnacal_type_t type;

    if ((scalar = get_key(vcp, mapping, key, 's')) == NULL) {
	return VNACAL_NOTYPE;
    }
    s = vnaproperty_get(scalar, ".");
    assert(s != NULL);
    if ((type = vnacal_name_to_type(s)) == -1) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"unknown calibration type: %s",
		vcp->vc_filename, get_line(scalar), s);
	return VNACAL_NOTYPE;
    }
    return type;
}

/*
 * parse_int_from_map: parse a required integer from a mapping
 *   @vcp: vnacal structure
 *   @mapping: mapping to parse
 *   @key: required key
 *   @min: minimum valid value for the int
 */
static int parse_int_from_map(vnacal_t *vcp, const vnaproperty_t *mapping,
	const char *key, int min)
{
    const vnaproperty_t *scalar;
    const char *s;
    char *e;
    int value;

    if ((scalar = get_key(vcp, mapping, key, 's')) == NULL) {
	return -1;
    }
    s = vnaproperty_get(scalar, ".");
    assert(s != NULL);
    value = strtol(s, &e, 0);
    if (*s == '\000' || *e != '\000') {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"%s: invalid integer: \"%s\"",
		vcp->vc_filename, get_line(scalar), key, s);
	return -1;
    }
    if (value < min) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"%s must be at least %d (found %d)",
		vcp->vc_filename, get_line(scalar),
		key, min, value);
	return -1;
    }
    return value;
}

/*
 * parse_double_from_map: parse a required double from a mapping
 *   @vcp: vnacal structure
 *   @mapping: mapping to parse
 *   @key: required key
 */
static double parse_double_from_map(vnacal_t *vcp,
	const vnaproperty_t *mapping, const char *key)
{
    const vnaproperty_t *scalar;
    const char *s;
    char *e;
    double value;

    if ((scalar = get_key(vcp, mapping, key, 's')) == NULL) {
	return HUGE_VAL;
    }
    s = vnaproperty_get(scalar, ".");
    assert(s != NULL);
    value = strtod(s, &e);
    if (*s == '\000' || *e != '\000') {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"%s: invalid floating point number: \"%s\"",
		vcp->vc_filename, get_line(scalar), key, s);
	return HUGE_VAL;
    }
    return value;
}

/*
 * parse_double_from_map: parse a required complex from a mapping
 *   @vcp: vnacal structure
 *   @mapping: mapping to parse
 *   @key: required key
 */
static double complex parse_complex_from_map(vnacal_t *vcp,
	const vnaproperty_t *mapping, const char *key)
{
    const vnaproperty_t *scalar;
    double complex value;

    if ((scalar = get_key(vcp, mapping, key, 's')) == NULL) {
	return HUGE_VAL;
    }
    if ((value = parse_complex(scalar, ".")) == HUGE_VAL) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"%s: invalid complex number: \"%s\"",
		vcp->vc_filename, get_line(scalar), key,
		vnaproperty_get(scalar, "."));
	return HUGE_VAL;
    }
    return value;
}

/*
 * parse_frequency_entry_v0_2: parse a single frequency entry in v0.2
 *   @calp: pointer to calibration structure
 *   @vprp_frequency: per-frequency mapping to parse
 *   @error_terms: array of El, Er, Em matrices, each with [cell][findex]
 *   @findex: frequency index
 */
static int parse_frequency_entry_v0_2(vnacal_calibration_t *calp,
	const vnaproperty_t *vprp_frequency,
	double complex ***error_terms, int findex)
{
    vnacal_t *vcp = calp->cal_vcp;
    const vnaproperty_t *vprp_error_terms;
    const int rows = calp->cal_rows;
    const int columns = calp->cal_columns;
    double f;
    int count;

    assert(vnaproperty_type(vprp_frequency, ".") == 'm');
    if (check_mapping(vcp, vprp_frequency, v0_2_frequency_keys) == -1) {
	return -1;
    }
    if ((f = parse_double_from_map(vcp, vprp_frequency, "f")) == HUGE_VAL) {
	return -1;
    }
    calp->cal_frequency_vector[findex] = f;
    if ((vprp_error_terms = get_key(vcp, vprp_frequency, "e", 'l')) == NULL) {
	return -1;
    }
    count = vnaproperty_count(vprp_frequency, "e[]");
    if (count != rows) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"expected %d rows but found %d",
		vcp->vc_filename, get_line(vprp_error_terms),
		rows, count);
	return -1;
    }
    for (int row = 0; row < rows; ++row) {
	vnaproperty_t *vprp_row;

	vprp_row = vnaproperty_get_subtree(vprp_error_terms, "[%d][]", row);
	if (vprp_row == NULL) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "row %d of matrix must be a sequence",
		    vcp->vc_filename, get_line(vprp_error_terms),
		    row);
	    return -1;
	}
	count = vnaproperty_count(vprp_row, "[]");
	if (count != columns) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected row %d of matrix to have %d columns "
		    "but found %d",
		    vcp->vc_filename, get_line(vprp_error_terms),
		    row, columns, count);
	    return -1;
	}
	for (int column = 0; column < columns; ++column) {
	    const vnaproperty_t *vprp_terms;
	    const int cell = row * columns + column;

	    vprp_terms = vnaproperty_get_subtree(vprp_row, "[%d][]", column);
	    if (vprp_terms == NULL ||
		    (count = vnaproperty_count(vprp_terms, "[]")) != 3) {
		_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
			"matrix[%d][%d] must be a sequence of 3 error terms",
			vcp->vc_filename, get_line(vprp_terms),
			row, column);
		return -1;
	    }
	    for (int term = 0; term < 3; ++term) {
		double complex clf;

		clf = parse_complex(vprp_terms, "[%d]", term);
		if (clf == HUGE_VAL) {
		    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
			    "invalid complex number at matrix[%d][%d][%d]",
			    vcp->vc_filename, get_line(vprp_terms),
			    row, column, term);
		    return -1;
		}
		error_terms[term][cell][findex] = clf;
	    }
	}
    }
    return 0;
}

/*
 * parse_error_term_matrix: parse a single error term vector or matrix
 *   @vcp: vnacal structure
 *   @vprp_matrix: matrix to parse
 *   @vetmp: description of matrix to parse
 *   @findex: frequency index
 */
static int parse_error_term_matrix(const vnaproperty_t *vprp_matrix,
	vnacal_error_term_matrix_t *vetmp, int findex)
{
    vnacal_calibration_t *calp = vetmp->vetm_calp;
    vnacal_t *vcp = calp->cal_vcp;
    int count;
    double complex **matrix = vetmp->vetm_matrix;
    const int rows = vetmp->vetm_rows;
    const int columns = vetmp->vetm_columns;

    assert(vnaproperty_type(vprp_matrix, ".") == 'l');
    count = vnaproperty_count(vprp_matrix, "[]");
    switch (vetmp->vetm_type) {
    case VETM_VECTOR:
	assert(rows == 1);
	if (count != columns) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected %s vector to have %d elements but found %d",
		    vcp->vc_filename, get_line(vprp_matrix),
		    vetmp->vetm_name, columns, count);
	    return -1;
	}
	for (int i = 0; i < vetmp->vetm_columns; ++i) {
	    if ((matrix[i][findex] = parse_complex(vprp_matrix,
			    "[%d]", i)) == HUGE_VAL) {
		_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
			"invalid complex number in %s vector",
			vcp->vc_filename, get_line(vprp_matrix),
			vetmp->vetm_name);
		return -1;
	    }
	}
	break;

    case VETM_MATRIX_ND:
    case VETM_MATRIX:
	if (count != rows) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected %s matrix to have %d rows but found %d",
		    vcp->vc_filename, get_line(vprp_matrix),
		    vetmp->vetm_name, rows, count);
	    return -1;
	}
	for (int row = 0; row < rows; ++row) {
	    const vnaproperty_t *vprp_row;

	    vprp_row = vnaproperty_get_subtree(vprp_matrix, "[%d][]", row);
	    if (vprp_row == NULL) {
		_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
			"row %d of %s matrix must be a sequence",
			vcp->vc_filename, get_line(vprp_matrix),
			row, vetmp->vetm_name);
		return -1;
	    }
	    count = vnaproperty_count(vprp_row, "[]");
	    if (count != columns) {
		_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected row %d of %s matrix to have %d columns "
			"but found %d",
			vcp->vc_filename, get_line(vprp_matrix),
			row, vetmp->vetm_name, columns, count);
		return -1;
	    }
	    for (int column = 0; column < columns; ++column) {
		if (row != column || vetmp->vetm_type != VETM_MATRIX_ND) {
		    double complex clf;

		    clf = parse_complex(vprp_row, "[%d]", column);
		    if (clf == HUGE_VAL) {
			_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
				"invalid complex number at matrix element "
				"%s[%d][%d]",
				vcp->vc_filename, get_line(vprp_row),
				vetmp->vetm_name, row, column);
			return -1;
		    }
		    (*matrix++)[findex] = clf;
		} else {
		    if (vnaproperty_get_subtree(vprp_row, "[%d]",
				column) != NULL) {
			_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
				"expected diagonal matrix element %s[%d][%d] "
				"to be null",
				vcp->vc_filename, get_line(vprp_row),
				vetmp->vetm_name, row, column);
			return -1;
		    }
		}
	    }
	}
	break;

    default:
	abort();
    }
    return 0;
}

/*
 * parse_frequency_entry: parse a frequency entry of a calibration
 *   @calp: calibration structure we're filling
 *   @vprp_frequency: frequency mapping to parse
 *   @matrix_list: list of error term matrix descriptors
 *   @findex: frequency index
 */
static int parse_frequency_entry(vnacal_calibration_t *calp,
	const vnaproperty_t *vprp_frequency,
	vnacal_error_term_matrix_t *matrix_list, int findex)
{
    vnacal_t *vcp = calp->cal_vcp;
    double f;

    assert(vnaproperty_type(vprp_frequency, ".") == 'm');
    if (check_mapping(vcp, vprp_frequency, frequency_keys) == -1) {
	return -1;
    }
    if (check_for_stray_matrices(calp, vprp_frequency) == -1) {
	return -1;
    }
    if ((f = parse_double_from_map(vcp, vprp_frequency, "f")) == HUGE_VAL) {
	return -1;
    }
    calp->cal_frequency_vector[findex] = f;
    for (vnacal_error_term_matrix_t *vetmp = matrix_list; vetmp != NULL;
	    vetmp = vetmp->vetm_next) {
	const vnaproperty_t *vprp_matrix;

	if ((vprp_matrix = get_key(vcp, vprp_frequency,
			vetmp->vetm_name, 'l')) == NULL) {
	    return -1;
	}
	if (parse_error_term_matrix(vprp_matrix, vetmp, findex) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_calibration: parse a single calibration entry
 *   @vcp: vnacal structure
 *   @vprp_calibration: calibration entry to parse
 *   @version: version code
 */
static int parse_calibration(vnacal_t *vcp,
	const vnaproperty_t *vprp_calibration,
	vnacal_version_t version)
{
    const char *name;
    vnacal_type_t type = VNACAL_NOTYPE;
    int rows, columns, frequencies;
    double complex z0;
    vnacal_layout_t vl;
    vnacal_calibration_t *calp = NULL;
    const vnaproperty_t *vprp_data;
    int count;
    vnacal_error_term_matrix_t *matrix_list = NULL;
    int rc = -1;

    assert(vnaproperty_type(vprp_calibration, ".") == 'm');
    if (check_mapping(vcp, vprp_calibration, version == V0_2 ?
		v0_2_calibration_keys : calibration_keys) == -1) {
	goto out;
    }
    if ((name = vnaproperty_get(vprp_calibration, "name")) == NULL) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"expected scalar \"name\"",
		vcp->vc_filename, get_line(vprp_calibration));
	goto out;
    }
    if (version != V0_2) {
	if ((type = parse_type_from_map(vcp, vprp_calibration, "type")) == -1) {
	    goto out;
	}
    } else {
	type = VNACAL_E12;
    }
    errno = 0;
    if ((rows = parse_int_from_map(vcp,
		    vprp_calibration, "rows", 1)) == -1) {
	goto out;
    }
    if ((columns = parse_int_from_map(vcp,
		    vprp_calibration, "columns", 1)) == -1) {
	goto out;
    }
    if ((frequencies = parse_int_from_map(vcp,
		    vprp_calibration, "frequencies", 0)) == -1) {
	goto out;
    }
    if ((z0 = parse_complex_from_map(vcp,
		    vprp_calibration, "z0")) == HUGE_VAL) {
	goto out;
    }
    _vnacal_layout(&vl, type, rows, columns);
    if ((calp = _vnacal_calibration_alloc(vcp, type, rows, columns,
		    frequencies, VL_ERROR_TERMS(&vl))) == NULL) {
	goto out;
    }
    if (version != V0_2) {
	if (vnaproperty_copy(&calp->cal_properties,
		vnaproperty_get_subtree(vprp_calibration,
		    "properties")) == -1) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "vnaproperty_copy: %s",
		    strerror(errno));
	    goto out;
	}
    }
    if ((vprp_data = get_key(vcp, vprp_calibration, "data", 'l')) == NULL) {
	goto out;
    }
    count = vnaproperty_count(vprp_data, "[]");
    if (count != frequencies) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		"expected %d frequency entries, but found %d",
		vcp->vc_filename, get_line(vprp_data),
		frequencies, count);
	goto out;
    }
    if (_vnacal_build_error_term_list(calp, &vl, &matrix_list) == -1) {
	goto out;
    }
    if (version == V0_2) {
	double complex **error_terms[3];
	static const char *error_term_names[3] = { "el", "er", "em" };

	count = 0;
	for (vnacal_error_term_matrix_t *vetmp = matrix_list;
		vetmp != NULL; vetmp = vetmp->vetm_next) {
	    assert(count < 3);
	    assert(vetmp->vetm_rows == rows);
	    assert(vetmp->vetm_columns == columns);
	    assert(strcmp(vetmp->vetm_name, error_term_names[count]) == 0);
	    error_terms[count++] = vetmp->vetm_matrix;
	}
	for (int findex = 0; findex < frequencies; ++findex) {
	    const vnaproperty_t *vprp_frequency;

	    vprp_frequency = vnaproperty_get_subtree(vprp_data, "[%d]", findex);
	    assert(vprp_frequency != NULL);
	    if (parse_frequency_entry_v0_2(calp, vprp_frequency,
			error_terms, findex) == -1) {
		goto out;
	    }
	}

    } else {
	for (int findex = 0; findex < frequencies; ++findex) {
	    const vnaproperty_t *vprp_frequency;

	    vprp_frequency = vnaproperty_get_subtree(vprp_data, "[%d]", findex);
	    assert(vprp_frequency != NULL);
	    if (parse_frequency_entry(calp, vprp_frequency,
			matrix_list, findex) == -1) {
		goto out;
	    }
	}
    }
    for (int findex = 1; findex < frequencies; ++findex) {
	if (calp->cal_frequency_vector[findex - 1] >=
	    calp->cal_frequency_vector[findex]) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "frequencies must be ascending",
		    vcp->vc_filename, get_line(vprp_calibration));
	    goto out;
	}
    }
    if (_vnacal_add_calibration_common("vnacal_load", vcp, calp, name) == -1) {
	goto out;
    }
    rc = 0;

out:
    _vnacal_free_error_term_matrices(&matrix_list);
    if (rc != 0) {
	_vnacal_calibration_free(calp);
    }
    return rc;
}

/*
 * vnacal_load: load the calibration from a file
 *   @pathname: calibration file name
 *   @error_fn: error reporting callback or NULL
 *   @error_arg: arbitrary argument passed through to error_fn or NULL
 *
 *   If error_fn is non-NULL, then vnacal_load and subsequent functions report
 *   error messages using error_fn before returning failure to the caller.
 */
vnacal_t *vnacal_load(const char *pathname,
	vnaerr_error_fn_t *error_fn, void *error_arg)
{
    vnacal_t *vcp = NULL;
    FILE *fp = NULL;
    vnacal_version_t version;
    vnaproperty_t *vprp_root = NULL;
    const vnaproperty_t *vprp_calibrations;
    const char *calibrations_name;
    int calibrations;
    char line_buf[81];

    /*
     * Allocate the vnacal_t structure.
     */
    if ((vcp = _vnacal_alloc("vnacal_load", error_fn, error_arg)) == NULL) {
	return NULL;
    }

    /*
     * Save the filename.
     */
    if ((vcp->vc_filename = strdup(pathname)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"strdup: %s", strerror(errno));
	goto error;
    }

    /*
     * Open the file and parse the version line.
     */
    if ((fp = fopen(pathname, "r")) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "fopen: %s: %s",
		vcp->vc_filename, strerror(errno));
	goto error;
    }
    if (fgets(line_buf, sizeof(line_buf), fp) == NULL) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line 1) error: "
		"expected #VNACal <major>.<minor>",
		vcp->vc_filename);
	goto error;
    }
    line_buf[sizeof(line_buf) - 1] = '\000';
    if ((version = parse_version(vcp, line_buf)) == -1) {
	goto error;
    }
    if (vnaproperty_import_yaml_from_file(&vprp_root, fp, pathname,
		error_fn, error_arg) == -1) {
	goto error;
    }
    fclose(fp);
    fp = NULL;
    if (vnaproperty_type(vprp_root, ".") != 'm') {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line 2) error: "
		"top-level object must be a mapping",
		vcp->vc_filename);
	goto error;
    }
    if (version == V0_2) {
	calibrations_name = "sets";
	if (check_mapping(vcp, vprp_root, v0_2_top_level_keys) == -1) {
	    goto error;
	}
    } else {
	calibrations_name = "calibrations";
	if (check_mapping(vcp, vprp_root, top_level_keys) == -1) {
	    goto error;
	}
	/*
	 * Copy global properties.  Treat lack of a properties line as
	 * the same as properties set to null.
	 */
	if (vnaproperty_copy(&vcp->vc_properties,
		vnaproperty_get_subtree(vprp_root, "properties")) == -1) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "vnaproperty_copy: %s",
		    strerror(errno));
	    goto error;
	}
    }

    /*
     * Get and iterate through the calibrations sequence.
     */
    if ((vprp_calibrations = get_key(vcp, vprp_root,
		    calibrations_name, 'l')) == NULL) {
	goto error;
    }
    calibrations = vnaproperty_count(vprp_calibrations, "[]");
    for (int calibration = 0; calibration < calibrations; ++calibration) {
	vnaproperty_t *vprp_calibration;

	vprp_calibration = vnaproperty_get_subtree(vprp_calibrations,
		"[%d]{}", calibration);
	if (vprp_calibration == NULL) {
	    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %d) error: "
		    "calibration[%d] must be a mapping",
		    vcp->vc_filename, get_line(vprp_calibrations),
		    calibration);
	    goto error;
	}
	if (parse_calibration(vcp, vprp_calibration, version) == -1) {
	    goto error;
	}
    }
    assert(fp == NULL);
    (void)vnaproperty_delete(&vprp_root, ".");
    return vcp;

error:
    if (fp != NULL) {
	(void)fclose(fp);
    }
    (void)vnaproperty_delete(&vprp_root, ".");
    vnacal_free(vcp);
    return NULL;
}
