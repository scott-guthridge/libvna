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
 * print_value: print a double in engineering form
 *   @fp:        output stream
 *   @precision: digits of precision to print
 *   @plus:	 include the plus sign
 *   @pad:       pad to consistent width
 *   @value:     value to print
 *
 * Return:
 *   The number of characters printed or -1 on error.
 */
static void print_value(FILE *fp, int precision, bool plus, bool pad,
	double value)
{
    char *cur;
    char *mantissa;
    int width, exponent, before, temp;
    char sign = '+';
    char buf1[MAX(precision, 1) + 8];
    char buf2[MAX(precision, 1) + 8];

    /*
     * Bound the minimum precision to 1 digit.  If precision is
     * VNAFILE_MAX_PRECISION, write hexadecimal floating point format.
     */
    if (precision < 1) {
	precision = 1;

    } else if (precision == VNAFILE_MAX_PRECISION) {
	fprintf(fp, "%a", value);
	return;
    }
    width = precision + 5; /* .e-EE */
    if (plus) {
	++width;
    }

    /*
     * Format the number and extract mantissa and exponent.  If the value
     * is nan or inf, use the value as-is.  Upon successfully reaching
     * the end of block, mantissa points to a string of all the digits,
     * decimal point removed, and the exponent is decoded in exponent.
     */
    (void)sprintf(buf1, "%.*e", precision - 1, value);
    cur = buf1;
    if (cur[0] == '-') {
	sign = '-';
	++cur;
    }
    if (!isdigit(*cur)) {
	(void)strcpy(buf2, buf1);
	goto finished;
    }
    mantissa = cur;
    if (cur[1] == '.') {
	cur[1] = cur[0];
	mantissa = &cur[1];
    }
    cur = strchr(mantissa, 'e');
    if (cur == NULL) {
	(void)strcpy(buf2, buf1);
	goto finished;
    }
    *cur++ = '\000';
    exponent = atoi(cur);

    /*
     * Add the sign, if needed.
     */
    cur = buf2;
    if (plus || sign == '-') {
	*cur++ = sign;
    }

    /*
     * Determine the number of digits to appear before the decial point
     * and adjust exponent.
     */
    switch (precision) {
    case 1:
	before = 1;
	break;

    case 2:
	temp = exponent + 1;
	before = ((temp >= 0) ? temp % 3 : 2 - (-temp - 1) % 3);
	assert(before >= 0 && before <= 2);
	break;

    default:
	temp = exponent;
	before = ((temp >= 0) ? temp % 3 : 2 - (-temp - 1) % 3) + 1;
	break;
    }
    exponent -= (before - 1);

    /*
     * Format the value and exponent and print.
     */
    (void)memcpy((void *)cur, (void *)mantissa, before);
    cur += before;
    if (precision - before > 0 || exponent == 0) {
	*cur++ = '.';
	(void)memcpy((void *)cur, (void *)&mantissa[before],
	    precision - before);
	cur += precision - before;
    }
    if (exponent != 0) {
	(void)sprintf(cur, "e%0+3d", exponent);
	cur += strlen(cur);
    } else if (pad) {
	(void)strcpy(cur, "    ");
	cur += 4;
    }
    *cur = '\000';

finished:
    if (pad) {
	fprintf(fp, "%-*s", width, buf2);
	return;
    }
    fprintf(fp, "%s", buf2);
}

/*
 * _vnafile_convert: convert the input matrix to the given type
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @conversions: cached conversions
 *   @vdp: input data
 *   @type: desired type
 *   @function: function name (for error messages)
 */
static int _vnafile_convert(vnafile_t *vfp,
	vnadata_t **conversions, const vnadata_t *vdp,
	vnadata_parameter_type_t type, const char *function)
{
    vnadata_t *target = NULL;

    if (conversions[type] != NULL) {
	return 0;
    }
    if (type == vnadata_get_type(vdp)) {
	return 0;
    }
    if ((target = vnadata_alloc()) == NULL) {
	_vnafile_error(vfp, VNAERR_SYSTEM,
		"vnadata_alloc: %s", strerror(errno));
	return -1;
    }
    if (vnadata_convert(vdp, target, type) == -1) {
	if (errno == EINVAL) {
	    _vnafile_error(vfp, VNAERR_USAGE,
		    "%s: cannot convert from %s to %s",
		function, vnadata_get_typename(vnadata_get_type(vdp)),
		vnadata_get_typename(type));
	    vnadata_free(target);
	    return -1;
	}
	_vnafile_error(vfp, VNAERR_SYSTEM,
		"vnadata_convert: %s", strerror(errno));
	vnadata_free(target);
	return -1;
    }
    conversions[type] = target;
    return 0;
}

/*
 * _vnafile_find_type: try to determine the file format from the filename
 *   @filename: input or output filename
 */
vnafile_type_t _vnafile_find_type(const char *filename, int *ports)
{
    const char *suffix;

    /*
     * Set *ports to the default of unknown (-1).
     */
    if (ports != NULL) {
	*ports = -1;
    }

    /*
     * If there's no file extension, default to native format.
     */
    if ((suffix = strrchr(filename, '.')) == NULL) {
	return VNAFILE_NATIVE;
    }
    ++suffix;

    /*
     * If the file has a suffix of .ts, assume touchstone 2.
     */
    if (strcasecmp(suffix, "ts") == 0) {
	return VNAFILE_TOUCHSTONE2;
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
		return VNAFILE_TOUCHSTONE1;
	    }
	}
    }

    /*
     * None of the above.  Default to native.
     */
    return VNAFILE_NATIVE;
}

/*
 * print_native_header: print header for native file format
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @fp: file pointer
 *   @vdp: input data
 */
static void print_native_header(vnafile_t *vfp, FILE *fp, const vnadata_t *vdp)
{
    int rows, columns, ports, diagonals;
    int output_fields = 1;
    int current_field = 0;
    int field_width;
    int port_width, port_pair_width;
    int parameter_width = 0;
    char parameter_buf[3 + 2 * 3 * sizeof(int) + 1];
    const double complex *z0_vector = vnadata_get_z0_vector(vdp);

    /*
     * Get the number of rows and columns.  We have to do a little
     * munging here if the vdp's parameter type is Zin because a Zin
     * vector is really a diagonal matrix stored as a vector.  We want
     * the dimensions of the matrix, not the vector here.
     */
    rows    = vnadata_get_rows(vdp);
    columns = vnadata_get_columns(vdp);
    if (vdp->vd_type == VPT_ZIN) {
	if (rows == 1) {
	    rows = columns;
	} else if (columns == 1) {
	    columns = rows;
	}
    }
    ports     = MAX(rows, columns);
    diagonals = MIN(rows, columns);

    /*
     * Determine the number of digits needed to display both a single
     * port number and a pair of port numbers.
     */
    if ((ports + 1) < 10) {
	port_width = 1;
    } else {
	port_width = (int)floor(log10(ports + 0.5)) + 1;
    }
    if (ports > 9) {
	port_pair_width = 2 * port_width + 1;
    } else {
	port_pair_width = 2 * port_width;
    }

    /*
     * Determine the number of output fields and from it determine the
     * number of digits needed to display the field number.  At the
     * same time find the number of digits needed for the parameter name
     * and ports.
     */
    if (z0_vector == NULL) {
	output_fields += 2 * ports;
    }
    for (int format = 0; format < vfp->vf_format_count; ++format) {
	const vnafile_format_t *vffp = &vfp->vf_format_vector[format];

	switch (vffp->vff_format) {
	case VNAFILE_FORMAT_DB_ANGLE:
	case VNAFILE_FORMAT_MAG_ANGLE:
	case VNAFILE_FORMAT_REAL_IMAG:
	    if (vffp->vff_parameter != VPT_ZIN) {
		output_fields += 2 * rows * columns;
		parameter_width = MAX(parameter_width, 1 + port_pair_width);
	    } else {
		output_fields += diagonals;
		parameter_width = MAX(parameter_width, 3 + port_width);
	    }
	    break;

	case VNAFILE_FORMAT_PRC:
	case VNAFILE_FORMAT_PRL:
	case VNAFILE_FORMAT_SRC:
	case VNAFILE_FORMAT_SRL:
	    output_fields += 2 * diagonals;
	    parameter_width = MAX(parameter_width, 3 + port_width);
	    break;

	case VNAFILE_FORMAT_IL:
	    output_fields += rows * columns - diagonals;
	    parameter_width = MAX(parameter_width, 2 + port_pair_width);
	    break;

	case VNAFILE_FORMAT_RL:
	    output_fields += diagonals;
	    parameter_width = MAX(parameter_width, 2 + port_width);
	    break;

	case VNAFILE_FORMAT_VSWR:
	    output_fields += diagonals;
	    parameter_width = MAX(parameter_width, 4 + port_width);
	    break;

	default:
	    abort();
	    /*NOTREACHED*/
	}
    }
    field_width = (int)floor(log10(output_fields + 0.5)) + 1;

    /*
     * Print the preamble.
     */
    (void)fprintf(fp, "# NPD\n");
    (void)fprintf(fp, "#:version 1.0\n");
    (void)fprintf(fp, "#:rows %d\n", rows);
    (void)fprintf(fp, "#:columns %d\n", columns);
    (void)fprintf(fp, "#:frequencies %d\n", vnadata_get_frequencies(vdp));
    (void)fprintf(fp, "#:parameters %s\n", vfp->vf_format_string);
    (void)fprintf(fp, "#:z0");
    if (z0_vector == NULL) {
	(void)fprintf(fp, " PER-FREQUENCY\n");
    } else {
	for (int port = 0; port < ports; ++port) {
	    double complex z0 = z0_vector[port];

	    (void)fputc(' ', fp);
	    print_value(fp, vfp->vf_dprecision, /*plus=*/false, /*pad=*/false,
		    creal(z0));
	    (void)fputc(' ', fp);
	    print_value(fp, vfp->vf_dprecision, /*plus=*/true, /*pad=*/false,
		    cimag(z0));
	    (void)fputc('j', fp);
	}
	(void)fputc('\n', fp);
    }
    (void)fprintf(fp, "#:fprecision %d\n", vfp->vf_fprecision);
    (void)fprintf(fp, "#:dprecision %d\n", vfp->vf_dprecision);
    (void)fprintf(fp, "#\n");

    /*
     * Print a key for each field.
     */
    (void)fprintf(fp, "# field %*d: %-*s (Hz)\n",
	    field_width, ++current_field, 10 + parameter_width, "frequency");
    if (z0_vector == NULL) {
	for (int port = 0; port < ports; ++port) {
	    (void)sprintf(parameter_buf, "Z%d", port + 1);
	    (void)fprintf(fp, "# field %*d: %-*s real      (ohms)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	    (void)fprintf(fp, "# field %*d: %-*s imaginary (ohms)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	}
    }
    for (int format = 0; format < vfp->vf_format_count; ++format) {
	const vnafile_format_t *vffp = &vfp->vf_format_vector[format];
	const char *name;
	const char *const *type_vector;
	static const char *st_types[] = {
	    "v-ratio", NULL
	};
	static const char *z_types[] = {
	    "ohms", NULL
	};
	static const char *y_types[] = {
	    "seimens", NULL
	};
	static const char *h_types[] = {
	    "ohms", "v-ratio", "i-ratio", "seimens", NULL
	};
	static const char *g_types[] = {
	    "seimens", "i-ratio", "v-ratio", "ohms", NULL
	};
	static const char *ab_types[] = {
	    "v-ratio", "ohms", "seimens", "i-ratio", NULL
	};
	switch (vffp->vff_format) {
	case VNAFILE_FORMAT_DB_ANGLE:
	case VNAFILE_FORMAT_MAG_ANGLE:
	case VNAFILE_FORMAT_REAL_IMAG:
	    /*
	     * Handle Zin.
	     */
	    if (vffp->vff_parameter == VPT_ZIN) {
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    (void)sprintf(parameter_buf, "Zin%d", diagonal + 1);
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vffp->vff_format) {
		    case VNAFILE_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " real      (ohms)\n");
			break;
		    case VNAFILE_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, " magnitude (ohms)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vffp->vff_format) {
		    case VNAFILE_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " imaginary (ohms)\n");
			break;
		    case VNAFILE_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, " angle     (degrees)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		}
		break;
	    }
	    /*
	     * Handle matrix types.
	     */
	    switch (vffp->vff_parameter) {
	    case VPT_S:
		name = "S";
		type_vector = st_types;
		break;

	    case VPT_Z:
		name = "Z";
		type_vector = z_types;
		break;

	    case VPT_Y:
		name = "Y";
		type_vector = y_types;
		break;

	    case VPT_T:
		name = "T";
		type_vector = st_types;
		break;

	    case VPT_H:
		name = "H";
		type_vector = h_types;
		break;

	    case VPT_G:
		name = "G";
		type_vector = g_types;
		break;

	    case VPT_A:
		name = "A";
		type_vector = ab_types;
		break;

	    case VPT_B:
		name = "B";
		type_vector = ab_types;
		break;

	    default:
		abort();
		/*NOTREACHED*/
	    }
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    const char *type = (type_vector[1] == NULL) ?
			type_vector[0] : type_vector[row * columns + column];

		    if (ports <= 9) {
			(void)sprintf(parameter_buf, "%s%d%d",
				name, row + 1, column + 1);
		    } else {
			(void)sprintf(parameter_buf, "%s%d,%d",
				name, row + 1, column + 1);
		    }
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vffp->vff_format) {
		    case VNAFILE_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " real      (%s)\n", type);
			break;
		    case VNAFILE_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, " magnitude (%s)\n", type);
			break;
		    case VNAFILE_FORMAT_DB_ANGLE:
			(void)fprintf(fp, " magnitude (dB)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vffp->vff_format) {
		    case VNAFILE_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " imaginary (%s)\n", type);
			break;
		    case VNAFILE_FORMAT_MAG_ANGLE:
		    case VNAFILE_FORMAT_DB_ANGLE:
			(void)fprintf(fp, " angle     (degrees)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		}
	    }
	    break;

	case VNAFILE_FORMAT_PRC:
	    assert(vffp->vff_parameter == VPT_ZIN);
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(parameter_buf, "PRC%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
                (void)fprintf(fp, "# field %*d: %-*s C         (farads)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
            }
            break;

	case VNAFILE_FORMAT_PRL:
	    assert(vffp->vff_parameter == VPT_ZIN);
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(parameter_buf, "PRL%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
                (void)fprintf(fp, "# field %*d: %-*s L         (henries)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
            }
            break;

	case VNAFILE_FORMAT_SRC:
	    assert(vffp->vff_parameter == VPT_ZIN);
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(parameter_buf, "SRC%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
                (void)fprintf(fp, "# field %*d: %-*s C         (farads)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
            }
            break;

	case VNAFILE_FORMAT_SRL:
	    assert(vffp->vff_parameter == VPT_ZIN);
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(parameter_buf, "SRL%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
                (void)fprintf(fp, "# field %*d: %-*s L         (henries)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
            }
            break;

	case VNAFILE_FORMAT_IL:
	    assert(vffp->vff_parameter == VPT_S);
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    if (row == column) {
			continue;
		    }
		    if (ports <= 9) {
			(void)sprintf(parameter_buf, "IL%d%d",
				row + 1, column + 1);
		    } else {
			(void)sprintf(parameter_buf, "IL%d,%d",
				row + 1, column + 1);
		    }
		    (void)fprintf(fp, "# field %*d: %-*s magnitude (dB)\n",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		}
	    }
            break;

	case VNAFILE_FORMAT_RL:
	    assert(vffp->vff_parameter == VPT_S);
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(parameter_buf, "RL%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: %-*s magnitude (dB)\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
            }
            break;

	case VNAFILE_FORMAT_VSWR:
	    assert(vffp->vff_parameter == VPT_S);
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(parameter_buf, "VSWR%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: %-*s\n",
                    field_width, ++current_field, parameter_width,
		    parameter_buf);
            }
            break;

	default:
	    abort();
	    /*NOTREACHED*/
	}
    }
    (void)fprintf(fp, "#\n");
}

/*
 * print_touchstone_header: print header for touchstone formats
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @fp: file pointer
 *   @vdp: input data
 *   @z0_touchstone: the first system impedance (before normalization)
 */
static void print_touchstone_header(vnafile_t *vfp, FILE *fp,
	const vnadata_t *vdp, double z0_touchstone)
{
    vnafile_format_t *vffp = &vfp->vf_format_vector[0];
    char parameter_name;
    const char *format;
    int ports = vnadata_get_rows(vdp);
    const double complex *z0_vector = vnadata_get_z0_vector(vdp);

    assert(vfp->vf_format_count == 1);
    assert(vfp->vf_type == VNAFILE_TOUCHSTONE1 ||
	   vfp->vf_type == VNAFILE_TOUCHSTONE2);
    if (vfp->vf_type == VNAFILE_TOUCHSTONE2) {
	(void)fprintf(fp, "[Version] 2.0\n");
    }
    switch (vffp->vff_parameter) {
    case VPT_S:
	parameter_name = 'S';
	break;
    case VPT_Z:
	parameter_name = 'Z';
	break;
    case VPT_Y:
	parameter_name = 'Y';
	break;
    case VPT_H:
	parameter_name = 'H';
	break;
    case VPT_G:
	parameter_name = 'G';
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    switch (vffp->vff_format) {
    case VNAFILE_FORMAT_DB_ANGLE:
	format = "DB";
	break;
    case VNAFILE_FORMAT_MAG_ANGLE:
	format = "MA";
	break;
    case VNAFILE_FORMAT_REAL_IMAG:
	format = "RI";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)fprintf(fp, "# Hz %c %s R ", parameter_name, format);
    (void)print_value(fp, vfp->vf_dprecision, /*plus=*/false, /*pad=*/false,
	    z0_touchstone);
    (void)fputc('\n', fp);
    if (vfp->vf_type == VNAFILE_TOUCHSTONE2) {
	bool mixed_z0 = false;

	(void)fprintf(fp, "[Number of Ports] %d\n", ports);
	if (ports == 2) {
	    (void)fprintf(fp, "[Two-Port Order] 12_21\n");
	}
	(void)fprintf(fp, "[Number of Frequencies] %d\n",
		vnadata_get_frequencies(vdp));
	for (int i = 1; i < ports; ++i) {
	    if (z0_vector[i] != z0_vector[0]) {
		mixed_z0 = true;
		break;
	    }
	}
	if (mixed_z0) {
	    (void)fprintf(fp, "[Reference]");
	    for (int i = 0; i < ports; ++i) {
		(void)fprintf(fp, " ");
		(void)print_value(fp, vfp->vf_dprecision, /*plus=*/false,
			/*pad=*/false, creal(z0_vector[i]));
	    }
	    (void)fputc('\n', fp);
	}
	(void)fprintf(fp, "[Network Data]\n");
    }
}

/*
 * Function Names
 */
static const char vnafile_check_name[] = "vnafile_check";
static const char vnafile_fsave_name[] = "vnafile_fsave";
static const char vnafile_save_name[] = "vnafile_save";

/*
 * vnafile_save_common: common save routine
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @vdp: input data
 *   @function: function name (for error messages)
 */
static int vnafile_save_common(vnafile_t *vfp, FILE *fp, const char *filename,
	const vnadata_t *vdp, const char *function)
{
    int rows, columns, ports, diagonals, frequencies;
    int aprecision;
    int rc = -1;
    bool auto_type = false;
    const double *frequency_vector;
    const double complex *z0_vector;
    double z0_touchstone;
    vnadata_t *conversions[VPT_NTYPES];

    /*
     * Validate pointer.
     */
    if (vfp == NULL || vfp->vf_magic != VF_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    aprecision = MAX(vfp->vf_dprecision, 3);

    /*
     * Init conversions to NULL.
     */
    (void)memset((void *)conversions, 0, sizeof(conversions));

    /*
     * Validate parameters.  These errors are for the application
     * developer, not the end user, so show the function name.
     */
    if (function == vnafile_fsave_name && fp == NULL) {
	_vnafile_error(vfp, VNAERR_USAGE,
		"%s: error: NULL file pointer", function);
	goto out;
    }
    if (filename == NULL) {
	_vnafile_error(vfp, VNAERR_USAGE,
		"%s: error: NULL filename", function);
	goto out;
    }
    if (vdp == NULL) {
	_vnafile_error(vfp, VNAERR_USAGE,
		"%s: error: NULL vdp", function);
	goto out;
    }

    /*
     * Get the characteristics of the underlying s-parameter matrix.
     */
    rows             = vnadata_get_rows(vdp);
    columns          = vnadata_get_columns(vdp);
    if (vdp->vd_type == VPT_ZIN) {
	if (rows == 1) {
	    rows = columns;
	} else if (columns == 1) {
	    columns = rows;
	}
    }
    ports            = MAX(rows, columns);
    diagonals        = MIN(rows, columns);
    frequencies      = vnadata_get_frequencies(vdp);
    frequency_vector = vnadata_get_frequency_vector(vdp);
    z0_vector        = vnadata_get_z0_vector(vdp);
    z0_touchstone    = z0_vector != NULL ? creal(z0_vector[0]) : 50.0;
    switch (vfp->vf_type) {
    case VNAFILE_AUTO:
	vfp->vf_type = _vnafile_find_type(filename, NULL);
	auto_type = true;
	break;

    case VNAFILE_NATIVE:
    case VNAFILE_TOUCHSTONE1:
    case VNAFILE_TOUCHSTONE2:
	break;

    default:
	_vnafile_error(vfp, VNAERR_USAGE, "%s: error: invalid type (%d)",
		function, vfp->vf_type);
	goto out;
    }

    /*
     * Enforce additional restrictions by file type.
     */
    switch (vfp->vf_type) {
    case VNAFILE_TOUCHSTONE1:
    case VNAFILE_TOUCHSTONE2:
	/*
	 * Touchstone format supports only square matrices with at
	 * least one port.
	 */
	if (rows != columns || ports < 1) {
	    _vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
		    "cannot save %d x %d matrix in touchstone format",
		    function, rows, columns);
	    goto out;
	}

	/*
	 * Touchstone format supports only one parameter type.
	 */
	if (vfp->vf_format_count > 1) {
	    _vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
		    "only a single parameter type may be used in "
		    "touchstone format", function);
	    goto out;
	}

	/*
	 * Touchstone format supports only S, Z, Y, H & G parameters
	 * and only DB, MA and RI formats.
	 */
	{
	    vnafile_format_t *vffp = &vfp->vf_format_vector[0];
	    vnadata_parameter_type_t parameter = vffp->vff_parameter;
	    bool bad = false;

	    if (parameter == VPT_UNDEF) {
		parameter = vdp->vd_type;
	    }
	    switch (parameter) {
	    case VPT_S:
	    case VPT_Z:
	    case VPT_Y:
	    case VPT_H:
	    case VPT_G:
		break;

	    default:
		bad = true;
	    }
	    switch (vffp->vff_format) {
	    case VNAFILE_FORMAT_DB_ANGLE:
	    case VNAFILE_FORMAT_MAG_ANGLE:
	    case VNAFILE_FORMAT_REAL_IMAG:
		break;

	    default:
		bad = true;
	    }
	    if (bad) {
		_vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
			"cannot save parameter %s in touchstone format",
			function, _vnafile_format_to_name(vffp));
		goto out;
	    }
	}

	/*
	 * Touchstone doesn't support frequency-dependent system impedances.
	 */
	if (z0_vector == NULL) {
	    _vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
		    "cannot save frequency-dependent system impedances in "
		    "touchstone format", function);
	    goto out;
	}

	/*
	 * Touchstone doesn't support complex system impedances.
	 */
	for (int i = 0; i < ports; ++i) {
	    if (cimag(z0_vector[i]) != 0.0) {
		_vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
			"cannot save complex system impedances in "
			"touchstone format", function);
		goto out;
	    }
	}

	/*
	 * Finished with Touchstone 2.  If Touchstone 1, apply yet
	 * more constraints.
	 */
	if (vfp->vf_type == VNAFILE_TOUCHSTONE2) {
	    break;
	}

	/*
	 * Touchstone 1 doesn't support more than four ports.
	 */
	if (ports > 4) {
	    if (auto_type) {
		vfp->vf_type = VNAFILE_TOUCHSTONE2;
		break;
	    }
	    _vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
		    "cannot save a system with more than four ports in "
		    "touchstone 1 format", function);
	    goto out;
	}

	/*
	 * Touchstone 1 format doesn't support ports having different
	 * system impedances.
	 */
	for (int i = 1; i < ports; ++i) {
	    if (z0_vector[i] != z0_vector[0]) {
		if (auto_type) {
		    vfp->vf_type = VNAFILE_TOUCHSTONE2;
		    break;
		}
		_vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
			"cannot save ports with different system impedances "
			"in touchstone 1 format", function);
		goto out;
	    }
	}
	break;

    case VNAFILE_NATIVE:
	/*
	 * For native format, disallow nonsensical use of dB for values
	 * that are neither power nor root power.  Unfortunately,
	 * Touchstone allows this, defining dB as 20*log10(cabs(value)),
	 * even for parameters in units of ohms or seimens.
	 */
	for (int i = 0; i < vfp->vf_format_count; ++i) {
	    vnafile_format_t *vffp = &vfp->vf_format_vector[i];
	    vnadata_parameter_type_t parameter = vffp->vff_parameter;

	    if (parameter == VPT_UNDEF) {
		parameter = vdp->vd_type;
	    }
	    if (vffp->vff_format == VNAFILE_FORMAT_DB_ANGLE &&
		    parameter != VPT_S && parameter != VPT_T) {
		_vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
			"%s: in native format, only power or root-power "
			"parameters can be displayed in dB",
			function, _vnafile_format_to_name(vffp));
		goto out;
	    }
	}
	break;

    default:
	abort();
	/*NOTREACHED*/
    }

    /*
     * If the matrix has no off diagonal elements, make sure return loss
     * is not requested.
     */
    if (ports < 2) {
	for (int i = 0; i < vfp->vf_format_count; ++i) {
	    vnafile_format_t *vffp = &vfp->vf_format_vector[i];

	    if (vffp->vff_format == VNAFILE_FORMAT_RL) {
		_vnafile_error(vfp, VNAERR_USAGE, "%s: error: "
			"return loss requires at least one "
			"off-diagonal element", function);
		goto out;
	    }
	}
    }

    /*
     * If touchstone 1, normalize all system impedances to 1
     */
    if (vfp->vf_type == VNAFILE_TOUCHSTONE1 && z0_vector[0] != 1.0) {
	vnadata_t *vdp_copy;
	vnadata_parameter_type_t target_type;

	/*
	 * If the input type is S or T, make a writeable copy.	Otherwise,
	 * convert to S using the existing z0.
	 */
	if ((vdp_copy = vnadata_alloc()) == NULL) {
	    _vnafile_error(vfp, VNAERR_SYSTEM,
		    "vnadata_alloc: %s", strerror(errno));
	    goto out;
	}
	target_type = vnadata_get_type(vdp) == VPT_T ? VPT_T : VPT_S;
	if (vnadata_convert(vdp, vdp_copy, target_type) == -1) {
	    if (errno == EINVAL) {
		_vnafile_error(vfp, VNAERR_USAGE,
			"%s: cannot convert type to %s",
			function, target_type == VPT_T ? "T" : "S");
	    } else {
		;
		_vnafile_error(vfp, VNAERR_SYSTEM,
			"vnadata_convert: %s", strerror(errno));
	    }
	    vnadata_free(vdp_copy);
	    goto out;

	}

	/*
	 * Set all z0's to 1.
	 */
	if (vnadata_set_all_z0(vdp_copy, 1.0) == -1) {
	    _vnafile_error(vfp, VNAERR_SYSTEM,
		    "vnadata_set_all_z0: %s", strerror(errno));
	    vnadata_free(vdp_copy);
	    goto out;
	}

	/*
	 * Save the copy in the conversions array so that it gets
	 * freed at out.  Replace vdp and z0_vector.
	 */
	conversions[target_type] = vdp_copy;
	vdp = vdp_copy;
	z0_vector = vnadata_get_z0_vector(vdp);
    }

    /*
     * Perform all the needed conversions.  We do this up-front so
     * that on error, we won't yet have touched the output file.
     */
    for (int i = 0; i < vfp->vf_format_count; ++i) {
	vnafile_format_t *vffp = &vfp->vf_format_vector[i];

	if (vffp->vff_parameter != VPT_UNDEF) {
	    if (_vnafile_convert(vfp, conversions, vdp, vffp->vff_parameter,
			function) == -1) {
		goto out;
	    }
	}
    }

    /*
     * If vnafile_check, we're done.
     */
    if (function == vnafile_check_name) {
	rc = 0;
	goto out;
    }

    /*
     * Go through the format vector and fix up any instances of "ri",
     * "ma" and "db" without parameter types, taking the parameter type
     * from the vnadata_t structure.
     */
    {
	bool changed = false;

	for (int i = 0; i < vfp->vf_format_count; ++i) {
	    vnafile_format_t *vffp = &vfp->vf_format_vector[i];

	    if (vffp->vff_parameter == VPT_UNDEF) {
		vffp->vff_parameter = vdp->vd_type;
		changed = true;
	    }
	}
	if (changed) {
	    if (_vnafile_update_format_string(vfp) == -1) {
		goto out;
	    }
	}
    }

    /*
     * If vnafile_save, open the output file.
     */
    if (function == vnafile_save_name) {
	if ((fp = fopen(filename, "w")) == NULL) {
	    _vnafile_error(vfp, VNAERR_SYSTEM, "fopen: %s: %s",
		    filename, strerror(errno));
	    goto out;
	}
    }

    /*
     * Print the file header.
     */
    switch (vfp->vf_type) {
    case VNAFILE_NATIVE:
	print_native_header(vfp, fp, vdp);
	break;

    case VNAFILE_TOUCHSTONE1:
    case VNAFILE_TOUCHSTONE2:
	print_touchstone_header(vfp, fp, vdp, z0_touchstone);
	break;

    default:
	abort();
	/*NOTREACHED*/
    }

    /*
     * For each frequency...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	/*
	 * Print the frequency.
	 */
	print_value(fp, vfp->vf_fprecision, /*plus=*/false, /*pad=*/true,
		vnadata_get_frequency(vdp, findex));

	/*
	 * Add frequency-dependent system impedances.
	 */
	if (z0_vector == NULL) {
	    const double complex *fz0_vector;

	    fz0_vector = vnadata_get_fz0_vector(vdp, findex);
	    for (int port = 0; port < ports; ++port) {
		double complex z0 = fz0_vector[port];

		(void)fputc(' ', fp);
		print_value(fp, vfp->vf_dprecision, /*plus=*/true,
			/*pad=*/true, creal(z0));
		(void)fputc(' ', fp);
		print_value(fp, vfp->vf_dprecision, /*plus=*/true,
			/*pad=*/true, cimag(z0));
	    }
	}

	/*
	 * For each parameter...
	 */
	for (int format = 0; format < vfp->vf_format_count; ++format) {
	    const vnafile_format_t *vffp = &vfp->vf_format_vector[format];
	    const vnadata_t *matrix = NULL;
	    const double complex *data = NULL;
	    bool last_arg;
	    bool done = false;

	    /*
	     * Get the needed matrix.
	     */
	    assert(vffp->vff_parameter != VPT_UNDEF);
	    if (vffp->vff_parameter == vdp->vd_type) {
		matrix = vdp;
	    } else {
		matrix = conversions[vffp->vff_parameter];
	    }
	    assert(matrix != NULL);

	    /*
	     * Print
	     */
	    switch (vffp->vff_parameter) {
	    case VPT_S:
		switch (vffp->vff_format) {
		case VNAFILE_FORMAT_IL:
		    for (int row = 0; row < rows; ++row) {
			for (int column = 0; column < columns; ++column) {
			    double complex value;

			    if (row == column) {
				continue;
			    }
			    value = vnadata_get_cell(matrix, findex, row,
				    column);
			    (void)fputc(' ', fp);
			    last_arg = format == vfp->vf_format_count - 1 &&
				row == rows - 1 && column == columns - 1;
			    print_value(fp, vfp->vf_dprecision,
				    /*plus=*/true, /*pad=*/!last_arg,
				    -20.0 * log10(cabs(value)));
			}
		    }
		    done = true;
		    break;

		case VNAFILE_FORMAT_RL:
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
			double complex value;

			value = vnadata_get_cell(matrix, findex, diagonal,
				diagonal);
			(void)fputc(' ', fp);
			last_arg = format == vfp->vf_format_count - 1 &&
			    diagonal == diagonals - 1;
			print_value(fp, vfp->vf_dprecision,
				/*plus=*/true, /*pad=*/!last_arg,
				-20.0 * log10(cabs(value)));
		    }
		    done = true;
		    break;

		case VNAFILE_FORMAT_VSWR:
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
			double complex sxx;
			double a;
			double vswr;

			sxx = vnadata_get_cell(matrix, findex,
				diagonal, diagonal);
			a = cabs(sxx);
			vswr = (1.0 + a) / fabs(1.0 - a);
			(void)fputc(' ', fp);
			last_arg = format == vfp->vf_format_count - 1 &&
			    diagonal == diagonals - 1;
			print_value(fp, vfp->vf_dprecision,
				/*plus=*/false, /*pad=*/!last_arg, vswr);
		    }
		    done = true;
		    break;

		default:
		    break;
		}
		if (done) {
		    break;
		}
		/*FALLTHROUGH*/

	    case VPT_Z:
	    case VPT_Y:
	    case VPT_T:
	    case VPT_H:
	    case VPT_G:
	    case VPT_A:
	    case VPT_B:
		for (int row = 0; row < rows; ++row) {
		    for (int column = 0; column < columns; ++column) {
			double complex value;

			/*
			 * In Touchstone formats, break the line after
			 * every four columns, and, except for two-port,
			 * after every row.
			 */
			if ((vfp->vf_type == VNAFILE_TOUCHSTONE1 ||
			     vfp->vf_type == VNAFILE_TOUCHSTONE2) &&
			     ((column  != 0 && column % 4 == 0) ||
			      (ports != 2 && row != 0 && column == 0))) {
			    (void)fputc('\n', fp);
			    for (int i = 0; i < vfp->vf_fprecision + 5; ++i)
				(void)fputc(' ', fp);
			}

			/*
			 * Format based on vff_format.
			 *
			 * Special case the 2x2 matrix in Touchstone 1
			 * which prints in column major order.
			 */
			if (vfp->vf_type == VNAFILE_TOUCHSTONE1 && ports == 2) {
			    assert(rows == columns);
			    value = vnadata_get_cell(matrix, findex,
				    column, row);
			} else {
			    value = vnadata_get_cell(matrix, findex,
				    row, column);
			}
			switch (vffp->vff_format) {
			case VNAFILE_FORMAT_DB_ANGLE:
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/true, 20.0 * log10(cabs(value)));
			    if (aprecision == VNAFILE_MAX_PRECISION) {
				(void)fprintf(fp, " %+a",
					180.0 / PI * carg(value));
			    } else {
				(void)fprintf(fp, " %+*.*f",
					aprecision + 2,
					aprecision - 3,
					180.0 / PI * carg(value));
			    }
			    break;

			case VNAFILE_FORMAT_MAG_ANGLE:
			    (void)fprintf(fp, "  ");
			    print_value(fp, vfp->vf_dprecision, /*plus=*/false,
				    /*pad=*/true, cabs(value));
			    if (aprecision == VNAFILE_MAX_PRECISION) {
				(void)fprintf(fp, " %+a",
					180.0 / PI * carg(value));
			    } else {
				(void)fprintf(fp, " %+*.*f",
					aprecision + 2,
					aprecision - 3,
					180.0 / PI * carg(value));
			    }
			    break;

			case VNAFILE_FORMAT_REAL_IMAG:
			    last_arg = format == vfp->vf_format_count - 1 &&
				row == rows - 1 && column == columns - 1;
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/true, creal(value));
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, cimag(value));
			    break;

			default:
			    abort();
			    /*NOTREACHED*/
			}
		    }
		}
		break;

	    case VPT_ZIN:
		data = vnadata_get_matrix(matrix, findex);
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    double complex value;

		    value = data[diagonal];
		    switch (vffp->vff_format) {
		    case VNAFILE_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, "  ");
			print_value(fp, vfp->vf_dprecision, /*plus=*/false,
				/*pad=*/true, cabs(value));
			if (aprecision == VNAFILE_MAX_PRECISION) {
			    (void)fprintf(fp, " %+a",
				    180.0 / PI * carg(value));
			} else {
			    (void)fprintf(fp, "  %+*.*f",
				    aprecision + 2,
				    aprecision - 3,
				    180.0 / PI * carg(value));
			}
			break;

		    case VNAFILE_FORMAT_REAL_IMAG:
			last_arg = format == vfp->vf_format_count - 1 &&
			    diagonal == diagonals - 1;
			(void)fputc(' ', fp);
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/true, creal(value));
			(void)fputc(' ', fp);
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, cimag(value));
			break;

		    case VNAFILE_FORMAT_PRC:
			{
			    double complex z;
			    double zr, zi;
			    double r, x, c;

			    z = data[diagonal];
			    zr = creal(z);
			    zi = cimag(z);
			    r = (zr*zr + zi*zi) / zr;
			    x = (zr*zr + zi*zi) / zi;
			    c = -1.0 /
				(2.0 * PI * frequency_vector[findex] * x);
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/true, r);
			    (void)fputc(' ', fp);
			    last_arg = format == vfp->vf_format_count - 1 &&
				diagonal == diagonals - 1;
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, c);
			}
			break;

		    case VNAFILE_FORMAT_PRL:
			{
			    double complex z;
			    double zr, zi;
			    double r, x, l;

			    z = data[diagonal];
			    zr = creal(z);
			    zi = cimag(z);
			    r = (zr*zr + zi*zi) / zr;
			    x = (zr*zr + zi*zi) / zi;
			    l = x / (2.0 * PI * frequency_vector[findex]);
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/true, r);
			    (void)fputc(' ', fp);
			    last_arg = format == vfp->vf_format_count - 1 &&
				diagonal == diagonals - 1;
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, l);
			}
			break;

		    case VNAFILE_FORMAT_SRC:
			{
			    double complex z;
			    double zr, zi;
			    double c;

			    z = data[diagonal];
			    zr = creal(z);
			    zi = cimag(z);
			    c = -1.0 /
				(2.0 * PI * frequency_vector[findex] * zi);
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/true, zr);
			    (void)fputc(' ', fp);
			    last_arg = format == vfp->vf_format_count - 1 &&
				diagonal == diagonals - 1;
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, c);
			}
			break;

		    case VNAFILE_FORMAT_SRL:
			{
			    double complex z;
			    double zr, zi;
			    double l;

			    z = data[diagonal];
			    zr = creal(z);
			    zi = cimag(z);
			    l = zi / (2.0 * PI * frequency_vector[findex]);
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/true, zr);
			    (void)fputc(' ', fp);
			    last_arg = format == vfp->vf_format_count - 1 &&
				diagonal == diagonals - 1;
			    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, l);
			}
			break;

		    default:
			abort();
			/*NOTREACHED*/
		    }
		}
		break;

	    default:
		abort();
		/*NOTREACHED*/
	    }
	}
	(void)fputc('\n', fp);
    }

    /*
     * If Touchstone 2, add the End keyword.
     */
    if (vfp->vf_type == VNAFILE_TOUCHSTONE2) {
	(void)fprintf(fp, "[End]\n");
    }

    /*
     * If vnafile_save, close the output file.
     */
    if (function == vnafile_save_name) {
	if (fclose(fp) == -1) {
	    _vnafile_error(vfp, VNAERR_SYSTEM, "fclose: %s: %s",
		    filename, strerror(errno));
	    fp = NULL;
	    goto out;
	}
	fp = NULL;
    }
    rc = 0;

out:
    if (function == vnafile_save_name && fp != NULL) {
	(void)fclose(fp);
	fp = NULL;
    }
    for (int i = 0; i < VPT_NTYPES; ++i) {
	vnadata_free(conversions[i]);
	conversions[i] = NULL;
    }
    return rc;
}

/*
 * vnafile_check: check that parameters are valid for save
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @filename: file to save
 *   @vdp: input data
 */
int vnafile_check(vnafile_t *vfp, const char *filename,
	const vnadata_t *vdp)
{
    return vnafile_save_common(vfp, NULL, filename, vdp, vnafile_check_name);
}

/*
 * vnafile_fsave: save network parameters to a file pointer
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @vdp: input data
 */
int vnafile_fsave(vnafile_t *vfp, FILE *fp, const char *filename,
	const vnadata_t *vdp)
{
    return vnafile_save_common(vfp, fp, filename, vdp, vnafile_fsave_name);
}

/*
 * vnafile_save: save network parameters to filename
 *   @vfp: pointer to the structure returned from vnafile_alloc
 *   @filename: file to save
 *   @vdp: input data
 */
int vnafile_save(vnafile_t *vfp, const char *filename,
	const vnadata_t *vdp)
{
    return vnafile_save_common(vfp, NULL, filename, vdp, vnafile_save_name);
}
