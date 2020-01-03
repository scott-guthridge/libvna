/*
 * Electrical Network Parameter Conversion Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A11 PARTICULAR PURPOSE.  See the GNU
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

#ifndef MIN
#define MIN(a, b)	((a) <= (b) ? (a) : (b))
#endif /* MIN */

#ifndef MAX
#define MAX(a, b)	((a) >= (b) ? (a) : (b))
#endif /* MAX */

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
 *   @vfp: pointer to the object returned from vnafile_alloc
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
	_vnafile_error(vfp, "%s: %s", function, strerror(errno));
	return -1;
    }
    if (vnadata_convert(vdp, target, type) == -1) {
	if (errno == EINVAL) {
	    _vnafile_error(vfp, "%s: cannot convert from %s to %s",
		function, vnadata_get_typename(vnadata_get_type(vdp)),
		vnadata_get_typename(type));
	    vnadata_free(target);
	    return -1;
	}
	_vnafile_error(vfp, "%s: vnadata_convert: %s",
		function, strerror(errno));
	vnadata_free(target);
	return -1;
    }
    conversions[type] = target;
    return 0;
}

/*
 * _vnafile_get_matrix: return the matrix of the requested type
 *   @conversions: cached conversions
 *   @vdp: input data
 *   @type: desired type
 */
static const vnadata_t *_vnafile_get_matrix(vnadata_t **conversions,
	const vnadata_t *vdp, vnadata_parameter_type_t type)
{
    if (conversions[type] != NULL) {
	return conversions[type];
    }
    if (type == vnadata_get_type(vdp)) {
	return vdp;
    }
    return NULL;
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
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @fp: file pointer
 *   @vdp: input data
 */
static void print_native_header(vnafile_t *vfp, FILE *fp, const vnadata_t *vdp)
{
    int rows, columns, ports, diagonals;
    int output_fields = 1;
    int current_field = 0;
    int field_width;
    int port_width;
    char port_buf[2 * 3 * sizeof(int) + 1];
    const double complex *z0_vector = vnadata_get_z0_vector(vdp);

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

    /*
     * If vdp is only a vector, fix rows and columns to reflect those
     * of the underlying s-parameter matrix.
     */
    if (vdp->vd_type == VPT_ZIN) {
	if (rows == 1) {
	    rows = columns;
	} else if (columns == 1) {
	    columns = rows;
	}
    }

    /*
     * Determine the number of digits needed to display the port number.
     */
    if ((ports + 1) < 10) {
	port_width = 1;
    } else {
	port_width = (int)floor(log10(ports + 0.5)) + 1;
    }

    /*
     * Determine the number of output fields and from it determine
     * the number of digits needed to display the field number.
     */
    if (z0_vector == NULL) {
	output_fields += 2 * ports;
    }
    for (int format = 0; format < vfp->vf_format_count; ++format) {
	const vnafile_format_t *vffp = &vfp->vf_format_vector[format];

	switch (vffp->vff_parameter) {
	case VNAFILE_PARAMETER_S:
	case VNAFILE_PARAMETER_Z:
	case VNAFILE_PARAMETER_Y:
	    output_fields += 2 * rows * columns;
	    break;

	case VNAFILE_PARAMETER_T:
	case VNAFILE_PARAMETER_H:
	case VNAFILE_PARAMETER_G:
	case VNAFILE_PARAMETER_A:
	case VNAFILE_PARAMETER_B:
	    output_fields += 8;
	    break;

	case VNAFILE_PARAMETER_ZIN:
	    output_fields += 2 * ports;
	    break;

	case VNAFILE_PARAMETER_IL:
	    output_fields += rows * columns - diagonals;
	    break;

	case VNAFILE_PARAMETER_RL:
	    output_fields += diagonals;
	    break;

	case VNAFILE_PARAMETER_PRC:
	case VNAFILE_PARAMETER_PRL:
	case VNAFILE_PARAMETER_SRC:
	case VNAFILE_PARAMETER_SRL:
	    output_fields += 2 * diagonals;
	    break;

	case VNAFILE_PARAMETER_VSWR:
	    output_fields += diagonals;
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
	    //(void)fputc('j', fp); ZZ
	}
	(void)fputc('\n', fp);
    }
    (void)fprintf(fp, "#:fprecision %d\n", vfp->vf_fprecision);
    (void)fprintf(fp, "#:dprecision %d\n", vfp->vf_dprecision);
    (void)fprintf(fp, "#\n");

    /*
     * Print a key for each field.
     */
    (void)fprintf(fp, "# field %*d: frequency     (Hz)\n",
	    field_width, ++current_field);
    if (z0_vector == NULL) {
	int pwidth;

	if (ports > 9) {
	    pwidth = 2 * port_width + 1;
	} else {
	    pwidth = 2 * port_width;
	}
	for (int port = 0; port < ports; ++port) {
	    (void)sprintf(port_buf, "%d", port + 1);
	    (void)fprintf(fp, "# field %*d: Z%-*s real      (ohms)\n",
		    field_width, ++current_field, pwidth, port_buf);
	    (void)fprintf(fp, "# field %*d: Z%-*s imaginary (ohms)\n",
		    field_width, ++current_field, pwidth, port_buf);
	}
    }
    for (int format = 0; format < vfp->vf_format_count; ++format) {
	const vnafile_format_t *vffp = &vfp->vf_format_vector[format];
	char name;
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

	switch (vffp->vff_parameter) {
	    do {
	case VNAFILE_PARAMETER_S:
		name = 'S';
		type_vector = st_types;
		break;

	case VNAFILE_PARAMETER_Z:
		name = 'Z';
		type_vector = z_types;
		break;

	case VNAFILE_PARAMETER_Y:
		name = 'Y';
		type_vector = y_types;
		break;

	case VNAFILE_PARAMETER_T:
		name = 'T';
		type_vector = st_types;
		break;

	case VNAFILE_PARAMETER_H:
		name = 'H';
		type_vector = h_types;
		break;

	case VNAFILE_PARAMETER_G:
		name = 'G';
		type_vector = g_types;
		break;

	case VNAFILE_PARAMETER_A:
		name = 'A';
		type_vector = ab_types;
		break;

	case VNAFILE_PARAMETER_B:
		name = 'B';
		type_vector = ab_types;
		break;
	    } while (false);
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    const char *type = (type_vector[1] == NULL) ?
			type_vector[0] : type_vector[row * columns + column];
		    int pwidth;

		    if (ports > 9) {
			(void)sprintf(port_buf, "%d,%d", row + 1, column + 1);
			pwidth = 2 * port_width + 1;
		    } else {
			(void)sprintf(port_buf, "%d%d", row + 1, column + 1);
			pwidth = 2 * port_width;
		    }

		    (void)fprintf(fp, "# field %*d: %c%-*s",
			field_width, ++current_field, name, pwidth,
			port_buf);
		    switch (vffp->vff_coordinates) {
		    case VNAFILE_COORDINATES_REAL_IMAG:
			(void)fprintf(fp, " real      (%s)\n", type);
			break;
		    case VNAFILE_COORDINATES_MAG_ANGLE:
			(void)fprintf(fp, " magnitude (%s)\n", type);
			break;
		    case VNAFILE_COORDINATES_DB_ANGLE:
			(void)fprintf(fp, " magnitude (dB)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		    (void)fprintf(fp, "# field %*d: %c%-*s",
			field_width, ++current_field, name, pwidth,
			port_buf);
		    switch (vffp->vff_coordinates) {
		    case VNAFILE_COORDINATES_REAL_IMAG:
			(void)fprintf(fp, " imaginary (%s)\n", type);
			break;
		    case VNAFILE_COORDINATES_MAG_ANGLE:
		    case VNAFILE_COORDINATES_DB_ANGLE:
			(void)fprintf(fp, " angle     (degrees)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		}
	    }
	    break;

	case VNAFILE_PARAMETER_ZIN:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: Zin%-*s",
                    field_width, ++current_field, port_width, port_buf);
		switch (vffp->vff_coordinates) {
		case VNAFILE_COORDINATES_REAL_IMAG:
		    (void)fprintf(fp, " real      (ohms)\n");
		    break;
		case VNAFILE_COORDINATES_MAG_ANGLE:
		    (void)fprintf(fp, " magnitude (ohms)\n");
		    break;
		default:
		    abort();
		    /*NOTREACHED*/
		}
                (void)fprintf(fp, "# field %*d: Zin%-*s",
                    field_width, ++current_field, port_width, port_buf);
		switch (vffp->vff_coordinates) {
		case VNAFILE_COORDINATES_REAL_IMAG:
		    (void)fprintf(fp, " imaginary (ohms)\n");
		    break;
		case VNAFILE_COORDINATES_MAG_ANGLE:
		    (void)fprintf(fp, " angle     (degrees)\n");
		    break;
		default:
		    abort();
		    /*NOTREACHED*/
		}
            }
            break;

	case VNAFILE_PARAMETER_IL:
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < columns; ++column) {
		    int pwidth;

		    if (row == column) {
			continue;
		    }
		    if (ports > 9) {
			(void)sprintf(port_buf, "%d,%d", row + 1, column + 1);
			pwidth = 2 * port_width + 1;
		    } else {
			(void)sprintf(port_buf, "%d%d", row + 1, column + 1);
			pwidth = 2 * port_width;
		    }
		    (void)fprintf(fp, "# field %*d: IL%-*s magnitude (dB)\n",
			field_width, ++current_field, pwidth, port_buf);
		}
	    }
            break;

	case VNAFILE_PARAMETER_RL:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: RL%-*s magnitude (dB)\n",
                    field_width, ++current_field, port_width, port_buf);
            }
            break;

	case VNAFILE_PARAMETER_PRC:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: PRC%-*s R         (ohms)\n",
                    field_width, ++current_field, port_width, port_buf);
                (void)fprintf(fp, "# field %*d: PRC%-*s C         (farads)\n",
                    field_width, ++current_field, port_width, port_buf);
            }
            break;

	case VNAFILE_PARAMETER_PRL:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: PRL%-*s R         (ohms)\n",
                    field_width, ++current_field, port_width, port_buf);
                (void)fprintf(fp, "# field %*d: PRL%-*s L         (henries)\n",
                    field_width, ++current_field, port_width, port_buf);
            }
            break;

	case VNAFILE_PARAMETER_SRC:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: SRC%-*s R         (ohms)\n",
                    field_width, ++current_field, port_width, port_buf);
                (void)fprintf(fp, "# field %*d: SRC%-*s C         (farads)\n",
                    field_width, ++current_field, port_width, port_buf);
            }
            break;

	case VNAFILE_PARAMETER_SRL:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: SRL%-*s R         (ohms)\n",
                    field_width, ++current_field, port_width, port_buf);
                (void)fprintf(fp, "# field %*d: SRL%-*s L         (henries)\n",
                    field_width, ++current_field, port_width, port_buf);
            }
            break;

	case VNAFILE_PARAMETER_VSWR:
            for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
                (void)sprintf(port_buf, "%d", diagonal + 1);
                (void)fprintf(fp, "# field %*d: VSWR%-*s\n",
                    field_width, ++current_field, port_width, port_buf);
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
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @fp: file pointer
 *   @vdp: input data
 *   @z0_touchstone: the first system impedance (before normalization)
 */
static void print_touchstone_header(vnafile_t *vfp, FILE *fp,
	const vnadata_t *vdp, double z0_touchstone)
{
    vnafile_format_t *vffp = &vfp->vf_format_vector[0];
    char parameter;
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
    case VNAFILE_PARAMETER_S:
	parameter = 'S';
	break;
    case VNAFILE_PARAMETER_Z:
	parameter = 'Z';
	break;
    case VNAFILE_PARAMETER_Y:
	parameter = 'Y';
	break;
    case VNAFILE_PARAMETER_H:
	parameter = 'H';
	break;
    case VNAFILE_PARAMETER_G:
	parameter = 'G';
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    switch (vffp->vff_coordinates) {
    case VNAFILE_COORDINATES_DB_ANGLE:
	format = "DB";
	break;
    case VNAFILE_COORDINATES_MAG_ANGLE:
	format = "MA";
	break;
    case VNAFILE_COORDINATES_REAL_IMAG:
	format = "RI";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)fprintf(fp, "# Hz %c %s R ", parameter, format);
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
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @vdp: input data
 *   @function: function name (for error messages)
 */
static int vnafile_save_common(vnafile_t *vfp, FILE *fp, const char *filename,
	const vnadata_t *vdp, const char *function)
{
    int rows, columns, ports, diagonals, frequencies;
    int aprecision = MAX(vfp->vf_dprecision, 3);
    int rc = -1;
    bool auto_type = false;
    const double *frequency_vector;
    const double complex *z0_vector;
    double z0_touchstone;
    vnadata_t *conversions[VPT_NTYPES];

    /*
     * Init conversions to NULL.
     */
    (void)memset((void *)conversions, 0, sizeof(conversions));

    /*
     * Validate parameters.  These errors are for the application
     * developer, not the end user, so show the function name.
     */
    if (vfp == NULL) {
	errno = EINVAL;
	goto out;
    }
    if (function == vnafile_fsave_name && fp == NULL) {
	_vnafile_error(vfp, "%s: error: NULL file pointer", function);
	errno = EINVAL;
	goto out;
    }
    if (filename == NULL) {
	_vnafile_error(vfp, "%s: error: NULL filename", function);
	errno = EINVAL;
	goto out;
    }
    if (vdp == NULL) {
	_vnafile_error(vfp, "%s: error: NULL vdp", function);
	errno = EINVAL;
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
	_vnafile_error(vfp, "%s: error: invalid type (%d)",
		function, vfp->vf_type);
	errno = EINVAL;
	goto out;
    }

    /*
     * Enforce additional restrictions by file type.
     */
    switch (vfp->vf_type) {
    case VNAFILE_TOUCHSTONE1:
    case VNAFILE_TOUCHSTONE2:
	/*
	 * Touchstone format supports only square matrices.
	 */
	if (rows != columns) {
	    _vnafile_error(vfp, "%s: error: cannot save %d x %d matrix in "
		    "touchstone format", function, rows, columns);
	    errno = EINVAL;
	    goto out;
	}

	/*
	 * Touchstone format supports only one parameter type.
	 */
	if (vfp->vf_format_count > 1) {
	    _vnafile_error(vfp, "%s: error: only a single parameter type "
		    "may be used in touchstone format", function);
	    errno = EINVAL;
	    goto out;
	}

	/*
	 * Touchstone format supports only S, Z, Y, H & G parameters.
	 */
	{
	    vnafile_format_t *vffp = &vfp->vf_format_vector[0];

	    switch (vffp->vff_parameter) {
	    case VNAFILE_PARAMETER_S:
	    case VNAFILE_PARAMETER_Z:
	    case VNAFILE_PARAMETER_Y:
	    case VNAFILE_PARAMETER_H:
	    case VNAFILE_PARAMETER_G:
		break;

	    default:
		_vnafile_error(vfp, "%s: error: cannot save parameter %s in "
			"touchstone format", function,
			_vnafile_format_to_name(vffp));
		errno = EINVAL;
		goto out;
	    }
	}

	/*
	 * Touchstone doesn't support frequency-dependent system impedances.
	 */
	if (z0_vector == NULL) {
	    _vnafile_error(vfp, "%s: error: cannot save frequency-dependent "
		    "system mpedances in touchstone format", function);
	    errno = EINVAL;
	    goto out;
	}

	/*
	 * Touchstone doesn't support complex system impedances.
	 */
	for (int i = 0; i < ports; ++i) {
	    if (cimag(z0_vector[i]) != 0.0) {
		_vnafile_error(vfp, "%s: error: cannot save complex system "
			"impedances in touchstone format", function);
		errno = EINVAL;
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
	    _vnafile_error(vfp, "%s: error: cannot save a system "
		    "with more than four ports in touchstone 1 "
		    "format", function);
	    errno = EINVAL;
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
		_vnafile_error(vfp, "%s: error: cannot save ports with "
			"different system impedances in touchstone 1 "
			"format", function);
		errno = EINVAL;
		goto out;
	    }
	}
	break;

    case VNAFILE_NATIVE:
	/*
	 * In native format, disallow nonsensical use of dB for values
	 * that are neither power nor root power.  Unfortunately,
	 * Touchstone allows this, defining dB as 20*log10(cabs(value)),
	 * even for parameters that are dimensionless or in units of
	 * ohms or seimens.
	 */
	for (int i = 0; i < vfp->vf_format_count; ++i) {
	    vnafile_format_t *vffp = &vfp->vf_format_vector[i];

	    switch (vffp->vff_parameter) {
	    case VNAFILE_PARAMETER_Z:
	    case VNAFILE_PARAMETER_Y:
	    case VNAFILE_PARAMETER_H:
	    case VNAFILE_PARAMETER_G:
	    case VNAFILE_PARAMETER_A:
	    case VNAFILE_PARAMETER_B:
	    case VNAFILE_PARAMETER_ZIN:
		if (vffp->vff_coordinates == VNAFILE_COORDINATES_DB_ANGLE) {
		    _vnafile_error(vfp, "%s: error: %s: in native format, only "
			    "power or root-power parameters can be displayed "
			    "in dB", function, _vnafile_format_to_name(vffp));
		    errno = EINVAL;
		    goto out;
		}
		break;

	    default:
		break;
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

	    if (vffp->vff_parameter == VNAFILE_PARAMETER_RL) {
		_vnafile_error(vfp, "%s: error: return loss requires at "
			"least one off-diagonal element", function);
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
	    _vnafile_error(vfp, "%s: %s", function, strerror(errno));
	    goto out;
	}
	target_type = vnadata_get_type(vdp) == VPT_T ? VPT_T : VPT_S;
	if (vnadata_convert(vdp, vdp_copy, target_type) == -1) {
	    _vnafile_error(vfp, "%s: vnadata_convert: %s",
		    function, strerror(errno));
	    vnadata_free(vdp_copy);
	    goto out;

	}

	/*
	 * Set all z0's to 1.
	 */
	if (vnadata_set_all_z0(vdp_copy, 1.0) == -1) {
	    _vnafile_error(vfp, "%s: vnadata_set_all_z0: %s",
		    function, strerror(errno));
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
	vnadata_parameter_type_t target_type;

	switch (vffp->vff_parameter) {
	case VNAFILE_PARAMETER_S:
	case VNAFILE_PARAMETER_IL:
	case VNAFILE_PARAMETER_RL:
	case VNAFILE_PARAMETER_VSWR:
	    target_type = VPT_S;
	    break;

	case VNAFILE_PARAMETER_Z:
	    target_type = VPT_Z;
	    break;

	case VNAFILE_PARAMETER_Y:
	    target_type = VPT_Y;
	    break;

	case VNAFILE_PARAMETER_T:
	    target_type = VPT_T;
	    break;

	case VNAFILE_PARAMETER_H:
	    target_type = VPT_H;
	    break;

	case VNAFILE_PARAMETER_G:
	    target_type = VPT_G;
	    break;

	case VNAFILE_PARAMETER_A:
	    target_type = VPT_A;
	    break;

	case VNAFILE_PARAMETER_B:
	    target_type = VPT_B;
	    break;

	case VNAFILE_PARAMETER_ZIN:
	case VNAFILE_PARAMETER_PRC:
	case VNAFILE_PARAMETER_PRL:
	case VNAFILE_PARAMETER_SRC:
	case VNAFILE_PARAMETER_SRL:
	    target_type = VPT_ZIN;
	    break;

	default:
	    abort();
	    /*NOTREACHED*/
	}
	if (_vnafile_convert(vfp, conversions, vdp, target_type,
			     function) == -1) {
	    goto out;
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
     * If vnafile_save, open the output file.
     */
    if (function == vnafile_save_name) {
	if ((fp = fopen(filename, "w")) == NULL) {
	    _vnafile_error(vfp, "fopen: %s: %s", filename, strerror(errno));
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
	    bool last_arg;
	    
	    /*
	     * Get the needed matrix.
	     */
	    switch (vffp->vff_parameter) {
	    case VNAFILE_PARAMETER_S:
	    case VNAFILE_PARAMETER_IL:
	    case VNAFILE_PARAMETER_RL:
	    case VNAFILE_PARAMETER_VSWR:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_S);
		break;

	    case VNAFILE_PARAMETER_Z:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_Z);
		break;

	    case VNAFILE_PARAMETER_Y:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_Y);
		break;

	    case VNAFILE_PARAMETER_T:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_T);
		break;

	    case VNAFILE_PARAMETER_H:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_H);
		break;

	    case VNAFILE_PARAMETER_G:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_G);
		break;

	    case VNAFILE_PARAMETER_A:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_A);
		break;

	    case VNAFILE_PARAMETER_B:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_B);
		break;

	    case VNAFILE_PARAMETER_ZIN:
	    case VNAFILE_PARAMETER_PRC:
	    case VNAFILE_PARAMETER_PRL:
	    case VNAFILE_PARAMETER_SRC:
	    case VNAFILE_PARAMETER_SRL:
		matrix = _vnafile_get_matrix(conversions, vdp, VPT_ZIN);
		break;

	    default:
		abort();
		/*NOTREACHED*/
	    }
	    assert(matrix != NULL);

	    /*
	     * Print
	     */
	    switch (vffp->vff_parameter) {
	    case VNAFILE_PARAMETER_S:
	    case VNAFILE_PARAMETER_Z:
	    case VNAFILE_PARAMETER_Y:
	    case VNAFILE_PARAMETER_T:
	    case VNAFILE_PARAMETER_H:
	    case VNAFILE_PARAMETER_G:
	    case VNAFILE_PARAMETER_A:
	    case VNAFILE_PARAMETER_B:
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
			 * Format according to given coordinates.
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
			switch (vffp->vff_coordinates) {
			case VNAFILE_COORDINATES_DB_ANGLE:
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

			case VNAFILE_COORDINATES_MAG_ANGLE:
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

			case VNAFILE_COORDINATES_REAL_IMAG:
			    last_arg = format == vfp->vf_format_count - 1 &&
				row == rows - 1 && column == columns - 1;
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_fprecision, /*plus=*/true,
				    /*pad=*/true, creal(value));
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_fprecision, /*plus=*/true,
				    /*pad=*/!last_arg, cimag(value));
			    break;

			default:
			    abort();
			    /*NOTREACHED*/
			}
		    }
		}
		break;

	    case VNAFILE_PARAMETER_ZIN:
		{
		    const double complex *data;

		    data = vnadata_get_matrix(matrix, findex);
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
			double complex value;

			value = data[diagonal];
			switch (vffp->vff_coordinates) {
			case VNAFILE_COORDINATES_MAG_ANGLE:
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

			case VNAFILE_COORDINATES_REAL_IMAG:
			    last_arg = format == vfp->vf_format_count - 1;
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_fprecision, /*plus=*/true,
				    /*pad=*/true, creal(value));
			    (void)fputc(' ', fp);
			    print_value(fp, vfp->vf_fprecision, /*plus=*/true,
				    /*pad=*/!last_arg, cimag(value));
			    break;

			default:
			    abort();
			    /*NOTREACHED*/
			}
		    }
		}
		break;

	    case VNAFILE_PARAMETER_IL:
		for (int row = 0; row < rows; ++row) {
		    for (int column = 0; column < columns; ++column) {
			double complex value;

			if (row == column) {
			    continue;
			}
			value = vnadata_get_cell(matrix, findex, row, column);
			(void)fputc(' ', fp);
			last_arg = format == vfp->vf_format_count - 1;
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, -20.0 * log10(cabs(value)));
		    }
		}
		break;

	    case VNAFILE_PARAMETER_RL:
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    double complex value;

		    value = vnadata_get_cell(matrix, findex, diagonal,
			    diagonal);
		    (void)fputc(' ', fp);
		    last_arg = format == vfp->vf_format_count - 1;
		    print_value(fp, vfp->vf_dprecision, /*plus=*/true,
			    /*pad=*/!last_arg, -20.0 * log10(cabs(value)));
		}
		break;

	    case VNAFILE_PARAMETER_PRC:
		{
		    const double complex *data;

		    data = vnadata_get_matrix(matrix, findex);
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
			double complex z;
			double zr, zi;
			double r, x, c;

			z = data[diagonal];
			zr = creal(z);
			zi = cimag(z);
			r = (zr*zr + zi*zi) / zr;
			x = (zr*zr + zi*zi) / zi;
			c = -1.0 / (2.0 * PI * frequency_vector[findex] * x);
			(void)fputc(' ', fp);
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/true, r);
			(void)fputc(' ', fp);
			last_arg = format == vfp->vf_format_count - 1;
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, c);
		    }
		}
		break;

	    case VNAFILE_PARAMETER_PRL:
		{
		    const double complex *data;

		    data = vnadata_get_matrix(matrix, findex);
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
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
			last_arg = format == vfp->vf_format_count - 1;
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, l);
		    }
		}
		break;

	    case VNAFILE_PARAMETER_SRC:
		{
		    const double complex *data;

		    data = vnadata_get_matrix(matrix, findex);
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
			double complex z;
			double zr, zi;
			double c;

			z = data[diagonal];
			zr = creal(z);
			zi = cimag(z);
			c = -1.0 / (2.0 * PI * frequency_vector[findex] * zi);
			(void)fputc(' ', fp);
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/true, zr);
			(void)fputc(' ', fp);
			last_arg = format == vfp->vf_format_count - 1;
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, c);
		    }
		}
		break;

	    case VNAFILE_PARAMETER_SRL:
		{
		    const double complex *data;

		    data = vnadata_get_matrix(matrix, findex);
		    for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
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
			last_arg = format == vfp->vf_format_count - 1;
			print_value(fp, vfp->vf_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, l);
		    }
		}
		break;

	    case VNAFILE_PARAMETER_VSWR:
		for (int diagonal = 0; diagonal < diagonals; ++diagonal) {
		    double complex sxx;
		    double a;
		    double vswr;

		    sxx = vnadata_get_cell(matrix, findex, diagonal, diagonal);
		    a = cabs(sxx);
		    vswr = (1.0 + a) / fabs(1.0 - a);
		    (void)fputc(' ', fp);
		    last_arg = format == vfp->vf_format_count - 1;
		    print_value(fp, vfp->vf_dprecision, /*plus=*/false,
			    /*pad=*/!last_arg, vswr);
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
	    _vnafile_error(vfp, "fclose: %s: %s", filename, strerror(errno));
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
 *   @vfp: pointer to the object returned from vnafile_alloc
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
 *   @vfp: pointer to the object returned from vnafile_alloc
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
 *   @vfp: pointer to the object returned from vnafile_alloc
 *   @filename: file to save
 *   @vdp: input data
 */
int vnafile_save(vnafile_t *vfp, const char *filename,
	const vnadata_t *vdp)
{
    return vnafile_save_common(vfp, NULL, filename, vdp, vnafile_save_name);
}
