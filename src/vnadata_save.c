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
     * VNADATA_MAX_PRECISION, write hexadecimal floating point format.
     */
    if (precision < 1) {
	precision = 1;

    } else if (precision == VNADATA_MAX_PRECISION) {
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
 * convert_input: convert the input matrix to the given type
 *   @function: function name (for error messages)
 *   @vdip:   internal parameter matrix
 *   @conversions: cached conversions
 *   @type: desired type
 */
static int convert_input(const char *function, vnadata_internal_t *vdip,
	vnadata_t **conversions, vnadata_parameter_type_t type)
{
    vnadata_t *vdp = &vdip->vdi_vd;
    vnadata_t *target = NULL;

    if (conversions[type] != NULL) {
	return 0;
    }
    if (type == vnadata_get_type(vdp)) {
	return 0;
    }
    if ((target = vnadata_alloc(vdip->vdi_error_fn,
		    vdip->vdi_error_arg)) == NULL) {
	return -1;
    }
    if (vnadata_convert(vdp, target, type) == -1) {
	vnadata_free(target);
	return -1;
    }
    conversions[type] = target;
    return 0;
}

/*
 * print_npd_header: print header for NPD format
 *   @vdip:   internal parameter matrix
 *   @fp: file pointer
 *   @vdp: input data
 */
static void print_npd_header(vnadata_internal_t *vdip, FILE *fp)
{
    vnadata_t *vdp = &vdip->vdi_vd;
    int rows, ports;
    int output_fields = 1;
    int current_field = 0;
    int field_width;
    int port_width, port_pair_width;
    int parameter_width = 0;
    char parameter_buf[3 + 2 * 3 * sizeof(int) + 1];
    const double complex *z0_vector = NULL;

    /*
     * Set z0_vector if we're not using per-frequency z0 values.
     */
    if (!(vdip->vdi_flags & VF_PER_F_Z0)) {
	z0_vector = vdip->vdi_z0_vector;
    }

    /*
     * Get the number of rows and columns.  Columns is the number
     * of ports.
     */
    rows  = vnadata_get_rows(vdp);
    ports = vnadata_get_columns(vdp);
    if (vdp->vd_type == VPT_ZIN) {
	assert(rows == 1);
    } else {
	assert(rows == ports);
    }

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
    for (int format = 0; format < vdip->vdi_format_count; ++format) {
	const vnadata_format_descriptor_t *vfdp =
	    &vdip->vdi_format_vector[format];

	switch (vfdp->vfd_format) {
	case VNADATA_FORMAT_DB_ANGLE:
	case VNADATA_FORMAT_MAG_ANGLE:
	case VNADATA_FORMAT_REAL_IMAG:
	    if (vfdp->vfd_parameter != VPT_ZIN) {
		output_fields += 2 * rows * ports;
		parameter_width = MAX(parameter_width, 1 + port_pair_width);
	    } else {
		output_fields += ports;
		parameter_width = MAX(parameter_width, 3 + port_width);
	    }
	    break;

	case VNADATA_FORMAT_PRC:
	case VNADATA_FORMAT_PRL:
	case VNADATA_FORMAT_SRC:
	case VNADATA_FORMAT_SRL:
	    output_fields += 2 * ports;
	    parameter_width = MAX(parameter_width, 3 + port_width);
	    break;

	case VNADATA_FORMAT_IL:
	    output_fields += ports * (ports - 1);
	    parameter_width = MAX(parameter_width, 2 + port_pair_width);
	    break;

	case VNADATA_FORMAT_RL:
	    output_fields += ports;
	    parameter_width = MAX(parameter_width, 2 + port_width);
	    break;

	case VNADATA_FORMAT_VSWR:
	    output_fields += ports;
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
    (void)fprintf(fp, "#NPD\n");
    (void)fprintf(fp, "#:version 1.0\n");
    (void)fprintf(fp, "#:ports %d\n", ports);
    (void)fprintf(fp, "#:frequencies %d\n", vnadata_get_frequencies(vdp));
    (void)fprintf(fp, "#:parameters %s\n", vdip->vdi_format_string);
    (void)fprintf(fp, "#:z0");
    if (z0_vector == NULL) {
	(void)fprintf(fp, " PER-FREQUENCY\n");
    } else {
	for (int port = 0; port < ports; ++port) {
	    double complex z0 = z0_vector[port];

	    (void)fputc(' ', fp);
	    print_value(fp, vdip->vdi_dprecision, /*plus=*/false, /*pad=*/false,
		    creal(z0));
	    (void)fputc(' ', fp);
	    print_value(fp, vdip->vdi_dprecision, /*plus=*/true, /*pad=*/false,
		    cimag(z0));
	    (void)fputc('j', fp);
	}
	(void)fputc('\n', fp);
    }
    (void)fprintf(fp, "#:fprecision %d\n", vdip->vdi_fprecision);
    (void)fprintf(fp, "#:dprecision %d\n", vdip->vdi_dprecision);
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
    for (int format = 0; format < vdip->vdi_format_count; ++format) {
	const vnadata_format_descriptor_t *vfdp =
	    &vdip->vdi_format_vector[format];
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
	switch (vfdp->vfd_format) {
	case VNADATA_FORMAT_DB_ANGLE:
	case VNADATA_FORMAT_MAG_ANGLE:
	case VNADATA_FORMAT_REAL_IMAG:
	    /*
	     * Handle Zin.
	     */
	    if (vfdp->vfd_parameter == VPT_ZIN) {
		for (int port = 0; port < ports; ++port) {
		    (void)sprintf(parameter_buf, "Zin%d", port + 1);
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vfdp->vfd_format) {
		    case VNADATA_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " real      (ohms)\n");
			break;
		    case VNADATA_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, " magnitude (ohms)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vfdp->vfd_format) {
		    case VNADATA_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " imaginary (ohms)\n");
			break;
		    case VNADATA_FORMAT_MAG_ANGLE:
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
	    switch (vfdp->vfd_parameter) {
	    case VPT_S:
		name = "S";
		type_vector = st_types;
		break;

	    case VPT_T:
		name = "T";
		type_vector = st_types;
		break;

	    case VPT_U:
		name = "U";
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
		for (int column = 0; column < ports; ++column) {
		    const char *type = (type_vector[1] == NULL) ?
			type_vector[0] : type_vector[row * ports + column];

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
		    switch (vfdp->vfd_format) {
		    case VNADATA_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " real      (%s)\n", type);
			break;
		    case VNADATA_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, " magnitude (%s)\n", type);
			break;
		    case VNADATA_FORMAT_DB_ANGLE:
			(void)fprintf(fp, " magnitude (dB)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		    (void)fprintf(fp, "# field %*d: %-*s",
			field_width, ++current_field,
			parameter_width, parameter_buf);
		    switch (vfdp->vfd_format) {
		    case VNADATA_FORMAT_REAL_IMAG:
			(void)fprintf(fp, " imaginary (%s)\n", type);
			break;
		    case VNADATA_FORMAT_MAG_ANGLE:
		    case VNADATA_FORMAT_DB_ANGLE:
			(void)fprintf(fp, " angle     (degrees)\n");
			break;
		    default:
			abort();
			/*NOTREACHED*/
		    }
		}
	    }
	    break;

	case VNADATA_FORMAT_PRC:
	    assert(vfdp->vfd_parameter == VPT_ZIN);
	    for (int port = 0; port < ports; ++port) {
		(void)sprintf(parameter_buf, "PRC%d", port + 1);
		(void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
		(void)fprintf(fp, "# field %*d: %-*s C         (farads)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	    }
	    break;

	case VNADATA_FORMAT_PRL:
	    assert(vfdp->vfd_parameter == VPT_ZIN);
	    for (int port = 0; port < ports; ++port) {
		(void)sprintf(parameter_buf, "PRL%d", port + 1);
		(void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
		(void)fprintf(fp, "# field %*d: %-*s L         (henries)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	    }
	    break;

	case VNADATA_FORMAT_SRC:
	    assert(vfdp->vfd_parameter == VPT_ZIN);
	    for (int port = 0; port < ports; ++port) {
		(void)sprintf(parameter_buf, "SRC%d", port + 1);
		(void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
		(void)fprintf(fp, "# field %*d: %-*s C         (farads)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	    }
	    break;

	case VNADATA_FORMAT_SRL:
	    assert(vfdp->vfd_parameter == VPT_ZIN);
	    for (int port = 0; port < ports; ++port) {
		(void)sprintf(parameter_buf, "SRL%d", port + 1);
		(void)fprintf(fp, "# field %*d: %-*s R         (ohms)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
		(void)fprintf(fp, "# field %*d: %-*s L         (henries)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	    }
	    break;

	case VNADATA_FORMAT_IL:
	    assert(vfdp->vfd_parameter == VPT_S);
	    for (int row = 0; row < rows; ++row) {
		for (int column = 0; column < ports; ++column) {
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

	case VNADATA_FORMAT_RL:
	    assert(vfdp->vfd_parameter == VPT_S);
	    for (int port = 0; port < ports; ++port) {
		(void)sprintf(parameter_buf, "RL%d", port + 1);
		(void)fprintf(fp, "# field %*d: %-*s magnitude (dB)\n",
		    field_width, ++current_field, parameter_width,
		    parameter_buf);
	    }
	    break;

	case VNADATA_FORMAT_VSWR:
	    assert(vfdp->vfd_parameter == VPT_S);
	    for (int port = 0; port < ports; ++port) {
		(void)sprintf(parameter_buf, "VSWR%d", port + 1);
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
 *   @vdip:   internal parameter matrix
 *   @fp: file pointer
 *   @vdp: input data
 *   @z0_touchstone: the first system impedance (before normalization)
 */
static void print_touchstone_header(vnadata_internal_t *vdip, FILE *fp,
	double z0_touchstone)
{
    vnadata_t *vdp = &vdip->vdi_vd;
    vnadata_format_descriptor_t *vfdp = &vdip->vdi_format_vector[0];
    char parameter_name;
    const char *format;
    int ports = vnadata_get_rows(vdp);
    const double complex *z0_vector;

    assert(vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE1 ||
	   vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE2);
    assert(vdip->vdi_format_count == 1);
    assert(!(vdip->vdi_flags & VF_PER_F_Z0));
    z0_vector = vdip->vdi_z0_vector;
    if (vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE2) {
	(void)fprintf(fp, "[Version] 2.0\n");
    }
    switch (vfdp->vfd_parameter) {
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
    switch (vfdp->vfd_format) {
    case VNADATA_FORMAT_DB_ANGLE:
	format = "DB";
	break;
    case VNADATA_FORMAT_MAG_ANGLE:
	format = "MA";
	break;
    case VNADATA_FORMAT_REAL_IMAG:
	format = "RI";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)fprintf(fp, "# Hz %c %s R ", parameter_name, format);
    (void)print_value(fp, vdip->vdi_dprecision, /*plus=*/false, /*pad=*/false,
	    z0_touchstone);
    (void)fputc('\n', fp);
    if (vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE2) {
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
		(void)print_value(fp, vdip->vdi_dprecision, /*plus=*/false,
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
static const char vnadata_check_name[] = "vnadata_cksave";
static const char vnadata_fsave_name[] = "vnadata_fsave";
static const char vnadata_save_name[] = "vnadata_save";

/*
 * vnadata_save_common: common save routine
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 *   @function: function name (for error messages)
 */
static int vnadata_save_common(vnadata_t *vdp, FILE *fp, const char *filename,
	const char *function)
{
    vnadata_internal_t *vdip;
    vnadata_parameter_type_t type;
    int rows, ports, frequencies;
    bool promote_ts2 = false;
    int aprecision;
    int rc = -1;
    const double *frequency_vector;
    const double complex *z0_vector = NULL;
    double z0_touchstone = 50.0;
    vnadata_t *conversions[VPT_NTYPES];

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
    aprecision = MAX(vdip->vdi_dprecision, 3);

    /*
     * Init conversions to NULL.
     */
    (void)memset((void *)conversions, 0, sizeof(conversions));

    /*
     * Validate parameters.  These errors are for the application
     * developer, not the end user, so show the function name.
     */
    if (function == vnadata_fsave_name && fp == NULL) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"%s: error: NULL file pointer", function);
	goto out;
    }
    if (filename == NULL) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"%s: error: NULL filename", function);
	goto out;
    }

    /*
     * Get the characteristics of network parameter data and make sure
     * the parameter type is known.  Note that for VPT_ZIN, the data
     * are stored as a row vector, but they really represent a diagonal
     * matrix.
     */
    type = vdp->vd_type;
    if (type == VPT_UNDEF) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"%s: cannot save with unknown network parameter data type",
		function);
	goto out;
    }
    rows             = vdp->vd_rows;
    ports            = vdp->vd_columns;
    frequencies      = vdp->vd_frequencies;
    frequency_vector = vdp->vd_frequency_vector;

    /*
     * If we don't have at least one port, fail.
     */
    if (ports < 1) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"%s: invalid data dimensions: %d x %d",
		function, rows, ports);
	goto out;
    }
    if (frequencies == 0) {
	_vnadata_error(vdip, VNAERR_USAGE,
		"%s: at least one frequency is required for save",
		function);
	goto out;
    }

    /*
     * Get simple z0 information.
     */
    if (!(vdip->vdi_flags & VF_PER_F_Z0)) {
	z0_vector = vdip->vdi_z0_vector;
	z0_touchstone = creal(z0_vector[0]);
    }

    /*
     * Set the file type.  If we can determine that the type is Touchstone
     * or NPD from the filename, use the determined type.  There's a
     * special-case here in that we allow Touchstone 1 format to be
     * saved with a ".ts" suffix, but will auto-promote it to Touchstone
     * 2 below if necessary.  If we can't determine the filetype from
     * the filename, but already have a filetype, keep the current type.
     * If all else fails, default to the native NPD format.
     */
    {
	vnadata_filetype_t filetype = _vnadata_parse_filename(filename, NULL);

	if (filetype == VNADATA_FILETYPE_TOUCHSTONE2 &&
		vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE1) {
	    promote_ts2 = true;
	} else if (filetype != VNADATA_FILETYPE_AUTO) {
	    vdip->vdi_filetype = filetype;
	} else if (vdip->vdi_filetype == VNADATA_FILETYPE_AUTO) {
	    vdip->vdi_filetype = VNADATA_FILETYPE_NPD;
	}
    }

    /*
     * If no formats have been given, default to "ri".
     */
    if (vdip->vdi_format_count == 0) {
	if (_vnadata_set_simple_format(vdip, type,
		    VNADATA_FORMAT_REAL_IMAG) == -1) {
	    goto out;
	}
    }

    /*
     * Check the compatibility of filetype, parameter type and format.
     */
    switch (vdip->vdi_filetype) {
    case VNADATA_FILETYPE_TOUCHSTONE1:
    case VNADATA_FILETYPE_TOUCHSTONE2:
	/*
	 * Touchstone format supports only one format.
	 */
	if (vdip->vdi_format_count > 1) {
	    _vnadata_error(vdip, VNAERR_USAGE, "%s: "
		    "only a single format may be specified in "
		    "Touchstone file type", function);
	    goto out;
	}

	/*
	 * Touchstone format supports only S, Z, Y, H & G parameters
	 * and only DB, MA and RI formats.
	 */
	{
	    vnadata_format_descriptor_t *vfdp = &vdip->vdi_format_vector[0];
	    vnadata_parameter_type_t fptype = vfdp->vfd_parameter;
	    bool bad = false;

	    if (fptype == VPT_UNDEF) {
		fptype = type;
	    }
	    switch (fptype) {
	    case VPT_UNDEF:
		break;

	    case VPT_S:
	    case VPT_Z:
	    case VPT_Y:
	    case VPT_H:
	    case VPT_G:
		break;

	    default:
		bad = true;
	    }
	    switch (vfdp->vfd_format) {
	    case VNADATA_FORMAT_DB_ANGLE:
	    case VNADATA_FORMAT_MAG_ANGLE:
	    case VNADATA_FORMAT_REAL_IMAG:
		break;

	    default:
		bad = true;
	    }
	    if (bad) {
		_vnadata_error(vdip, VNAERR_USAGE, "%s: %s format "
			"cannot be saved in Touchstone file type",
			function, _vnadata_format_to_name(vfdp));
		goto out;
	    }
	}

	/*
	 * Touchstone doesn't support frequency-dependent system
	 * impedances.
	 */
	if (z0_vector == NULL) {
	    _vnadata_error(vdip, VNAERR_USAGE, "%s: "
		    "cannot save frequency-dependent reference impedances "
		    "in Touchstone file type", function);
	    goto out;
	}

	/*
	 * Touchstone requires references to be real and positive.
	 */
	for (int i = 0; i < ports; ++i) {
	    if (cimag(z0_vector[i]) != 0.0 || creal(z0_vector[i]) <= 0.0) {
		_vnadata_error(vdip, VNAERR_USAGE, "%s: "
			"references must be be real and positive in "
			"Touchstone file type", function);
		goto out;
	    }
	}

	/*
	 * Finished with Touchstone 2.	If Touchstone 1, apply yet
	 * more constraints.
	 */
	if (vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE2) {
	    break;
	}

	/*
	 * Touchstone 1 doesn't support more than four ports.
	 */
	if (ports > 4) {
	    if (promote_ts2) {
		vdip->vdi_filetype = VNADATA_FILETYPE_TOUCHSTONE2;
		break;
	    }
	    _vnadata_error(vdip, VNAERR_USAGE, "%s: "
		    "cannot save a system with more than four ports in "
		    "Touchstone 1 file type", function);
	    goto out;
	}

	/*
	 * Touchstone 1 format doesn't support ports having different
	 * system impedances.
	 */
	for (int i = 1; i < ports; ++i) {
	    if (z0_vector[i] != z0_vector[0]) {
		if (promote_ts2) {
		    vdip->vdi_filetype = VNADATA_FILETYPE_TOUCHSTONE2;
		    i = ports;
		    break;
		}
		_vnadata_error(vdip, VNAERR_USAGE, "%s: "
			"cannot save ports with different reference "
			"impedances in touchstone 1 format", function);
		goto out;
	    }
	}
	break;

    case VNADATA_FILETYPE_NPD:
	for (int i = 0; i < vdip->vdi_format_count; ++i) {
	    const vnadata_format_descriptor_t *vfdp =
		&vdip->vdi_format_vector[i];
	    vnadata_parameter_type_t fptype = vfdp->vfd_parameter;

	    /*
	     * For NPD format, disallow nonsensical use of dB
	     * for values that are neither power nor root power.
	     * Touchstone, on the other hand, allows this, defining dB
	     * as 20*log10(cabs(value)), even for parameters in units
	     * of ohms or seimens.
	     */
	    if (fptype == VPT_UNDEF) {
		fptype = type;
	    }
	    if (vfdp->vfd_format == VNADATA_FORMAT_DB_ANGLE &&
		    !_VNADATA_IS_POWER(fptype)) {
		_vnadata_error(vdip, VNAERR_USAGE, "%s: "
			"%s: in NPD format, only power or root-power "
			"parameters can be displayed in dB",
			function, _vnadata_format_to_name(vfdp));
		goto out;
	    }

	    /*
	     * If insertion loss is requested, make sure the matrix has
	     * off-diagonal elements.
	     */
	    if (vfdp->vfd_format == VNADATA_FORMAT_IL && ports < 2) {
		_vnadata_error(vdip, VNAERR_USAGE, "%s: "
			"return loss requires at least one "
			"off-diagonal element", function);
		goto out;
	    }
	}
	break;

    case VNADATA_FILETYPE_AUTO:
	abort();
	/*NOTREACHED*/
    }

    /*
     * For each parameter, make sure the parameter data is convertible
     * to the required type.
     */
    for (int i = 0; i < vdip->vdi_format_count; ++i) {
	const vnadata_format_descriptor_t *vfdp = &vdip->vdi_format_vector[i];
	vnadata_parameter_type_t fptype = vfdp->vfd_parameter;

	if (fptype == VPT_UNDEF) {
	    fptype = type;
	}
	if (_VNADATA_IS_MATRIX(fptype) && !_VNADATA_IS_MATRIX(type)) {
	    _vnadata_error(vdip, VNAERR_USAGE, "%s: cannot convert "
		    "%s parameters for format %s",
		    function, vnadata_get_type_name(type),
		    _vnadata_format_to_name(vfdp));
	    goto out;
	}
    }

    /*
     * If vnadata_cksave, we're done.
     */
    if (function == vnadata_check_name) {
	rc = 0;
	goto out;
    }

    /*
     * If touchstone 1, normalize all system impedances to 1.
     */
    if (vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE1 &&
	    z0_vector[0] != 1.0) {
	vnadata_t *vdp_copy;
	vnadata_parameter_type_t target_type;

	/*
	 * If the input type is S or T, make a writeable copy.	Otherwise,
	 * convert to S using the existing z0.
	 */
	if ((vdp_copy = vnadata_alloc(vdip->vdi_error_fn,
			vdip->vdi_error_arg)) == NULL) {
	    goto out;
	}
	switch (vnadata_get_type(vdp)) {
	case VPT_S:
	default:
	    target_type = VPT_S;
	    break;

	case VPT_T:
	    target_type = VPT_T;
	    break;

	case VPT_U:
	    target_type = VPT_U;
	    break;
	}
	if (vnadata_convert(vdp, vdp_copy, target_type) == -1) {
	    vnadata_free(vdp_copy);
	    goto out;
	}

	/*
	 * Set all z0's to 1.
	 */
	if (vnadata_set_all_z0(vdp_copy, 1.0) == -1) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "vnadata_set_all_z0: %s", strerror(errno));
	    vnadata_free(vdp_copy);
	    goto out;
	}

	/*
	 * Save the copy in the conversions array so that it gets freed
	 * at out.  Replace vdp and z0_vector.	Even if we need the
	 * original type, we have to convert it from out new matrix in
	 * order to normalize impedances.
	 */
	conversions[target_type] = vdp_copy;
	vdp = vdp_copy;
	type = vdp->vd_type;
	vdip = VDP_TO_VDIP(vdp);
	if (!(vdip->vdi_flags & VF_PER_F_Z0)) {
	    z0_vector = vdip->vdi_z0_vector;
	} else {
	    assert(z0_vector == NULL);
	}
    }

    /*
     * Perform all the needed conversions.
     */
    for (int i = 0; i < vdip->vdi_format_count; ++i) {
	vnadata_format_descriptor_t *vfdp = &vdip->vdi_format_vector[i];

	if (vfdp->vfd_parameter != VPT_UNDEF) {
	    if (convert_input(function, vdip, conversions,
			vfdp->vfd_parameter) == -1) {
		goto out;
	    }
	}
    }

    /*
     * Go through the format vector and fix up any instances of "ri",
     * "ma" and "db" without parameter types, taking the parameter type
     * from the vnadata_t structure.
     */
    {
	bool changed = false;

	for (int i = 0; i < vdip->vdi_format_count; ++i) {
	    vnadata_format_descriptor_t *vfdp = &vdip->vdi_format_vector[i];

	    if (vfdp->vfd_parameter == VPT_UNDEF) {
		vfdp->vfd_parameter = type;
		changed = true;
	    }
	}
	if (changed) {
	    if (_vnadata_update_format_string(vdip) == -1) {
		goto out;
	    }
	}
    }

    /*
     * If vnadata_save, open the output file.
     */
    if (function == vnadata_save_name) {
	if ((fp = fopen(filename, "w")) == NULL) {
	    _vnadata_error(vdip, VNAERR_SYSTEM, "fopen: %s: %s",
		    filename, strerror(errno));
	    goto out;
	}
    }

    /*
     * Print the file header.
     */
    switch (vdip->vdi_filetype) {
    case VNADATA_FILETYPE_TOUCHSTONE1:
    case VNADATA_FILETYPE_TOUCHSTONE2:
	print_touchstone_header(vdip, fp, z0_touchstone);
	break;

    case VNADATA_FILETYPE_NPD:
	print_npd_header(vdip, fp);
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
	print_value(fp, vdip->vdi_fprecision, /*plus=*/false, /*pad=*/true,
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
		print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
			/*pad=*/true, creal(z0));
		(void)fputc(' ', fp);
		print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
			/*pad=*/true, cimag(z0));
	    }
	}

	/*
	 * For each parameter...
	 */
	for (int format = 0; format < vdip->vdi_format_count; ++format) {
	    const vnadata_format_descriptor_t *vfdp =
		&vdip->vdi_format_vector[format];
	    const vnadata_t *matrix = NULL;
	    const double complex *data = NULL;
	    bool last_arg;
	    bool done = false;

	    /*
	     * Get the required matrix.
	     */
	    assert(vfdp->vfd_parameter != VPT_UNDEF);
	    if (vfdp->vfd_parameter == type) {
		matrix = vdp;
	    } else {
		matrix = conversions[vfdp->vfd_parameter];
	    }
	    assert(matrix != NULL);

	    /*
	     * Print
	     */
	    switch (vfdp->vfd_parameter) {
	    case VPT_S:
		switch (vfdp->vfd_format) {
		case VNADATA_FORMAT_IL:
		    for (int row = 0; row < rows; ++row) {
			for (int column = 0; column < ports; ++column) {
			    double complex value;

			    if (row == column) {
				continue;
			    }
			    value = vnadata_get_cell(matrix, findex, row,
				    column);
			    (void)fputc(' ', fp);
			    last_arg = format == vdip->vdi_format_count - 1 &&
				row == rows - 1 && column == ports - 1;
			    print_value(fp, vdip->vdi_dprecision,
				    /*plus=*/true, /*pad=*/!last_arg,
				    -20.0 * log10(cabs(value)));
			}
		    }
		    done = true;
		    break;

		case VNADATA_FORMAT_RL:
		    for (int port = 0; port < ports; ++port) {
			double complex value;

			value = vnadata_get_cell(matrix, findex, port, port);
			(void)fputc(' ', fp);
			last_arg = format == vdip->vdi_format_count - 1 &&
			    port == ports - 1;
			print_value(fp, vdip->vdi_dprecision,
				/*plus=*/true, /*pad=*/!last_arg,
				-20.0 * log10(cabs(value)));
		    }
		    done = true;
		    break;

		case VNADATA_FORMAT_VSWR:
		    for (int port = 0; port < ports; ++port) {
			double complex sxx;
			double a;
			double vswr;

			sxx = vnadata_get_cell(matrix, findex, port, port);
			a = cabs(sxx);
			vswr = (1.0 + a) / fabs(1.0 - a);
			(void)fputc(' ', fp);
			last_arg = format == vdip->vdi_format_count - 1 &&
			    port == ports - 1;
			print_value(fp, vdip->vdi_dprecision,
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

	    case VPT_T:
	    case VPT_U:
	    case VPT_Z:
	    case VPT_Y:
	    case VPT_H:
	    case VPT_G:
	    case VPT_A:
	    case VPT_B:
		for (int row = 0; row < rows; ++row) {
		    for (int column = 0; column < ports; ++column) {
			double complex value;
			vnadata_filetype_t filetype = vdip->vdi_filetype;

			/*
			 * In Touchstone formats, break the line after
			 * every four columns, and, except for two-port,
			 * after every row.
			 */
			if ((filetype == VNADATA_FILETYPE_TOUCHSTONE1 ||
			     filetype == VNADATA_FILETYPE_TOUCHSTONE2) &&
			     ((column  != 0 && column % 4 == 0) ||
			      (ports != 2 && row != 0 && column == 0))) {
			    (void)fputc('\n', fp);
			    for (int i = 0; i < vdip->vdi_fprecision + 5; ++i)
				(void)fputc(' ', fp);
			}

			/*
			 * Format based on vfd_format.
			 *
			 * Special case the 2x2 matrix in Touchstone 1
			 * which prints in column major order.
			 */
			if (filetype == VNADATA_FILETYPE_TOUCHSTONE1 &&
				ports == 2) {
			    assert(rows == ports);
			    value = vnadata_get_cell(matrix, findex,
				    column, row);
			} else {
			    value = vnadata_get_cell(matrix, findex,
				    row, column);
			}
			switch (vfdp->vfd_format) {
			case VNADATA_FORMAT_DB_ANGLE:
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/true, 20.0 * log10(cabs(value)));
			    if (aprecision == VNADATA_MAX_PRECISION) {
				(void)fprintf(fp, " %+a",
					180.0 / M_PI * carg(value));
			    } else {
				(void)fprintf(fp, " %+*.*f",
					aprecision + 4,
					aprecision - 1,
					180.0 / M_PI * carg(value));
			    }
			    break;

			case VNADATA_FORMAT_MAG_ANGLE:
			    (void)fprintf(fp, "  ");
			    print_value(fp, vdip->vdi_dprecision,
				    /*plus=*/false,
				    /*pad=*/true, cabs(value));
			    if (aprecision == VNADATA_MAX_PRECISION) {
				(void)fprintf(fp, " %+a",
					180.0 / M_PI * carg(value));
			    } else {
				(void)fprintf(fp, " %+*.*f",
					aprecision + 4,
					aprecision - 1,
					180.0 / M_PI * carg(value));
			    }
			    break;

			case VNADATA_FORMAT_REAL_IMAG:
			    last_arg = format == vdip->vdi_format_count - 1 &&
				row == rows - 1 && column == ports - 1;
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/true, creal(value));
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
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
		for (int port = 0; port < ports; ++port) {
		    double complex value;

		    value = data[port];
		    switch (vfdp->vfd_format) {
		    case VNADATA_FORMAT_MAG_ANGLE:
			(void)fprintf(fp, "  ");
			print_value(fp, vdip->vdi_dprecision, /*plus=*/false,
				/*pad=*/true, cabs(value));
			if (aprecision == VNADATA_MAX_PRECISION) {
			    (void)fprintf(fp, " %+a",
				    180.0 / M_PI * carg(value));
			} else {
			    (void)fprintf(fp, "  %+*.*f",
				    aprecision + 2,
				    aprecision - 3,
				    180.0 / M_PI * carg(value));
			}
			break;

		    case VNADATA_FORMAT_REAL_IMAG:
			last_arg = format == vdip->vdi_format_count - 1 &&
			    port == ports - 1;
			(void)fputc(' ', fp);
			print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				/*pad=*/true, creal(value));
			(void)fputc(' ', fp);
			print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				/*pad=*/!last_arg, cimag(value));
			break;

		    case VNADATA_FORMAT_PRC:
			{
			    double complex z;
			    double zr, zi;
			    double r, x, c;

			    z = data[port];
			    zr = creal(z);
			    zi = cimag(z);
			    r = (zr*zr + zi*zi) / zr;
			    x = (zr*zr + zi*zi) / zi;
			    c = -1.0 /
				(2.0 * M_PI * frequency_vector[findex] * x);
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/true, r);
			    (void)fputc(' ', fp);
			    last_arg = format == vdip->vdi_format_count - 1 &&
				port == ports - 1;
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, c);
			}
			break;

		    case VNADATA_FORMAT_PRL:
			{
			    double complex z;
			    double zr, zi;
			    double r, x, l;

			    z = data[port];
			    zr = creal(z);
			    zi = cimag(z);
			    r = (zr*zr + zi*zi) / zr;
			    x = (zr*zr + zi*zi) / zi;
			    l = x / (2.0 * M_PI * frequency_vector[findex]);
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/true, r);
			    (void)fputc(' ', fp);
			    last_arg = format == vdip->vdi_format_count - 1 &&
				port == ports - 1;
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, l);
			}
			break;

		    case VNADATA_FORMAT_SRC:
			{
			    double complex z;
			    double zr, zi;
			    double c;

			    z = data[port];
			    zr = creal(z);
			    zi = cimag(z);
			    c = -1.0 /
				(2.0 * M_PI * frequency_vector[findex] * zi);
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/true, zr);
			    (void)fputc(' ', fp);
			    last_arg = format == vdip->vdi_format_count - 1 &&
				port == ports - 1;
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/!last_arg, c);
			}
			break;

		    case VNADATA_FORMAT_SRL:
			{
			    double complex z;
			    double zr, zi;
			    double l;

			    z = data[port];
			    zr = creal(z);
			    zi = cimag(z);
			    l = zi / (2.0 * M_PI * frequency_vector[findex]);
			    (void)fputc(' ', fp);
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
				    /*pad=*/true, zr);
			    (void)fputc(' ', fp);
			    last_arg = format == vdip->vdi_format_count - 1 &&
				port == ports - 1;
			    print_value(fp, vdip->vdi_dprecision, /*plus=*/true,
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
    if (vdip->vdi_filetype == VNADATA_FILETYPE_TOUCHSTONE2) {
	(void)fprintf(fp, "[End]\n");
    }

    /*
     * If vnadata_save, close the output file.
     */
    if (function == vnadata_save_name) {
	if (fclose(fp) == -1) {
	    _vnadata_error(vdip, VNAERR_SYSTEM, "fclose: %s: %s",
		    filename, strerror(errno));
	    fp = NULL;
	    goto out;
	}
	fp = NULL;
    }
    rc = 0;

out:
    if (function == vnadata_save_name && fp != NULL) {
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
 * vnadata_cksave: check that the given parameters and format are valid for save
 *   @vdp: a pointer to the vnadata_t structure
 *   @filename: file to save
 */
int vnadata_cksave(vnadata_t *vdp, const char *filename)
{
    return vnadata_save_common(vdp, NULL, filename, vnadata_check_name);
}

/*
 * vnadata_fsave: save network parameters to a file pointer
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 */
int vnadata_fsave(vnadata_t *vdp, FILE *fp, const char *filename)
{
    return vnadata_save_common(vdp, fp, filename, vnadata_fsave_name);
}

/*
 * vnadata_save: save network parameters to filename
 *   @vdp: a pointer to the vnadata_t structure
 *   @filename: file to save
 */
int vnadata_save(vnadata_t *vdp, const char *filename)
{
    return vnadata_save_common(vdp, NULL, filename, vnadata_save_name);
}
