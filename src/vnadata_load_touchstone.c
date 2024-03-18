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

#define VNADATA_NO_BOUNDS_CHECK

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "vnadata_internal.h"


/*
 * Constants
 */
#define PI		3.14159265358979323846264338327950288419716939937508
#define LOG10		2.30258509299404568401799145468436420760110148862877
#define RAD_PER_DEG	(PI / 180.0)

/*
 * ts_token_t: touchstone tokens
 */
typedef enum ts_token {
    /*
     * Version 2 keywords
     */
    T_KW_BEGIN_INFORMATION,
    T_KW_END_INFORMATION,
    T_KW_MATRIX_FORMAT,
    T_KW_MIXED_MODE_ORDER,
    T_KW_NETWORK_DATA,
    T_KW_NOISE_DATA,
    T_KW_NUMBER_OF_FREQUENCIES,
    T_KW_NUMBER_OF_NOISE_FREQUENCIES,
    T_KW_NUMBER_OF_PORTS,
    T_KW_REFERENCE,
    T_KW_TWO_PORT_ORDER,
    T_KW_VERSION,
    T_KW_END,

    /*
     * Option keywords
     */
    T_OP_HZ,
    T_OP_KHZ,
    T_OP_MHZ,
    T_OP_GHZ,
    T_OP_THZ,	/* non-standard */
    T_OP_S,
    T_OP_Y,
    T_OP_Z,
    T_OP_H,
    T_OP_G,
    T_OP_DB,
    T_OP_MA,
    T_OP_RI,
    T_OP_R,

    /*
     * Basic tokens
     */
    T_OPTION,
    T_WORD,
    T_INT,
    T_DOUBLE,
    T_EOL,
    T_EOF,
    T_ERROR,
} ts_token_t;

/*
 * ts_flags_t: flags to next_token
 */
typedef enum ts_flags {
    F_NONE	= 0x0000,		/* no special flags */
    F_NOCONV	= 0x0001,		/* don't convert numbers */
    F_INT	= 0x0002,		/* convert integer */
    F_EOL	= 0x0004,		/* scan newlines */
} ts_flags_t;

/*
 * ts_parser_state_t: touchstone parser state
 */
typedef struct ts_parser_state {
    vnadata_internal_t *tps_vdip;
    FILE *tps_fp;
    const char *tps_filename;
    int tps_line;
    int tps_char;
    bool tps_in_option_line;
    ts_token_t tps_token;
    size_t tps_text_length;
    size_t tps_text_allocation;
    char *tps_text;
    union {
	int    tps_int;
	double tps_double;
    } u;
    double tps_frequency_multiplier;
    vnadata_parameter_type_t tps_parameter_type;
    char tps_data_format;			/* (D)B, (M)A, (R)I */
    double tps_z0;
    int tps_ports;
    size_t tps_value_count;
    size_t tps_value_allocation;
    double *tps_value_vector;
} ts_parser_state_t;

/*
 * get_token_name
 */
static const char *get_token_name(ts_parser_state_t *tpsp)
{
    switch (tpsp->tps_token) {
    case T_KW_BEGIN_INFORMATION:
	return "[Begin Information]";
    case T_KW_END_INFORMATION:
	return "[End Information]";
    case T_KW_MATRIX_FORMAT:
	return "[Matrix Format]";
    case T_KW_MIXED_MODE_ORDER:
	return "[Mixed-Mode Order]";
    case T_KW_NETWORK_DATA:
	return "[Network Data]";
    case T_KW_NOISE_DATA:
	return "[Noise Data]";
    case T_KW_NUMBER_OF_FREQUENCIES:
	return "[Number of Frequencies]";
    case T_KW_NUMBER_OF_NOISE_FREQUENCIES:
	return "[Number of Noise Frequencies]";
    case T_KW_NUMBER_OF_PORTS:
	return "[Number of Ports]";
    case T_KW_REFERENCE:
	return "[Reference]";
    case T_KW_TWO_PORT_ORDER:
	return "[Two-Port Order]";
    case T_KW_VERSION:
	return "[Version]";
    case T_KW_END:
	return "[End]";
    case T_OP_HZ:
	return "Hz";
    case T_OP_KHZ:
	return "KHz";
    case T_OP_MHZ:
	return "MHz";
    case T_OP_GHZ:
	return "GHz";
    case T_OP_THZ:
	return "THz";
    case T_OP_S:
	return "S";
    case T_OP_Y:
	return "Y";
    case T_OP_Z:
	return "Z";
    case T_OP_H:
	return "H";
    case T_OP_G:
	return "G";
    case T_OP_R:
	return "R";
    case T_OP_DB:
	return "DB";
    case T_OP_MA:
	return "MA";
    case T_OP_RI:
	return "RI";
    case T_OPTION:
	return "#";
    case T_WORD:
	return tpsp->tps_text;
    case T_INT:
	return tpsp->tps_text;
    case T_DOUBLE:
	return tpsp->tps_text;
    case T_EOL:
	return "<EOL>";
    case T_EOF:
	return "<EOF>";
    case T_ERROR:
	return "<ERROR>";
    default:
	break;
    }
    abort();
    /*NOTREACHED*/
}

/*
 * start_text: start accumulating text
 *   @tpsp: touchstone parser state structure
 */
static void start_text(ts_parser_state_t *tpsp)
{
    tpsp->tps_text_length = 0;
}

/*
 * add_char: add a character to the text buffer
 *   @tpsp: touchstone parser state structure
 *   @c: character to add
 *
 * Return:
 *	 0: success
 *	-1: out of memory error
 */
static int add_char(ts_parser_state_t *tpsp, char c)
{
    assert(tpsp->tps_text_allocation != 0);
    if (tpsp->tps_text_length + 1 >= tpsp->tps_text_allocation) {
	char *cp;
	size_t new_allocation = 2 * tpsp->tps_text_allocation;

	if ((cp = realloc(tpsp->tps_text, new_allocation)) == NULL) {
	    return -1;
	}
	tpsp->tps_text = cp;
	tpsp->tps_text_allocation = new_allocation;
    }
    tpsp->tps_text[tpsp->tps_text_length++] = c;
    return 0;
}

/*
 * end_text: stop accumulating text
 *   @tpsp: touchstone parser state structure
 */
static void end_text(ts_parser_state_t *tpsp)
{
    tpsp->tps_text[tpsp->tps_text_length] = '\000';
}

/*
 * next_char: read the next character from the input
 *   @tpsp: touchstone parser state structure
 */
static void next_char(ts_parser_state_t *tpsp)
{
    tpsp->tps_char = getc(tpsp->tps_fp);
    if (islower(tpsp->tps_char)) {
	tpsp->tps_char = toupper(tpsp->tps_char);
    }
}

/*
 * is_in_word_char: return true if ch can occur within a word
 *   @ch: character to test
 */
static int is_in_word_char(int ch)
{
    if (isalnum(ch)) {
	return 1;
    }
    switch (ch) {
    case '+':
    case ',':
    case '-':
    case '.':
    case '_':
	return 1;

    default:
	break;
    }
    return 0;
}

/*
 * convert_int: convert an integer
 *   @tpsp: touchstone parser state structure
 */
static bool convert_int(ts_parser_state_t *tpsp)
{
    char *end;

    tpsp->u.tps_int = strtol(tpsp->tps_text, &end, 0);
    return end > tpsp->tps_text && *end == '\000';
}

/*
 * convert_double: convert a double
 *   @tpsp: touchstone parser state structure
 */
static bool convert_double(ts_parser_state_t *tpsp)
{
    char *end;

    tpsp->u.tps_double = strtod(tpsp->tps_text, &end);
    return end > tpsp->tps_text && *end == '\000';
}

/*
 * next_token: scan the next token from the input
 *   @tpsp: touchstone parser state structure
 *   @flags: bitwise OR of ts_flags_t values
 */
static int next_token(ts_parser_state_t *tpsp, uint32_t flags)
{
    vnadata_internal_t *vdip = tpsp->tps_vdip;
    int saved_errno;

    for (;;) {
	switch (tpsp->tps_char) {
	case EOF:	/* end of file */
	    tpsp->tps_token = T_EOF;
	    return 0;

	case '\n':	/* newline */
	    ++tpsp->tps_line;
	    next_char(tpsp);
	    if ((flags & F_EOL) || tpsp->tps_in_option_line) {
		tpsp->tps_in_option_line = false;
		tpsp->tps_token = T_EOL;
		return 0;
	    }
	    continue;

	case '!':	/* comment */
	    do {
		next_char(tpsp);
	    } while (tpsp->tps_char != '\n' && tpsp->tps_char != EOF);
	    continue;

	case '+':	/* things that start "words" */
	case '-':
	case '.':
	    goto word;

	case '#':	/* start of option line */
	    next_char(tpsp);
	    tpsp->tps_in_option_line = true;
	    tpsp->tps_token = T_OPTION;
	    return 0;

	case '[':	/* start of keyword */
	    next_char(tpsp);
	    start_text(tpsp);
	    while (tpsp->tps_char != ']') {
		if (tpsp->tps_char == '\n' || tpsp->tps_char == EOF) {
		    break;
		}
		if (add_char(tpsp, (char)tpsp->tps_char) == -1) {
		    end_text(tpsp);
		    tpsp->tps_token = T_ERROR;
		    return -1;
		}
		next_char(tpsp);
	    }
	    end_text(tpsp);
	    if (tpsp->tps_char != ']') {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"missing closing brace of keyword",
			tpsp->tps_filename, tpsp->tps_line);
		tpsp->tps_token = T_ERROR;
		return -1;
	    }
	    next_char(tpsp);
	    switch (tpsp->tps_text_length) {
	    case 3:
		if (strcmp(tpsp->tps_text, "END") == 0) {
		    tpsp->tps_token = T_KW_END;
		    return 0;
		}
		break;

	    case 7:
		if (strcmp(tpsp->tps_text, "VERSION") == 0) {
		    tpsp->tps_token = T_KW_VERSION;
		    return 0;
		}
		break;

	    case 9:
		if (strcmp(tpsp->tps_text, "REFERENCE") == 0) {
		    tpsp->tps_token = T_KW_REFERENCE;
		    return 0;
		}
		break;

	    case 10:
		if (strcmp(tpsp->tps_text, "NOISE DATA") == 0) {
		    tpsp->tps_token = T_KW_NOISE_DATA;
		    return 0;
		}
		break;

	    case 12:
		if (strcmp(tpsp->tps_text, "NETWORK DATA") == 0) {
		    tpsp->tps_token = T_KW_NETWORK_DATA;
		    return 0;
		}
		break;

	    case 13:
		if (strcmp(tpsp->tps_text, "MATRIX FORMAT") == 0) {
		    tpsp->tps_token = T_KW_MATRIX_FORMAT;
		    return 0;
		}
		break;

	    case 14:
		if (strcmp(tpsp->tps_text, "TWO-PORT ORDER") == 0) {
		    tpsp->tps_token = T_KW_TWO_PORT_ORDER;
		    return 0;
		}
		break;

	    case 15:
		if (strcmp(tpsp->tps_text, "NUMBER OF PORTS") == 0) {
		    tpsp->tps_token = T_KW_NUMBER_OF_PORTS;
		    return 0;
		}
		if (strcmp(tpsp->tps_text, "END INFORMATION") == 0) {
		    tpsp->tps_token = T_KW_END_INFORMATION;
		    return 0;
		}
		break;

	    case 16:
		if (strcmp(tpsp->tps_text, "MIXED-MODE ORDER") == 0) {
		    tpsp->tps_token = T_KW_MIXED_MODE_ORDER;
		    return 0;
		}
		break;

	    case 17:
		if (strcmp(tpsp->tps_text, "BEGIN INFORMATION") == 0) {
		    tpsp->tps_token = T_KW_BEGIN_INFORMATION;
		    return 0;
		}
		break;

	    case 21:
		if (strcmp(tpsp->tps_text, "NUMBER OF FREQUENCIES") == 0) {
		    tpsp->tps_token = T_KW_NUMBER_OF_FREQUENCIES;
		    return 0;
		}
		break;

	    case 27:
		if (strcmp(tpsp->tps_text,
			    "NUMBER OF NOISE FREQUENCIES") == 0) {
		    tpsp->tps_token = T_KW_NUMBER_OF_NOISE_FREQUENCIES;
		    return 0;
		}
		break;

	    default:
		break;
	    }
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "unknown keyword [%s]",
		    tpsp->tps_filename, tpsp->tps_line,
		    tpsp->tps_text);
	    tpsp->tps_token = T_ERROR;
	    return 0;

	default:
	    if (isspace(tpsp->tps_char)) {	/* whitespace (not newline) */
		next_char(tpsp);
		continue;
	    }
	    break;
	}

	/*
	 * Scan words and numbers.
	 */
	if (isalnum(tpsp->tps_char)) {
	word:
	    start_text(tpsp);
	    do {
		if (add_char(tpsp, (char)tpsp->tps_char) == -1) {
		    end_text(tpsp);
		    tpsp->tps_token = T_ERROR;
		    return 0;
		}
		next_char(tpsp);
	    } while (is_in_word_char(tpsp->tps_char));
	    end_text(tpsp);

	    /*
	     * Convert numbers.
	     */
	    if (!(flags & F_NOCONV)) {
		if (flags & F_INT) {
		    if (convert_int(tpsp)) {
			tpsp->tps_token = T_INT;
			return 0;
		    }
		}
		if (convert_double(tpsp)) {
		    tpsp->tps_token = T_DOUBLE;
		    return 0;
		}
	    }

	    /*
	     * If we're on the option line, return option keywords.
	     */
	    if (tpsp->tps_in_option_line) {
		switch ((tpsp->tps_text_length << 8) | tpsp->tps_text[0]) {
		case 0x100 | 'G':
		    tpsp->tps_token = T_OP_G;
		    return 0;

		case 0x100 | 'H':
		    tpsp->tps_token = T_OP_H;
		    return 0;

		case 0x100 | 'R':
		    tpsp->tps_token = T_OP_R;
		    return 0;

		case 0x100 | 'S':
		    tpsp->tps_token = T_OP_S;
		    return 0;

		case 0x100 | 'Y':
		    tpsp->tps_token = T_OP_Y;
		    return 0;

		case 0x100 | 'Z':
		    tpsp->tps_token = T_OP_Z;
		    return 0;

		case 0x200 | 'D':
		    if (tpsp->tps_text[1] == 'B') {
			tpsp->tps_token = T_OP_DB;
			return 0;
		    }
		    break;

		case 0x200 | 'H':
		    if (tpsp->tps_text[1] == 'Z') {
			tpsp->tps_token = T_OP_HZ;
			return 0;
		    }
		    break;

		case 0x200 | 'M':
		    if (tpsp->tps_text[1] == 'A') {
			tpsp->tps_token = T_OP_MA;
			return 0;
		    }
		    break;

		case 0x200 | 'R':
		    if (tpsp->tps_text[1] == 'I') {
			tpsp->tps_token = T_OP_RI;
			return 0;
		    }
		    break;

		case 0x300 | 'G':
		    if (tpsp->tps_text[1] == 'H' && tpsp->tps_text[2] == 'Z') {
			tpsp->tps_token = T_OP_GHZ;
			return 0;
		    }
		    break;

		case 0x300 | 'K':
		    if (tpsp->tps_text[1] == 'H' && tpsp->tps_text[2] == 'Z') {
			tpsp->tps_token = T_OP_KHZ;
			return 0;
		    }
		    break;

		case 0x300 | 'M':
		    if (tpsp->tps_text[1] == 'H' && tpsp->tps_text[2] == 'Z') {
			tpsp->tps_token = T_OP_MHZ;
			return 0;
		    }
		    break;

		case 0x300 | 'T':
		    if (tpsp->tps_text[1] == 'H' && tpsp->tps_text[2] == 'Z') {
			tpsp->tps_token = T_OP_THZ;
			return 0;
		    }
		    break;
		}
	    }
	    tpsp->tps_token = T_WORD;
	    return 0;
	}

	/*
	 * Otherwise, we have an unexpected character.
	 */
	if (isprint(tpsp->tps_char)) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "unexpected character '%c'",
		    tpsp->tps_filename, tpsp->tps_line, tpsp->tps_char);
	} else {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "unexpected character '\\x%02x'",
		    tpsp->tps_filename, tpsp->tps_line, tpsp->tps_char);
	}
	saved_errno = errno;
	next_char(tpsp);
	errno = saved_errno;
	tpsp->tps_token = T_ERROR;
	return -1;
    }
}

/*
 * parse_data_line: parse a line of floating point numbers
 *   @tpsp: touchstone parser state structure
 */
static int parse_data_line(ts_parser_state_t *tpsp)
{
    vnadata_internal_t *vdip = tpsp->tps_vdip;
    tpsp->tps_value_count = 0;

    assert(tpsp->tps_token == T_DOUBLE);
    while (tpsp->tps_token == T_DOUBLE) {
	if (tpsp->tps_value_count >= tpsp->tps_value_allocation) {
	    size_t new_allocation;
	    double *lfp;

	    if (tpsp->tps_value_allocation == 0) {
		new_allocation = 9;
	    } else {
		new_allocation = 2 * tpsp->tps_value_allocation;
	    }
	    if ((lfp = realloc(tpsp->tps_value_vector,
			    new_allocation * sizeof(double))) == NULL) {
		_vnadata_error(vdip, VNAERR_SYSTEM,
			"realloc: %s", strerror(errno));
		return -1;
	    }
	    tpsp->tps_value_vector = lfp;
	    tpsp->tps_value_allocation = new_allocation;
	}
	tpsp->tps_value_vector[tpsp->tps_value_count++] = tpsp->u.tps_double;
	if (next_token(tpsp, F_EOL) == -1) {
	    return -1;
	}
    }
    if (tpsp->tps_token == T_EOL) {
	if (next_token(tpsp, F_NONE) == -1) {
	    return -1;
	}
	return 0;
    }
    if (tpsp->tps_token == T_EOF) {
	return 0;
    }
    _vnadata_error(vdip, VNAERR_SYNTAX,
	    "%s (line %d) error: unexpected token %s",
	    tpsp->tps_filename, tpsp->tps_line, get_token_name(tpsp));
    return -1;
}

/*
 * convert_value_pair: convert a pair of double to complex based on format
 *   @tpsp: touchstone parser state structure
 *   @value_pair: pointer to vector of two doubles
 *   @result: complex result
 */
static void convert_value_pair(ts_parser_state_t *tpsp,
	const double *value_pair, double complex *result)
{
    switch (tpsp->tps_data_format) {
    case 'D':	/* DB */
	*result = cexp(LOG10 * value_pair[0] / 20.0 +
			I * RAD_PER_DEG * value_pair[1]);
	break;

    case 'M':	/* MA */
	*result = value_pair[0] * cexp(I * RAD_PER_DEG * value_pair[1]);
	break;

    case 'R':	/* RI */
	*result = value_pair[0] + I * value_pair[1];
	break;

    default:
	abort();
    }
}

/*
 * load_touchstone1: load Touchstone version 1
 *   @tpsp: touchstone parser state structure
 */
static int load_touchstone1(ts_parser_state_t *tpsp)
{
    vnadata_internal_t *vdip = tpsp->tps_vdip;
    vnadata_t *vdp = &vdip->vdi_vd;
    int findex, row, column;
    bool maybe4ports = false;

    if (tpsp->tps_token != T_DOUBLE) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"expected a frequency value; found %s",
		tpsp->tps_filename, tpsp->tps_line,
		get_token_name(tpsp));
	return -1;
    }
    if (parse_data_line(tpsp) == -1)
	return -1;

    if (tpsp->tps_value_count % 2 == 0 || tpsp->tps_value_count < 3) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"first Touchstone V1 data line must have a odd number "
		"greater than 1 of fields",
	    tpsp->tps_filename, tpsp->tps_line);
	return -1;
    }
    if (tpsp->tps_value_count == 5)
	goto parse_noise_data;

    if (tpsp->tps_parameter_type == VPT_H ||
	    tpsp->tps_parameter_type == VPT_G) {
	if (tpsp->tps_value_count != 9) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected 9 fields; found %d",
		tpsp->tps_filename, tpsp->tps_line, (int)tpsp->tps_value_count);
	    return -1;
	}
	tpsp->tps_ports = 2;

    } else if (tpsp->tps_value_count == 9) {
	tpsp->tps_ports = 2;
	maybe4ports = true;

    } else {
	tpsp->tps_ports = (tpsp->tps_value_count - 1) / 2;
    }
    if (vnadata_init(vdp, tpsp->tps_parameter_type,
		tpsp->tps_ports, tpsp->tps_ports, /*frequencies*/0) == -1) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"realloc: %s", strerror(errno));
	return -1;
    }
    (void)vnadata_set_all_z0(vdp, tpsp->tps_z0);

    /*
     * 2x2
     */
    if (tpsp->tps_ports == 2) {
	for (;;) {
	    /*
	     * Validate and load the frequency.
	     */
	    findex = vdp->vd_frequencies;
	    if (findex != 0 &&
		    tpsp->tps_frequency_multiplier *
		    tpsp->tps_value_vector[0] <= vnadata_get_frequency(vdp,
			findex - 1)) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"frequencies must be in increasing order",
			tpsp->tps_filename, tpsp->tps_line);
		return -1;
	    }
	    if (vnadata_add_frequency(vdp, tpsp->tps_frequency_multiplier *
			tpsp->tps_value_vector[0]) == -1) {
		_vnadata_error(vdip, VNAERR_SYSTEM,
			"realloc: %s", strerror(errno));
		return -1;
	    }

	    /*
	     * Convert fields 1..8 and store transposed.
	     */
	    assert(tpsp->tps_value_count == 9);
	    convert_value_pair(tpsp, &tpsp->tps_value_vector[1],
		    &vdp->vd_data[findex][0]);
	    convert_value_pair(tpsp, &tpsp->tps_value_vector[3],
		    &vdp->vd_data[findex][2]);
	    convert_value_pair(tpsp, &tpsp->tps_value_vector[5],
		    &vdp->vd_data[findex][1]);
	    convert_value_pair(tpsp, &tpsp->tps_value_vector[7],
		    &vdp->vd_data[findex][3]);

	    /*
	     * If no more values, break.
	     */
	    if (tpsp->tps_token != T_DOUBLE)
		break;

	    /*
	     * Read the next line.
	     */
	    if (parse_data_line(tpsp) == -1)
		return -1;

	    /*
	     * If the new line has 9 fields, then our guess of a 2x2
	     * matrix was right.  If it has 8, then our guess was wrong --
	     * what we parsed so far was only the first line of a
	     * 4x4 matrix.  If it has 5 fields, then we've finished the
	     * data block and found noise parameters.  Any other number
	     * of fields is an error.
	     */
	    if (tpsp->tps_value_count != 9) {
		if (tpsp->tps_value_count == 5)
		    goto parse_noise_data;

		if (maybe4ports && tpsp->tps_value_count == 8) {
		    double complex d21, d12;

		    /*
		     * 4x4: transpose first row back into the correct order
		     * and resize the vnadata_t structure.
		     */
		    d21 = vdp->vd_data[findex][1];
		    d12 = vdp->vd_data[findex][2];
		    vdp->vd_data[findex][1] = d12;
		    vdp->vd_data[findex][2] = d21;
		    tpsp->tps_ports = 4;
		    if (vnadata_resize(vdp, tpsp->tps_parameter_type,
				tpsp->tps_ports, tpsp->tps_ports,
				findex + 1) == -1) {
			_vnadata_error(vdip, VNAERR_SYSTEM,
				"realloc: %s", strerror(errno));
			return -1;
		    }
		    (void)vnadata_set_all_z0(vdp, tpsp->tps_z0);
		    row = 1;
		    goto nxn_next_row;
		}
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected 9 fields; found %d",
		    tpsp->tps_filename, tpsp->tps_line,
		    (int)tpsp->tps_value_count);
		return -1;
	    }
	    maybe4ports = false;
	}

    /*
     * NxN (not 2x2)
     */
    } else {
	for (;;) {
	    /* first row */
	    findex = vdp->vd_frequencies;
	    if (findex != 0 &&
		    tpsp->tps_frequency_multiplier *
		    tpsp->tps_value_vector[0] <= vnadata_get_frequency(vdp,
			findex - 1)) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"frequencies must be in increasing order",
			tpsp->tps_filename, tpsp->tps_line);
		return -1;
	    }
	    if (vnadata_add_frequency(vdp, tpsp->tps_frequency_multiplier *
			tpsp->tps_value_vector[0]) == -1) {
		_vnadata_error(vdip, VNAERR_SYSTEM,
			"realloc: %s", strerror(errno));
		return -1;
	    }
	    for (column = 0; column < tpsp->tps_ports; ++column) {
		int idx = 1 + 2 * column;

		convert_value_pair(tpsp, &tpsp->tps_value_vector[idx],
			&vdp->vd_data[findex][column]);
	    }

	    /* remaining rows */
	    for (row = 1; row < tpsp->tps_ports; ++row) {
		if (tpsp->tps_token != T_DOUBLE) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "unexpected token %s",
			    tpsp->tps_filename, tpsp->tps_line,
			    get_token_name(tpsp));
		    return -1;
		}
		if (parse_data_line(tpsp) == -1)
		    return -1;

	    nxn_next_row:
		if (tpsp->tps_value_count != 2 * tpsp->tps_ports) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "expected %d fields; found %d",
			tpsp->tps_filename, tpsp->tps_line,
			2 * tpsp->tps_ports, (int)tpsp->tps_value_count);
		    return -1;
		}
		for (column = 0; column < tpsp->tps_ports; ++column) {
		    int idx = 2 * column;

		    convert_value_pair(tpsp, &tpsp->tps_value_vector[idx],
			    &vdp->vd_data[findex][tpsp->tps_ports * row +
			    column]);
		}
	    }
	    if (tpsp->tps_token != T_DOUBLE)
		break;

	    if (parse_data_line(tpsp) == -1)
		return -1;

	    if (tpsp->tps_value_count != 1 + 2 * tpsp->tps_ports) {
		if (tpsp->tps_value_count == 5)
		    goto parse_noise_data;

		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected %d fields; found %d",
		    tpsp->tps_filename, tpsp->tps_line,
		    1 + 2 * tpsp->tps_ports, (int)tpsp->tps_value_count);
		return -1;
	    }
	}
    }
    return 0;

    /*
     * Parse and discard noise data.  TODO: extend vnadata_t to store it.
     */
parse_noise_data:
    while (tpsp->tps_token == T_DOUBLE) {
	if (parse_data_line(tpsp) == -1)
	    return -1;

	if (tpsp->tps_value_count != 5) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected 5 noise fields; found %d",
		tpsp->tps_filename, tpsp->tps_line,
		(int)tpsp->tps_value_count);
	    return -1;
	}
    }
    return 0;
}

#define VNADATA_LOAD_INITIAL_TEXT_ALLOCATION	64

/*
 * Two-port order
 */
#define T12_21	1
#define T21_12	2

/*
 * parse_value_pair: parse two values and convert to complex
 *   @tpsp: touchstone parser state structure
 *   @nexpected: number of value pairs expected (for error messages)
 *   @result: address to receive complex result
 */
static int parse_value_pair(ts_parser_state_t *tpsp, int nexpected,
	double complex *result)
{
    vnadata_internal_t *vdip = tpsp->tps_vdip;
    double v[2];

    for (int i = 0; i < 2; ++i) {
	if (tpsp->tps_token != T_DOUBLE) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected %d value pairs",
		    tpsp->tps_filename, tpsp->tps_line, nexpected);
	    return -1;
	}
	v[i] = tpsp->u.tps_double;
	if (next_token(tpsp, F_NONE) == -1) {
	    return -1;
	}
    }
    convert_value_pair(tpsp, v, result);
    return 0;
}

/*
 * _vnadata_load_touchstone
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 */
int _vnadata_load_touchstone(vnadata_internal_t *vdip, FILE *fp,
	const char *filename)
{
    vnadata_t *vdp = &vdip->vdi_vd;
    ts_parser_state_t tps;
    int version = 1;
    int two_port_order  = -1;
    int two_port_order_line = -1;
    int number_of_frequencies = -1;
    int number_of_noise_frequencies = -1;
    int expected_pairs = 0;
    char matrix_format = 'F';	/* (F)ull, (L)ower, (U)pper */
    double complex *reference = NULL;
    int rc = -1;

    /*
     * Initialize the parser
     */
    (void)memset((void *)&tps, 0, sizeof(tps));
    tps.tps_vdip			= vdip;
    tps.tps_fp				= fp;
    tps.tps_filename			= filename;
    tps.tps_line			= 1;
    tps.tps_char			= '\000';
    tps.tps_in_option_line		= false;
    tps.tps_token			= T_EOL;
    tps.tps_text_length			= 0;
    tps.tps_text_allocation		= 0;
    tps.tps_text			= NULL;
    tps.tps_frequency_multiplier	= 1.0e+9;
    tps.tps_parameter_type		= VPT_S;
    tps.tps_data_format			= 'M';
    tps.tps_z0				= 50.0;	/* Touchstone default */
    tps.tps_ports			= -1;
    if ((tps.tps_text = malloc(VNADATA_LOAD_INITIAL_TEXT_ALLOCATION)) == NULL) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	goto out;
    }
    tps.tps_text_allocation = VNADATA_LOAD_INITIAL_TEXT_ALLOCATION;
    next_char(&tps);
    if (next_token(&tps, F_NONE) == -1)
	goto out;

    /*
     * Parse the [Version] line if present.
     */
    if (tps.tps_token == T_KW_VERSION) {
	double value;
	char *endptr;

	if (next_token(&tps, F_NOCONV) == -1) {
	    goto out;
	}
	if (tps.tps_token != T_WORD) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected version number; found %s",
		tps.tps_filename, tps.tps_line, get_token_name(&tps));
	    goto out;
	}
	value = strtod(tps.tps_text, &endptr);
	if (endptr != tps.tps_text && *endptr != '\000') {
	    value = 0.0;	/* invalid if, e.g. "2.0a" */
	}
	if (value == 2.0) {
	    version = 2;

	} else if (value >= 1.0 && value < 2.0) {
	    _vnadata_error(vdip, VNAERR_WARNING, "%s (line %d) warning: "
		    "Touchstone file contains dubious [Version] 1.x line",
		    tps.tps_filename, tps.tps_line);
	    version = 1;

	} else {
	    _vnadata_error(vdip, VNAERR_VERSION, "%s (line %d) error: "
		    "unsupported Touchstone version %s",
		    tps.tps_filename, tps.tps_line, tps.tps_text);
	    goto out;
	}
	if (next_token(&tps, F_NONE) == -1) {
	    goto out;
	}
    }

    /*
     * Parse the option line.
     */
    if (tps.tps_token != T_OPTION) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"expected # option line; found %s",
	    tps.tps_filename, tps.tps_line, get_token_name(&tps));
	goto out;
    }
    if (next_token(&tps, F_NONE) == -1) {
	goto out;
    }
    while (tps.tps_token != T_EOL) {
	switch (tps.tps_token) {
	case T_OP_HZ:
	    tps.tps_frequency_multiplier = 1.0;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_KHZ:
	    tps.tps_frequency_multiplier = 1.0e+3;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_MHZ:
	    tps.tps_frequency_multiplier = 1.0e+6;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_GHZ:
	    tps.tps_frequency_multiplier = 1.0e+9;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_THZ:	/* non-standard */
	    tps.tps_frequency_multiplier = 1.0e+12;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_S:
	    tps.tps_parameter_type = VPT_S;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_Y:
	    tps.tps_parameter_type = VPT_Y;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_Z:
	    tps.tps_parameter_type = VPT_Z;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_H:
	    tps.tps_parameter_type = VPT_H;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_G:
	    tps.tps_parameter_type = VPT_G;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_DB:
	    tps.tps_data_format = 'D';
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_MA:
	    tps.tps_data_format = 'M';
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_RI:
	    tps.tps_data_format = 'R';
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_OP_R:
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    if (tps.tps_token != T_DOUBLE) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected an impedance value after R",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    tps.tps_z0 = tps.u.tps_double;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_EOF:
	    break;

	default:
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "unexpected token \"%s\" in option line",
		tps.tps_filename, tps.tps_line, get_token_name(&tps));
	    goto out;
	}
    }
    if (tps.tps_token == T_EOL) {
	if (next_token(&tps, F_NONE) == -1) {
	    goto out;
	}
    }

    /*
     * Parse additional V2 keywords.
     */
    for (;;) {
	switch (tps.tps_token) {
	case T_KW_NUMBER_OF_PORTS:
	    if (next_token(&tps, F_INT) == -1) {
		goto out;
	    }
	    if (tps.tps_token != T_INT || tps.u.tps_int < 0) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected a positive integer after [Number of Ports]",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    tps.tps_ports = tps.u.tps_int;
	    if (tps.tps_ports != 2 &&
		    (tps.tps_parameter_type == VPT_G ||
		     tps.tps_parameter_type == VPT_H)) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"parameter type %s is incompatible with "
			"[Number of Ports] %d",
			tps.tps_filename, tps.tps_line,
			vnadata_get_type_name(tps.tps_parameter_type),
			tps.tps_ports);
		goto out;
	    }
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_KW_TWO_PORT_ORDER:
	    two_port_order_line = tps.tps_line;
	    if (next_token(&tps, F_NOCONV) == -1) {
		goto out;
	    }
	    if (tps.tps_token == T_WORD &&
		    strcmp(tps.tps_text, "12_21") == 0) {
		two_port_order = T12_21;
	    } else if (tps.tps_token == T_WORD &&
		    strcmp(tps.tps_text, "21_12") == 0) {
		two_port_order = T21_12;
	    } else {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected 12_21 or 21_12 after [Two-Port Order]",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_KW_NUMBER_OF_FREQUENCIES:
	    if (next_token(&tps, F_INT) == -1) {
		goto out;
	    }
	    if (tps.tps_token != T_INT) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected a positive integer after "
			"[Number of Frequencies]",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    number_of_frequencies = tps.u.tps_int;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_KW_NUMBER_OF_NOISE_FREQUENCIES:
	    if (next_token(&tps, F_INT) == -1) {
		goto out;
	    }
	    if (tps.tps_token != T_INT || tps.u.tps_int < 0) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected a positive integer after "
			"[Number of Noise Frequencies]",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    number_of_noise_frequencies = tps.u.tps_int;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_KW_REFERENCE:
	    if (tps.tps_ports < 0) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"[Number of Ports] must appear before [Reference]",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    if ((reference = calloc(tps.tps_ports,
			    sizeof(double complex))) == NULL) {
		_vnadata_error(vdip, VNAERR_SYSTEM,
			"calloc: %s", strerror(errno));
		goto out;
	    }
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    for (int i = 0; i < tps.tps_ports; ++i) {
		if (tps.tps_token != T_DOUBLE) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "expected %d values(s) after [Reference]",
			tps.tps_filename, tps.tps_line, tps.tps_ports);
		    goto out;
		}
		reference[i] = tps.u.tps_double;
		if (next_token(&tps, F_NONE) == -1) {
		    goto out;
		}
	    }
	    continue;

	case T_KW_MATRIX_FORMAT:
	    if (next_token(&tps, F_NOCONV) == -1) {
		goto out;
	    }
	    if (tps.tps_token == T_WORD &&
		    strcmp(tps.tps_text, "FULL") == 0) {
		matrix_format = 'F';
	    } else if (tps.tps_token == T_WORD &&
		    strcmp(tps.tps_text, "UPPER") == 0) {
		matrix_format = 'U';
	    } else if (tps.tps_token == T_WORD &&
		    strcmp(tps.tps_text, "LOWER") == 0) {
		matrix_format = 'L';
	    } else {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected Full, Upper or Lower after [Matrix Format]",
		    tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    continue;

	case T_KW_MIXED_MODE_ORDER:
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "[Mixed-Mode Order] not yet supported",
		tps.tps_filename, tps.tps_line);
	    goto out;

	case T_KW_BEGIN_INFORMATION:
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    if (tps.tps_token == T_KW_END_INFORMATION) {	/* optional */
		if (next_token(&tps, F_NONE) == -1) {
		    goto out;
		}
	    }
	    continue;

	default:
	    break;
	}
	break;
    }

    /*
     * Update the vnadata structure.
     */
    {
	vnadata_format_t format_type;

	/* set the file type */
	vdip->vdi_filetype = (version == 2) ?
	    VNADATA_FILETYPE_TOUCHSTONE2 : VNADATA_FILETYPE_TOUCHSTONE1;

	/* set the format string */
	switch (tps.tps_data_format) {
	case 'D':
	    format_type = VNADATA_FORMAT_DB_ANGLE;
	    break;
	case 'M':
	    format_type = VNADATA_FORMAT_MAG_ANGLE;
	    break;
	case 'R':
	    format_type = VNADATA_FORMAT_REAL_IMAG;
	    break;
	default:
	    abort();
	}
	if (_vnadata_set_simple_format(vdip, tps.tps_parameter_type,
		    format_type) == -1) {
	    _vnadata_error(vdip, VNAERR_SYSTEM,
		    "malloc: %s", strerror(errno));
	    goto out;
	}
    }

    /*
     * If version 1, call the V1 parser.  We've seen examples of files
     * that begin with an illegal [Version] 1.0 keyword followed by
     * other V2 keywords.  We'll tolerate those and parse V1 using the
     * V2 parser if all the required V2 keywords are present.  The most
     * significant difference between this hybrid format from V2 is that
     * we unnormalize the data at the end.
     */
    if (version == 1 && tps.tps_ports == -1 && number_of_frequencies == -1
	    && two_port_order == -1) {
	rc = load_touchstone1(&tps);
	goto expect_eof;
    }

    /*
     * Expect [Network Data].
     */
    if (tps.tps_token != T_KW_NETWORK_DATA) {
	_vnadata_error(vdip, VNAERR_SYNTAX,
		"%s (line %d) error: unexpected token %s",
	    tps.tps_filename, tps.tps_line, get_token_name(&tps));
	goto out;
    }
    if (next_token(&tps, F_NONE) == -1) {
	goto out;
    }

    /*
     * Make sure all required parameters were given.
     */
    if (tps.tps_ports < 0) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"[Number of Ports] must appear before [Network Data]",
		tps.tps_filename, tps.tps_line);
	goto out;
    }
    if (number_of_frequencies < 0) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"[Number of Frequencies] must appear before [Network Data]",
		tps.tps_filename, tps.tps_line);
	goto out;
    }
    if (tps.tps_ports == 2 && two_port_order == -1) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"[Two-Port Order] must appear before [Network Data]",
		tps.tps_filename, tps.tps_line);
	goto out;

    } else if (tps.tps_ports != 2 && two_port_order != -1) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"[Two-Port Order] may not be used with [Number of Ports] %d",
		tps.tps_filename, two_port_order_line, tps.tps_ports);
	goto out;
    }

    /*
     * Set up the output matrix.
     */
    if (vnadata_init(vdp, tps.tps_parameter_type, tps.tps_ports,
		tps.tps_ports, number_of_frequencies) == -1) {
	_vnadata_error(vdip, VNAERR_SYSTEM,
		"realloc: %s", strerror(errno));
	goto out;
    }

    /*
     * Set the reference impedances.
     */
    if (reference != NULL) {
	(void)vnadata_set_z0_vector(vdp, reference);
    } else {
	(void)vnadata_set_all_z0(vdp, tps.tps_z0);
    }

    /*
     * Parse [Network Data]
     */
    if (matrix_format == 'F') {
	expected_pairs = tps.tps_ports * tps.tps_ports;
    } else {
	expected_pairs = tps.tps_ports * (tps.tps_ports + 1) / 2;
    }
    for (int findex = 0; findex < number_of_frequencies; ++findex) {
	if (tps.tps_token != T_DOUBLE) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected frequency",
		    tps.tps_filename, tps.tps_line);
	    goto out;
	}
	if (findex != 0 &&
		tps.u.tps_double <= vnadata_get_frequency(vdp, findex - 1)) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "frequencies must be in increasing order",
		    tps.tps_filename, tps.tps_line);
	    goto out;
	}
	(void)vnadata_set_frequency(vdp, findex,
		tps.tps_frequency_multiplier * tps.u.tps_double);
	if (next_token(&tps, F_NONE) == -1) {
	    goto out;
	}
	switch (matrix_format) {
	case 'F':	/* Full */
	    for (int row = 0; row < tps.tps_ports; ++row) {
		for (int column = 0; column < tps.tps_ports; ++column) {
		    double complex x;

		    if (parse_value_pair(&tps, expected_pairs, &x) == -1) {
			goto out;
		    }
		    if (two_port_order == T21_12) {
			(void)vnadata_set_cell(vdp, findex, column, row, x);
		    } else {
			(void)vnadata_set_cell(vdp, findex, row, column, x);
		    }
		}
	    }
	    break;

	case 'U':	/* Upper */
	    for (int row = 0; row < tps.tps_ports; ++row) {
		for (int column = row; column < tps.tps_ports; ++column) {
		    double complex x;

		    if (parse_value_pair(&tps, expected_pairs, &x) == -1) {
			goto out;
		    }
		    (void)vnadata_set_cell(vdp, findex, row, column, x);
		    (void)vnadata_set_cell(vdp, findex, column, row, x);
		}
	    }
	    break;

	case 'L':	/* Lower */
	    for (int row = 0; row < tps.tps_ports; ++row) {
		for (int column = 0; column <= row; ++column) {
		    double complex x;

		    if (parse_value_pair(&tps, expected_pairs, &x) == -1) {
			goto out;
		    }
		    (void)vnadata_set_cell(vdp, findex, row, column, x);
		    (void)vnadata_set_cell(vdp, findex, column, row, x);
		}
	    }
	    break;

	default:
	    abort();
	}
    }

    /*
     * Parse and discard noise data.
     */
    if (number_of_noise_frequencies >= 0) {
	double f_prev = -1.0;

	if (tps.tps_token != T_KW_NOISE_DATA) {
	    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		    "expected [Noise Data]",
		    tps.tps_filename, tps.tps_line);
	    goto out;
	}
	if (next_token(&tps, F_NONE) == -1) {
	    goto out;
	}
	for (int i = 0; i < number_of_noise_frequencies; ++i) {
	    if (tps.tps_token != T_DOUBLE || tps.u.tps_double < 0.0) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"expected non-negative noise frequency",
			tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    if (i > 0 && tps.u.tps_double < f_prev) {
		_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			"noise frequencies be increasing",
			tps.tps_filename, tps.tps_line);
		goto out;
	    }
	    f_prev = tps.u.tps_double;
	    if (next_token(&tps, F_NONE) == -1) {
		goto out;
	    }
	    for (int j = 0; j < 4; ++j) {
		if (tps.tps_token != T_DOUBLE) {
		    _vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
			    "expected five noise parameters",
			    tps.tps_filename, tps.tps_line);
		    goto out;
		}
		if (next_token(&tps, F_NONE) == -1) {
		    goto out;
		}
	    }
	}
    }

    /*
     * Expect [End]
     */
    if (tps.tps_token == T_KW_END) {
	if (next_token(&tps, F_NONE) == -1)
	    goto out;

    } else {
	_vnadata_error(vdip, VNAERR_WARNING,
		"%s (line %d) warning: expected [End] keyword",
		tps.tps_filename, tps.tps_line);
    }

    /*
     * Expect <EOF>
     */
expect_eof:
    if (tps.tps_token != T_EOF) {
	_vnadata_error(vdip, VNAERR_SYNTAX, "%s (line %d) error: "
		"extra token(s) at end of file: %s",
		tps.tps_filename, tps.tps_line, get_token_name(&tps));
	goto out;
    }

    /*
     * If V1, unnormalize the data.
     */
    if (version == 1) {
	switch (tps.tps_parameter_type) {
	case VPT_S:
	default:
	    break;

	case VPT_Z:
	    for (int findex = 0; findex < vdp->vd_frequencies; ++findex) {
		for (int cell = 0; cell < tps.tps_ports * tps.tps_ports;
			++cell) {
		    vdp->vd_data[findex][cell] *= tps.tps_z0;
		}
	    }
	    break;

	case VPT_Y:
	    for (int findex = 0; findex < vdp->vd_frequencies; ++findex) {
		for (int cell = 0; cell < tps.tps_ports * tps.tps_ports;
			++cell) {
		    vdp->vd_data[findex][cell] /= tps.tps_z0;
		}
	    }
	    break;

	case VPT_H:
	    for (int findex = 0; findex < vdp->vd_frequencies; ++findex) {
		vdp->vd_data[findex][0] *= tps.tps_z0;
		vdp->vd_data[findex][3] /= tps.tps_z0;
	    }
	    break;

	case VPT_G:
	    for (int findex = 0; findex < vdp->vd_frequencies; ++findex) {
		vdp->vd_data[findex][0] /= tps.tps_z0;
		vdp->vd_data[findex][3] *= tps.tps_z0;
	    }
	    break;
	}
    }
    rc = 0;

out:
    free((void *)reference);
    free((void *)tps.tps_text);
    free((void *)tps.tps_value_vector);
    return rc;
}
