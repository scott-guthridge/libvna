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

#include <arpa/inet.h>	/* for htonl, ntohl */
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

/*
 * matrix_id_t: which error term matrix
 */
typedef enum matrix_id {
    E,
    EL, ER, EM,
    TS, TI, TX, TM,
    UM, UI, UX, US,
    MATRIX_IDS
} matrix_id_t;

/*
 * matrix_names: names of the matrices
 */
static const char *matrix_names[] = {
    "e",
    "el", "em", "er",
    "ts", "ti", "tx", "tm",
    "um", "ui", "ux", "us",
};

#define E_MASK	((1U << EL) | (1U << ER) | (1U << EM))
#define T_MASK	((1U << TS) | (1U << TI) | (1U << TX) | (1U << TM))
#define U_MASK	((1U << UM) | (1U << UI) | (1U << UX) | (1U << US))

/*
 * vnacal_load_state_t: parser state
 */
typedef struct vnacal_load_state {
    vnacal_t *vls_vcp;
    yaml_document_t vls_document;
    int vls_major_version;
    int vls_minor_version;
    int vls_findex;
} vnacal_load_state_t;

/*
 * parse_int: parse an integer
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: yaml node containing text
 *   @result: returned integer
 */
static int parse_int(vnacal_load_state_t *vlsp,
	yaml_node_t *node, int *result)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    char extra;

    if (node->type != YAML_SCALAR_NODE) {
	goto error;
    }
    if (sscanf((const char *)node->data.scalar.value, "%d %c",
		result, &extra) != 1) {
	goto error;
    }
    return 0;

error:
    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: expected integer",
	    vcp->vc_filename, node->start_mark.line + 1);
    return -1;
}

/*
 * parse_double: parse a double
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: yaml node containing text
 *   @result: returned double
 */
static int parse_double(vnacal_load_state_t *vlsp,
	yaml_node_t *node, double *result)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    char extra;

    if (node->type != YAML_SCALAR_NODE) {
	goto error;
    }
    if (sscanf((const char *)node->data.scalar.value, "%lf %c",
		result, &extra) != 1) {
	goto error;
    }
    return 0;

error:
    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: expected number",
	    vcp->vc_filename, node->start_mark.line + 1);
    return -1;
}

/*
 * parse_complex: parse a complex number
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: yaml node containing text
 *   @result: returned double
 */
static int parse_complex(vnacal_load_state_t *vlsp,
	yaml_node_t *node, double complex *result)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    const char *cur;
    char *end;
    double value1 = 0.0, value2 = 0.0;
    int code = 0;

    if (node->type != YAML_SCALAR_NODE) {
	goto error;
    }
    cur = (const char *)node->data.scalar.value;
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
	goto error;
    }
    switch (code) {
    case 1:	/* number */
	*result = value1;
	break;
    case 4:	/* j */
	*result = I;
	break;
    case 5:	/* number j */
	*result = value1 * I;
	break;
    case 6:	/* number number j */
	*result = value1 + value2 * I;
	break;
    case 12:	/* +j */
	*result = I;
	break;
    case 13:	/* number + j */
	*result = value1 + I;
	break;
    case 20:	/* -j */
	*result = -I;
	break;
    case 21:	/* number - j */
	*result = value1 - I;
	break;
    default:
	goto error;
    }
    return 0;

error:
    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: expected number",
	    vcp->vc_filename, node->start_mark.line + 1);
    return -1;
}

/*
 * parse_type:
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: yaml node containing text
 *   @result: returned double
 */
static int parse_type(vnacal_load_state_t *vlsp,
	yaml_node_t *node, vnacal_type_t *result)
{
    vnacal_t *vcp = vlsp->vls_vcp;

    if (node->type != YAML_SCALAR_NODE) {
	goto error;
    }
    *result = vnacal_name_to_type((const char *)node->data.scalar.value);
    if (*result == (vnacal_type_t)-1) {
	goto error;
    }
    return 0;

error:
    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: expected type",
	    vcp->vc_filename, node->start_mark.line + 1);
    return -1;
}


/*
 * parse_properties: parse user properties
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: user properties
 */
static vnaproperty_t *parse_properties(vnacal_load_state_t
	*vlsp, yaml_node_t *node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    vnaproperty_t *root = NULL;
    vnaproperty_t *subtree = NULL;

    switch (node->type) {
    case YAML_SCALAR_NODE:
	if ((root = vnaproperty_scalar_alloc((const char *)node->data.
			scalar.value)) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "vnaproperty_scalar_alloc: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto out;
	}
	return root;

    case YAML_MAPPING_NODE:
	{
	    yaml_node_pair_t *pair;

	    if ((root = vnaproperty_map_alloc()) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"vnaproperty_map_alloc: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto out;
	    }
	    for (pair = node->data.mapping.pairs.start;
		 pair < node->data.mapping.pairs.top; ++pair) {
		yaml_node_t *key, *value;

		key = yaml_document_get_node(&vlsp->vls_document, pair->key);
		if (key->type != YAML_SCALAR_NODE) {
		    _vnacal_error(vcp, VNAERR_WARNING,
			    "%s (line %ld) warning: "
			    "non-scalar property key ignored\n",
			    vcp->vc_filename, key->start_mark.line + 1);
		    continue;
		}
		value = yaml_document_get_node(&vlsp->vls_document,
			pair->value);
		subtree = parse_properties(vlsp, value);
		if (subtree == NULL) {
		    goto out;
		}
		if (vnaproperty_map_set(root,
			    (const char *)key->data.scalar.value,
			    subtree) == -1) {
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "vnaproperty_map_set: %s: %s",
			    vcp->vc_filename, strerror(errno));
		    goto out;
		}
		subtree = NULL;
	    }
	    return root;
	}

    case YAML_SEQUENCE_NODE:
	{
	    yaml_node_item_t *item;

	    if ((root = vnaproperty_list_alloc()) == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"vnaproperty_list_alloc: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto out;
	    }
	    for (item = node->data.sequence.items.start;
		 item < node->data.sequence.items.top; ++item) {
		int term = item - node->data.sequence.items.start;
		yaml_node_t *value;

		value = yaml_document_get_node(&vlsp->vls_document, *item);
		subtree = parse_properties(vlsp, value);
		if (subtree == NULL) {
		    goto out;
		}
		if (vnaproperty_list_set(root, term, subtree) == -1) {
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "vnaproperty_list_set: %s: %s",
			    vcp->vc_filename, strerror(errno));
		    goto out;
		}
	    }
	    return root;
	}

    default:
	abort();
    }

out:
    vnaproperty_free(subtree);
    vnaproperty_free(root);
    return NULL;
}

/*
 * parse_vector: parse a vector
 *   @vlsp: pointer to vnacal_load info structure
 *   @vector: where to place data
 *   @length: expected length
 *   @node: yaml node to parse
 */
static int parse_vector(vnacal_load_state_t *vlsp, double complex **vector,
	int length, yaml_node_t *node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    int findex = vlsp->vls_findex;
    yaml_node_item_t *item;
    int items;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected sequence",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    items = node->data.sequence.items.top - node->data.sequence.items.start;
    if (items != length) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected %d terms but found %d",
		vcp->vc_filename, node->start_mark.line + 1, length, items);
	return -1;
    }
    for (item = node->data.sequence.items.start;
	 item < node->data.sequence.items.top; ++item) {
	int term = item - node->data.sequence.items.start;
	yaml_node_t *term_node;

	term_node = yaml_document_get_node(&vlsp->vls_document, *item);
	if (parse_complex(vlsp, term_node, &vector[term][findex]) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_old_e_matrix: parse the "VNACAL 2.0" format e matrix
 *   @vlsp: pointer to vnacal_load info structure
 *   @e_matrices: array of El, Er, Em matrices, each with [cell][findex]
 *   @rows: number of expected rows
 *   @columns: number of expected columns
 *   @matrix_node: yaml node to parse
 */
static int parse_old_e_matrix(vnacal_load_state_t *vlsp,
	double complex ***e_matrices, int rows, int columns,
	yaml_node_t *matrix_node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    int findex = vlsp->vls_findex;
    int items;
    int cell = 0;

    if (matrix_node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected sequence",
		vcp->vc_filename, matrix_node->start_mark.line + 1);
	return -1;
    }
    items = matrix_node->data.sequence.items.top -
	matrix_node->data.sequence.items.start;
    if (items != rows) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected %d rows but found %d",
		vcp->vc_filename, matrix_node->start_mark.line + 1,
		rows, items);
	return -1;
    }
    for (yaml_node_item_t *row_item = matrix_node->data.sequence.items.start;
	 row_item < matrix_node->data.sequence.items.top; ++row_item) {
	yaml_node_t *row_node;
	yaml_node_item_t *column_item;

	row_node = yaml_document_get_node(&vlsp->vls_document, *row_item);
	if (row_node->type != YAML_SEQUENCE_NODE) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %ld) error: expected sequence",
		    vcp->vc_filename, row_node->start_mark.line + 1);
	    return -1;
	}
	items = row_node->data.sequence.items.top -
	    row_node->data.sequence.items.start;
	if (items != columns) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %ld) error: expected %d columns but found %d",
		    vcp->vc_filename, row_node->start_mark.line + 1,
		    columns, items);
	    return -1;
	}
	for (column_item = row_node->data.sequence.items.start;
	     column_item < row_node->data.sequence.items.top; ++column_item) {
	    yaml_node_t *cell_node;
	    yaml_node_item_t *term_item;

	    cell_node = yaml_document_get_node(&vlsp->vls_document,
		    *column_item);
	    if (cell_node->type != YAML_SEQUENCE_NODE) {
		_vnacal_error(vcp, VNAERR_SYNTAX,
			"%s (line %ld) error: expected sequence",
			vcp->vc_filename, row_node->start_mark.line + 1);
		return -1;
	    }
	    items = cell_node->data.sequence.items.top -
		cell_node->data.sequence.items.start;
	    if (items != 3) {
		_vnacal_error(vcp, VNAERR_SYNTAX,
			"%s (line %ld) error: expected 3 terms but found %d",
			vcp->vc_filename, row_node->start_mark.line + 1,
			items);
		return -1;
	    }
	    for (term_item = cell_node->data.sequence.items.start;
		    term_item < cell_node->data.sequence.items.top;
		    ++term_item) {
		int term = term_item - cell_node->data.sequence.items.start;
		yaml_node_t *term_node;

		term_node = yaml_document_get_node(&vlsp->vls_document,
			*term_item);
		if (term_node->type != YAML_SCALAR_NODE) {
		    _vnacal_error(vcp, VNAERR_SYNTAX,
			    "%s (line %ld) error: expected scalar",
			    vcp->vc_filename, row_node->start_mark.line + 1);
		    return -1;
		}
		if (parse_complex(vlsp, term_node,
			    &e_matrices[term][cell][findex]) == -1) {
		    return -1;
		}
	    }
	    ++cell;
	}
    }
    return 0;
}

/*
 * parse_matrix: parse a matrix
 *   @vlsp: pointer to vnacal_load info structure
 *   @matrix: where to place data
 *   @rows: number of expected rows
 *   @columns: number of expected columns
 *   @matrix_node: yaml node to parse
 *   @no_diagonal: matrix should have ~'s on its major diagonal
 */
static int parse_matrix(vnacal_load_state_t *vlsp, double complex **matrix,
	int rows, int columns, yaml_node_t *matrix_node, bool no_diagonal)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    int findex = vlsp->vls_findex;
    int items;
    int cell = 0;

    if (matrix_node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected sequence",
		vcp->vc_filename, matrix_node->start_mark.line + 1);
	return -1;
    }
    items = matrix_node->data.sequence.items.top -
	matrix_node->data.sequence.items.start;
    if (items != rows) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected %d rows but found %d",
		vcp->vc_filename, matrix_node->start_mark.line + 1,
		rows, items);
	return -1;
    }
    for (yaml_node_item_t *row_item = matrix_node->data.sequence.items.start;
	 row_item < matrix_node->data.sequence.items.top; ++row_item) {
	int row = row_item - matrix_node->data.sequence.items.start;
	yaml_node_t *row_node;
	yaml_node_item_t *column_item;

	row_node = yaml_document_get_node(&vlsp->vls_document, *row_item);
	if (row_node->type != YAML_SEQUENCE_NODE) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %ld) error: expected sequence",
		    vcp->vc_filename, row_node->start_mark.line + 1);
	    return -1;
	}
	items = row_node->data.sequence.items.top -
	    row_node->data.sequence.items.start;
	if (items != columns) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %ld) error: expected %d columns but found %d",
		    vcp->vc_filename, row_node->start_mark.line + 1,
		    columns, items);
	    return -1;
	}
	for (column_item = row_node->data.sequence.items.start;
	     column_item < row_node->data.sequence.items.top; ++column_item) {
	    int column = column_item - row_node->data.sequence.items.start;
	    yaml_node_t *cell_node;

	    cell_node = yaml_document_get_node(&vlsp->vls_document,
		    *column_item);
	    if (no_diagonal && row == column) {
		const char *value;

		if (cell_node->type != YAML_SCALAR_NODE) {
		    _vnacal_error(vcp, VNAERR_SYNTAX,
			    "%s (line %ld) error: expected scalar",
			    vcp->vc_filename, row_node->start_mark.line + 1);
		    return -1;
		}
		value = (const char *)cell_node->data.scalar.value;
		if ((value[0] != '~' || value[1] != '\000') &&
			strcmp(value, "null") != 0) {
		    _vnacal_error(vcp, VNAERR_SYNTAX,
			    "%s (line %ld) error: expected ~ in diagonal",
			    vcp->vc_filename, row_node->start_mark.line + 1);
		    return -1;
		}
	    } else {
		if (parse_complex(vlsp, cell_node,
			    &matrix[cell++][findex]) == -1) {
		    return -1;
		}
	    }
	}

    }
    return 0;
}

/*
 * parse_matrices: parse the error terms matrices
 *   @vlsp: pointer to vnacal_load info structure
 *   @vlp: pointer to vnacal_layout_t structure
 *   @calp: pointer to calibration structure
 *   @matrices: yaml_node pointers of matrices found
 */
static int parse_matrices(vnacal_load_state_t *vlsp, const vnacal_layout_t *vlp,
	vnacal_calibration_t *calp, yaml_node_t **matrices)
{
    double complex **e = calp->cal_error_term_vector;

    switch (calp->cal_type) {
    case VNACAL_T8:
    case VNACAL_TE10:
	{
	    const int ts_terms  = VL_TS_TERMS(vlp);
	    const int ts_offset = VL_TS_OFFSET(vlp);
	    const int ti_terms  = VL_TI_TERMS(vlp);
	    const int ti_offset = VL_TI_OFFSET(vlp);
	    const int tx_terms  = VL_TX_TERMS(vlp);
	    const int tx_offset = VL_TX_OFFSET(vlp);
	    const int tm_terms  = VL_TM_TERMS(vlp);
	    const int tm_offset = VL_TM_OFFSET(vlp);
	    double complex **ts = &e[ts_offset];
	    double complex **ti = &e[ti_offset];
	    double complex **tx = &e[tx_offset];
	    double complex **tm = &e[tm_offset];

	    if (parse_vector(vlsp, ts, ts_terms, matrices[TS]) == -1) {
		return -1;
	    }
	    if (parse_vector(vlsp, ti, ti_terms, matrices[TI]) == -1) {
		return -1;
	    }
	    if (parse_vector(vlsp, tx, tx_terms, matrices[TX]) == -1) {
		return -1;
	    }
	    if (parse_vector(vlsp, tm, tm_terms, matrices[TM]) == -1) {
		return -1;
	    }
	}
	if (calp->cal_type == VNACAL_TE10) {
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex **el = &e[el_offset];

	    if (parse_matrix(vlsp, el, el_rows, el_columns,
			matrices[EL], true) == -1) {
		return -1;
	    }
	}
	return 0;

    case VNACAL_U8:
    case VNACAL_UE10:
	{
	    const int um_terms  = VL_UM_TERMS(vlp);
	    const int um_offset = VL_UM_OFFSET(vlp);
	    const int ui_terms  = VL_UI_TERMS(vlp);
	    const int ui_offset = VL_UI_OFFSET(vlp);
	    const int ux_terms  = VL_UX_TERMS(vlp);
	    const int ux_offset = VL_UX_OFFSET(vlp);
	    const int us_terms  = VL_US_TERMS(vlp);
	    const int us_offset = VL_US_OFFSET(vlp);
	    double complex **um = &e[um_offset];
	    double complex **ui = &e[ui_offset];
	    double complex **ux = &e[ux_offset];
	    double complex **us = &e[us_offset];

	    if (parse_vector(vlsp, um, um_terms, matrices[UM]) == -1) {
		return -1;
	    }
	    if (parse_vector(vlsp, ui, ui_terms, matrices[UI]) == -1) {
		return -1;
	    }
	    if (parse_vector(vlsp, ux, ux_terms, matrices[UX]) == -1) {
		return -1;
	    }
	    if (parse_vector(vlsp, us, us_terms, matrices[US]) == -1) {
		return -1;
	    }
	}
	if (calp->cal_type == VNACAL_UE10) {
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex **el = &e[el_offset];

	    if (parse_matrix(vlsp, el, el_rows, el_columns,
			matrices[EL], true) == -1) {
		return -1;
	    }
	}
	return 0;

    case VNACAL_T16:
        {
            const int ts_rows    = VL_TS_ROWS(vlp);
            const int ts_columns = VL_TS_COLUMNS(vlp);
            const int ts_offset  = VL_TS_OFFSET(vlp);
            const int ti_rows    = VL_TI_ROWS(vlp);
            const int ti_columns = VL_TI_COLUMNS(vlp);
            const int ti_offset  = VL_TI_OFFSET(vlp);
            const int tx_rows    = VL_TX_ROWS(vlp);
            const int tx_columns = VL_TX_COLUMNS(vlp);
            const int tx_offset  = VL_TX_OFFSET(vlp);
            const int tm_rows    = VL_TM_ROWS(vlp);
            const int tm_columns = VL_TM_COLUMNS(vlp);
            const int tm_offset  = VL_TM_OFFSET(vlp);
            double complex **ts = &e[ts_offset];
            double complex **ti = &e[ti_offset];
            double complex **tx = &e[tx_offset];
            double complex **tm = &e[tm_offset];

	    if (parse_matrix(vlsp, ts, ts_rows, ts_columns,
			matrices[TS], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, ti, ti_rows, ti_columns,
			matrices[TI], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, tx, tx_rows, tx_columns,
			matrices[TX], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, tm, tm_rows, tm_columns,
			matrices[TM], false) == -1) {
		return -1;
	    }
	}
	return 0;

    case VNACAL_U16:
	{
	    const int um_rows    = VL_UM_ROWS(vlp);
	    const int um_columns = VL_UM_COLUMNS(vlp);
	    const int um_offset  = VL_UM_OFFSET(vlp);
	    const int ui_rows    = VL_UI_ROWS(vlp);
	    const int ui_columns = VL_UI_COLUMNS(vlp);
	    const int ui_offset  = VL_UI_OFFSET(vlp);
	    const int ux_rows    = VL_UX_ROWS(vlp);
	    const int ux_columns = VL_UX_COLUMNS(vlp);
	    const int ux_offset  = VL_UX_OFFSET(vlp);
	    const int us_rows    = VL_US_ROWS(vlp);
	    const int us_columns = VL_US_COLUMNS(vlp);
	    const int us_offset  = VL_US_OFFSET(vlp);
	    double complex **um = &e[um_offset];
	    double complex **ui = &e[ui_offset];
	    double complex **ux = &e[ux_offset];
	    double complex **us = &e[us_offset];

	    if (parse_matrix(vlsp, um, um_rows, um_columns,
			matrices[UM], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, ui, ui_rows, ui_columns,
			matrices[UI], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, ux, ux_rows, ux_columns,
			matrices[UX], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, us, us_rows, us_columns,
			matrices[US], false) == -1) {
		return -1;
	    }
	}
	return 0;

    case VNACAL_UE14:
	{
	    const int m_columns  = VL_M_COLUMNS(vlp);
	    const int um_terms   = VL_UM14_TERMS(vlp);
	    const int ui_terms   = VL_UI14_TERMS(vlp);
	    const int ux_terms   = VL_UX14_TERMS(vlp);
	    const int us_terms   = VL_US14_TERMS(vlp);
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex *packed_um[um_terms][m_columns];
	    double complex *packed_ui[ui_terms][m_columns];
	    double complex *packed_ux[ux_terms][m_columns];
	    double complex *packed_us[us_terms][m_columns];
	    double complex **el = &e[el_offset];

	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		const int um_offset = VL_UM14_OFFSET(vlp, m_column);
		const int ui_offset = VL_UI14_OFFSET(vlp, m_column);
		const int ux_offset = VL_UX14_OFFSET(vlp, m_column);
		const int us_offset = VL_US14_OFFSET(vlp, m_column);
		double complex **um = &e[um_offset];
		double complex **ui = &e[ui_offset];
		double complex **ux = &e[ux_offset];
		double complex **us = &e[us_offset];

		for (int term = 0; term < um_terms; ++term) {
		    packed_um[term][m_column] = um[term];
		}
		for (int term = 0; term < ui_terms; ++term) {
		    packed_ui[term][m_column] = ui[term];
		}
		for (int term = 0; term < ux_terms; ++term) {
		    packed_ux[term][m_column] = ux[term];
		}
		for (int term = 0; term < us_terms; ++term) {
		    packed_us[term][m_column] = us[term];
		}
	    }
	    if (parse_matrix(vlsp, &packed_um[0][0], um_terms, m_columns,
			matrices[UM], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, &packed_ui[0][0], ui_terms, m_columns,
			matrices[UI], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, &packed_ux[0][0], ux_terms, m_columns,
			matrices[UX], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, &packed_us[0][0], us_terms, m_columns,
			matrices[US], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, el, el_rows, el_columns,
			matrices[EL], true) == -1) {
		return -1;
	    }
	}
	return 0;

    case VNACAL_E12:
	{
	    const int m_rows    = VL_M_ROWS(vlp);
	    const int m_columns = VL_M_COLUMNS(vlp);
	    const int el_terms  = VL_EL12_TERMS(vlp);
	    const int er_terms  = VL_ER12_TERMS(vlp);
	    const int em_terms  = VL_EM12_TERMS(vlp);
	    double complex *packed_el[el_terms][m_columns];
	    double complex *packed_er[er_terms][m_columns];
	    double complex *packed_em[em_terms][m_columns];

	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		const int el_offset = VL_EL12_OFFSET(vlp, m_column);
		const int er_offset = VL_ER12_OFFSET(vlp, m_column);
		const int em_offset = VL_EM12_OFFSET(vlp, m_column);
		double complex **el = &e[el_offset];
		double complex **er = &e[er_offset];
		double complex **em = &e[em_offset];

		for (int term = 0; term < el_terms; ++term) {
		    packed_el[term][m_column] = el[term];
		}
		for (int term = 0; term < er_terms; ++term) {
		    packed_er[term][m_column] = er[term];
		}
		for (int term = 0; term < em_terms; ++term) {
		    packed_em[term][m_column] = em[term];
		}
	    }
	    if (vlsp->vls_major_version == 2) {
		double complex **e_matrices[3] = {
		    &packed_el[0][0], &packed_er[0][0], &packed_em[0][0]
		};
		assert(el_terms == m_rows);
		assert(er_terms == m_rows);
		assert(em_terms == m_rows);
		return parse_old_e_matrix(vlsp, e_matrices,
			m_rows, m_columns, matrices[E]);
	    }
	    if (parse_matrix(vlsp, &packed_el[0][0], el_terms, m_columns,
			matrices[EL], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, &packed_er[0][0], er_terms, m_columns,
			matrices[ER], false) == -1) {
		return -1;
	    }
	    if (parse_matrix(vlsp, &packed_em[0][0], em_terms, m_columns,
			matrices[EM], false) == -1) {
		return -1;
	    }
	}
	return 0;

    default:
	abort();
    }
}

/*
 * parse_data: parse the error term vector
 *   @vlsp: pointer to vnacal_load info structure
 *   @vlp: pointer to vnacal_layout_t structure
 *   @calp: pointer to calibration structure
 *   @node: error term matrix
 */
static int parse_data(vnacal_load_state_t *vlsp, const vnacal_layout_t *vlp,
	vnacal_calibration_t *calp, yaml_node_t *node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    int count;
    yaml_node_item_t *item;
    yaml_node_pair_t *pair;
    yaml_node_t *matrices[MATRIX_IDS];
    uint32_t required_matrices = 0;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected sequence for \"data\"",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    count = node->data.sequence.items.top - node->data.sequence.items.start;
    if (count != calp->cal_frequencies) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: "
		"expected %d frequency entries but found %d",
		vcp->vc_filename, node->start_mark.line + 1,
		calp->cal_frequencies, count);
	return -1;
    }
    for (item = node->data.sequence.items.start;
         item < node->data.sequence.items.top; ++item) {
	int findex = item - node->data.sequence.items.start;
	yaml_node_t *child = yaml_document_get_node(&vlsp->vls_document,
		*item);
	double frequency = -1.0;

	for (pair = child->data.mapping.pairs.start;
	     pair < child->data.mapping.pairs.top; ++pair) {
	    yaml_node_t *key, *value;
	    char prefix[4];
#define PREFIX(c1, c2, c3, c4) \
	((c1) | ((c2) << 8) | ((c3) << 16) | ((c4) << 24))

	    /*
	     * Get key and value nodes.
	     */
	    key   = yaml_document_get_node(&vlsp->vls_document, pair->key);
	    value = yaml_document_get_node(&vlsp->vls_document, pair->value);

	    /*
	     * Match keys, ignoring non-scalar keys.
	     */
	    if (key->type != YAML_SCALAR_NODE) {
		continue;
	    }

	    /*
	     * Get the first four bytes of the key, padded on the right
	     * with zeros, convert to integer in an endian neutral
	     * way and switch.  To match a string longer than three
	     * characters, simply switch on the first four, then use
	     * strcmp to check the rest.
	     */
	    (void)memset((void *)prefix, 0, sizeof(prefix));
	    (void)strncpy((void *)prefix, (void *)key->data.scalar.value,
		    sizeof(prefix));
	    switch (PREFIX(prefix[0], prefix[1], prefix[2], prefix[3])) {
	    case PREFIX('e',   0,   0,   0):
		matrices[E] = value;
		break;

	    case PREFIX('e', 'l',   0,   0):
		matrices[EL] = value;
		break;

	    case PREFIX('e', 'm',   0,   0):
		matrices[EM] = value;
		break;

	    case PREFIX('e', 'r',   0,   0):
		matrices[ER] = value;
		break;

	    case PREFIX('f',   0,   0,   0):
		if (key->data.scalar.value[1] == '\000') {
		    if (parse_double(vlsp, value, &frequency) == -1) {
			return -1;
		    }
		    break;
		}
		break;

	    case PREFIX('t', 'i',   0,   0):
		matrices[TI] = value;
		break;

	    case PREFIX('t', 'm',   0,   0):
		matrices[TM] = value;
		break;

	    case PREFIX('t', 's',   0,   0):
		matrices[TS] = value;
		break;

	    case PREFIX('t', 'x',   0,   0):
		matrices[TX] = value;
		break;

	    case PREFIX('u', 'i',   0,   0):
		matrices[UI] = value;
		break;

	    case PREFIX('u', 'm',   0,   0):
		matrices[UM] = value;
		break;

	    case PREFIX('u', 's',   0,   0):
		matrices[US] = value;
		break;

	    case PREFIX('u', 'x',   0,   0):
		matrices[UX] = value;
		break;

	    default:
		break;
	    }
	}
#undef PREFIX

	/*
	 * Make sure we have the required error terms matrices.
	 */
	if (vlsp->vls_major_version == 2) {
	    required_matrices |= 1 << E;

	} else switch (calp->cal_type) {
	case VNACAL_TE10:
	    required_matrices |= 1 << EL;
	    /*FALLTHROUGH*/

	case VNACAL_T8:
	case VNACAL_T16:
	    required_matrices |= T_MASK;
	    break;

	case VNACAL_UE10:
	case VNACAL_UE14:
	case _VNACAL_E12_UE14:
	    required_matrices |= 1 << EL;
	    /*FALLTHROUGH*/

	case VNACAL_U8:
	case VNACAL_U16:
	    required_matrices |= U_MASK;
	    break;

	case VNACAL_E12:
	    required_matrices |= E_MASK;
	    break;

	default:
	    abort();
	}
	for (int i = 0; i < MATRIX_IDS; ++i) {
	    if (((1 << i) & required_matrices) && matrices[i] == NULL) {
		_vnacal_error(vcp, VNAERR_SYNTAX,
			"%s (line %ld) error: missing required matrix \"%s\"",
			vcp->vc_filename, child->start_mark.line + 1,
			matrix_names[i]);
		return -1;
	    }
	}

	/*
	 * Make sure we have the frequency and that it's ascending.
	 */
	if (frequency < 0.0) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %ld) error: missing required field \"f\"",
		    vcp->vc_filename, child->start_mark.line + 1);
	    return -1;
	}
	if (findex > 1 &&
		frequency <= calp->cal_frequency_vector[findex - 1]) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %ld) error: frequencies are not in "
		    "ascending order",
		    vcp->vc_filename, child->start_mark.line + 1);
	    return -1;
	}

	/*
	 * Parse the error terms matrices.
	 */
	calp->cal_frequency_vector[findex] = frequency;
	vlsp->vls_findex = findex;
	if (parse_matrices(vlsp, vlp, calp, matrices) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_set: parse a calibration
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: the calibration
 */
static int parse_set(vnacal_load_state_t *vlsp, yaml_node_t *node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    yaml_node_pair_t *pair;
    const char *name = NULL;
    vnacal_type_t type = (vnacal_type_t)-1;
    int type_line = -1;
    int rows = -1, columns = -1, frequencies = -1;
    double complex z0 = VNADATA_DEFAULT_Z0;
    yaml_node_t *properties = NULL;
    yaml_node_t *data = NULL;
    vnacal_layout_t vl;
    vnacal_calibration_t *calp = NULL;

    if (node->type != YAML_MAPPING_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected mapping "
		"for \"calibration\"",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    for (pair = node->data.mapping.pairs.start;
         pair < node->data.mapping.pairs.top; ++pair) {
	yaml_node_t *key, *value;

	/*
	 * Get key and value nodes.
	 */
	key   = yaml_document_get_node(&vlsp->vls_document, pair->key);
	value = yaml_document_get_node(&vlsp->vls_document, pair->value);

	/*
	 * Ignore non-scalar keys.
	 */
	if (key->type != YAML_SCALAR_NODE) {
	    continue;
	}

	switch (key->data.scalar.value[0]) {
	case 'c':
	    if (strcmp((const char *)key->data.scalar.value, "columns") == 0) {
		if (parse_int(vlsp, value, &columns) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	case 'd':
	    if (strcmp((const char *)key->data.scalar.value, "data") == 0) {
		data = value;
		break;
	    }
	    break;

	case 'f':
	    if (strcmp((const char *)key->data.scalar.value,
			"frequencies") == 0) {
		if (parse_int(vlsp, value, &frequencies) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	case 'n':
	    if (strcmp((const char *)key->data.scalar.value, "name") == 0) {
		if (value->type != YAML_SCALAR_NODE) {
		    _vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: "
			    "expected scalar for name",
			    vcp->vc_filename, value->start_mark.line + 1);
		    return -1;
		}
		name = (const char *)value->data.scalar.value;
		break;
	    }
	    break;

	case 'p':
	    if (strcmp((const char *)key->data.scalar.value,
			"properties") == 0) {
		properties = value;
		break;
	    }
	    break;

	case 'r':
	    if (strcmp((const char *)key->data.scalar.value, "rows") == 0) {
		if (parse_int(vlsp, value, &rows) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	case 't':
	    if (strcmp((const char *)key->data.scalar.value, "type") == 0) {
		if (parse_type(vlsp, value, &type) == -1) {
		    return -1;
		}
		type_line = (int)value->start_mark.line + 1;
		break;
	    }
	    break;

	case 'z':
	    if (strcmp((const char *)key->data.scalar.value, "z0") == 0) {
		if (parse_complex(vlsp, value, &z0) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	default:
	    break;
	}
    }
    if (name == NULL || rows < 0 || columns < 0 || frequencies < 0 ||
	    data == NULL) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: missing required field \"name\", "
		"\"rows\", \"columns\", \"frequencies\" or \"data\"",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    if (vlsp->vls_major_version == 2) {
	if (type != (vnacal_type_t)-1 && type != VNACAL_E12) {
	    _vnacal_error(vcp, VNAERR_SYNTAX,
		    "%s (line %d) error: type unexpected in version %d.%d",
		    vcp->vc_filename, type_line,
		    vlsp->vls_major_version, vlsp->vls_minor_version);
	    return -1;
	}
	type = VNACAL_E12;

    } else if (type == (vnacal_type_t)-1) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: missing required field \"type\"",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    _vnacal_layout(&vl, type, rows, columns);
    if ((calp = _vnacal_calibration_alloc(vcp, type, rows, columns,
		    frequencies, VL_ERROR_TERMS(&vl))) == NULL) {
	return -1;
    }

    calp->cal_z0 = z0;
    if (properties != NULL) {
	calp->cal_properties = parse_properties(vlsp, properties);
	if (calp->cal_properties == NULL) {
	    _vnacal_calibration_free(calp);
	    return -1;
	}
    }
    if (parse_data(vlsp, &vl, calp, data) == -1) {
	_vnacal_calibration_free(calp);
	return -1;
    }
    if (_vnacal_add_calibration_common("vnacal_load", vcp, calp, name) == -1) {
	_vnacal_calibration_free(calp);
	return -1;
    }
    return 0;
}

/*
 * parse_calibrations: parse a sequence of calibrations
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: the calibration
 */
static int parse_calibrations(vnacal_load_state_t *vlsp, yaml_node_t *node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    yaml_node_item_t *item;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected sequence for \"calibrations\"",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    for (item = node->data.sequence.items.start;
         item < node->data.sequence.items.top; ++item) {
	yaml_node_t *child = yaml_document_get_node(&vlsp->vls_document,
		*item);

	if (parse_set(vlsp, child) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_document: parse the YAML document
 *   @vlsp: pointer to vnacal_load info structure
 *   @node: top-level YAML node
 */
static int parse_document(vnacal_load_state_t *vlsp, yaml_node_t *node)
{
    vnacal_t *vcp = vlsp->vls_vcp;
    yaml_node_pair_t *pair;

    if (node->type != YAML_MAPPING_NODE) {
	_vnacal_error(vcp, VNAERR_SYNTAX,
		"%s (line %ld) error: expected map at top level",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    for (pair = node->data.mapping.pairs.start;
         pair < node->data.mapping.pairs.top; ++pair) {
	yaml_node_t *key, *value;

	/*
	 * Get key and value.  Ignore non-scalar keys.
	 */
	key = yaml_document_get_node(&vlsp->vls_document, pair->key);
	value = yaml_document_get_node(&vlsp->vls_document, pair->value);
	if (key->type != YAML_SCALAR_NODE) {
	    continue;
	}

	/*
	 * Process global properties.
	 */
	if (strcmp((const char *)key->data.scalar.value, "properties") == 0) {
	    if ((vcp->vc_properties = parse_properties(vlsp, value)) == NULL) {
		return -1;
	    }
	}

	/*
	 * Process calibrations.
	 */
	if (strcmp((const char *)key->data.scalar.value, "calibrations") == 0 ||
		(vlsp->vls_major_version == 2 &&
		 strcmp((const char *)key->data.scalar.value, "sets") == 0)) {
	    if (parse_calibrations(vlsp, value) == -1) {
		return -1;
	    }
	}
    }
    return 0;
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
    vnacal_t *vcp;
    FILE *fp;
    vnacal_load_state_t vls;
    yaml_parser_t parser;
    yaml_node_t *root;
    bool delete_document = false;

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
	vnacal_free(vcp);
	return NULL;
    }

    /*
     * Init the parser state.
     */
    (void)memset((void *)&vls, 0, sizeof(vls));
    vls.vls_vcp = vcp;

    /*
     * Open and parse the file.
     */
    if ((fp = fopen(pathname, "r")) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "fopen: %s: %s",
		vcp->vc_filename, strerror(errno));
	vnacal_free(vcp);
	return NULL;
    }
    if (fscanf(fp, "#VNACAL %d.%d",
		&vls.vls_major_version, &vls.vls_minor_version) != 2) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line 1) error: "
		"expected #VNACAL <major>.<minor>",
		vcp->vc_filename);
	goto error;
    }
    if (vls.vls_major_version < 2 || vls.vls_major_version >= 4) {
	_vnacal_error(vcp, VNAERR_VERSION, "%s (line 1) error: "
		"unsupported version %d.%d",
		vcp->vc_filename, vls.vls_major_version,
		vls.vls_minor_version);
	goto error;
    }
    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, fp);
    if (!yaml_parser_load(&parser, &vls.vls_document)) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s (line %ld) error: %s",
		vcp->vc_filename, (long)parser.problem_mark.line + 1,
		parser.problem);
	goto error;
    }
    delete_document = true;
    (void)fclose(fp);
    fp = NULL;
    if ((root = yaml_document_get_root_node(&vls.vls_document)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYNTAX, "%s error: empty YAML document",
		vcp->vc_filename);
	goto error;
    }
    if (parse_document(&vls, root) == -1) {
	goto error;
    }
    yaml_document_delete(&vls.vls_document);
    return vcp;

error:
    if (delete_document) {
	yaml_document_delete(&vls.vls_document);
    }
    if (fp != NULL) {
	(void)fclose(fp);
    }
    vnacal_free(vcp);
    return NULL;
}
