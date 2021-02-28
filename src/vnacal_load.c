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
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "vnacal_internal.h"


/*
 * vnacal_load_t: parser state
 */
typedef struct vnacal_load {
    vnacal_t *vl_vcp;			/* back pointer to vnacal_t structure */
    yaml_document_t vl_document;	/* yaml document */
} vnacal_load_t;

/*
 * parse_int: parse an integer
 *   @vlp: parser state
 *   @result: returned integer
 */
static int parse_int(vnacal_load_t *vlp, yaml_node_t *node, int *result)
{
    vnacal_t *vcp = vlp->vl_vcp;
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
    _vnacal_error(vcp, "%s (line %ld) error: expected integer",
	    vcp->vc_filename, node->start_mark.line + 1);
    errno = EBADMSG;
    return -1;
}

/*
 * parse_double: parse a double
 *   @vlp: parser state
 *   @result: returned double
 */
static int parse_double(vnacal_load_t *vlp, yaml_node_t *node, double *result)
{
    vnacal_t *vcp = vlp->vl_vcp;
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
    _vnacal_error(vcp, "%s (line %ld) error: expected number",
	    vcp->vc_filename, node->start_mark.line + 1);
    errno = EBADMSG;
    return -1;
}

/*
 * parse_complex: parse a complex number
 *   @vlp: parser state
 *   @result: returned double
 */
static int parse_complex(vnacal_load_t *vlp, yaml_node_t *node,
	double complex *result)
{
    vnacal_t *vcp = vlp->vl_vcp;
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
    if (*cur == '-') {
	code += 6;
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
	code += 3;
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
    case 3:	/* j */
	*result = I;
	break;
    case 4:	/* number j */
	*result = value1 * I;
	break;
    case 5:	/* number number j */
	*result = value1 + value2 * I;
	break;
    case 9:	/* -j */
	*result = -I;
	break;
    case 10:	/* number - j */
	*result = value1 - I;
	break;
    default:
	goto error;
    }
    return 0;

error:
    _vnacal_error(vcp, "%s (line %ld) error: expected number",
	    vcp->vc_filename, node->start_mark.line + 1);
    errno = EBADMSG;
    return -1;
}

/*
 * parse_properties: parse user properties
 *   @vlp: parser state
 *   @node: user properties
 */
static vnaproperty_t *parse_properties(vnacal_load_t *vlp, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    vnaproperty_t *root = NULL;
    vnaproperty_t *subtree = NULL;

    switch (node->type) {
    case YAML_SCALAR_NODE:
	if ((root = vnaproperty_scalar_alloc((const char *)node->data.
			scalar.value)) == NULL) {
	    _vnacal_error(vcp, "%s (line %ld) error: "
		    "vnaproperty_scalar_alloc: %s",
		    vcp->vc_filename, (long)node->start_mark.line + 1,
		    strerror(errno));
	    goto out;
	}
	return root;

    case YAML_MAPPING_NODE:
	{
	    yaml_node_pair_t *pair;

	    if ((root = vnaproperty_map_alloc()) == NULL) {
		_vnacal_error(vcp, "%s (line %ld) error: "
			"vnaproperty_map_alloc: %s",
			vcp->vc_filename, (long)node->start_mark.line + 1,
			strerror(errno));
		goto out;
	    }
	    for (pair = node->data.mapping.pairs.start;
		 pair < node->data.mapping.pairs.top; ++pair) {
		yaml_node_t *key, *value;

		key = yaml_document_get_node(&vlp->vl_document, pair->key);
		if (key->type != YAML_SCALAR_NODE) {
		    _vnacal_error(vcp, "%s (line %ld) warning: non-scalar "
			    "property key ignored\n",
			    vcp->vc_filename, key->start_mark.line + 1);
		    continue;
		}
		value = yaml_document_get_node(&vlp->vl_document, pair->value);
		subtree = parse_properties(vlp, value);
		if (subtree == NULL) {
		    goto out;
		}
		if (vnaproperty_map_set(root,
			    (const char *)key->data.scalar.value,
			    subtree) == -1) {
		    _vnacal_error(vcp, "%s (line %ld) error: "
			    "vnaproperty_map_set: %s",
			    vcp->vc_filename, (long)key->start_mark.line + 1,
			    strerror(errno));
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
		_vnacal_error(vcp, "%s (line %ld) error: "
			"vnaproperty_list_alloc: %s",
			vcp->vc_filename, (long)node->start_mark.line + 1,
			strerror(errno));
		goto out;
	    }
	    for (item = node->data.sequence.items.start;
		 item < node->data.sequence.items.top; ++item) {
		int term = item - node->data.sequence.items.start;
		yaml_node_t *value;

		value = yaml_document_get_node(&vlp->vl_document, *item);
		subtree = parse_properties(vlp, value);
		if (subtree == NULL) {
		    goto out;
		}
		if (vnaproperty_list_set(root, term, subtree) == -1) {
		    _vnacal_error(vcp, "%s (line %ld) error: "
			    "vnaproperty_list_set: %s",
			    vcp->vc_filename, (long)node->start_mark.line + 1,
			    strerror(errno));
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
 * parse_matrix_column: parse a row of the error term matrix
 *   @vlp: parser state
 *   @etsp: error terms set
 *   @findex: frequency index
 *   @node: error term matrix
 */
static int parse_matrix_column(vnacal_load_t *vlp, vnacal_etermset_t *etsp,
	int row, int column, int findex, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    int cell = row * etsp->ets_columns + column;
    int terms;
    yaml_node_item_t *item;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected sequence "
		"value", vcp->vc_filename,
		node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    terms = node->data.sequence.items.top - node->data.sequence.items.start;
    if (terms != 3) {
	_vnacal_error(vcp, "%s (line %ld) error: expected 3 terms "
		"but found %d",
		vcp->vc_filename, node->start_mark.line + 1,
		terms);
	errno = EBADMSG;
	return -1;
    }
    for (item = node->data.sequence.items.start;
	 item < node->data.sequence.items.top; ++item) {
	int term = item - node->data.sequence.items.start;
	yaml_node_t *term_node;
	vnacal_error_terms_t *etp = &etsp->ets_error_term_matrix[cell];

	term_node = yaml_document_get_node(&vlp->vl_document, *item);
	if (parse_complex(vlp, term_node,
		    &etp->et_data_vectors[term][findex]) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_matrix_row: parse a row of the error term matrix
 *   @vlp: parser state
 *   @etsp: error terms set
 *   @row: row of matrix
 *   @findex: frequency index
 *   @node: error term matrix
 */
static int parse_matrix_row(vnacal_load_t *vlp, vnacal_etermset_t *etsp,
	int row, int findex, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    int columns;
    yaml_node_item_t *item;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected sequence "
		"value", vcp->vc_filename, node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    columns = node->data.sequence.items.top - node->data.sequence.items.start;
    if (columns != etsp->ets_columns) {
	_vnacal_error(vcp, "%s (line %ld) error: expected %d columns "
		"but found %d",
		vcp->vc_filename, node->start_mark.line + 1,
		etsp->ets_columns, columns);
	errno = EBADMSG;
	return -1;
    }
    for (item = node->data.sequence.items.start;
	 item < node->data.sequence.items.top; ++item) {
	int column = item - node->data.sequence.items.start;
	yaml_node_t *column_node;

	column_node = yaml_document_get_node(&vlp->vl_document, *item);
	if (parse_matrix_column(vlp, etsp, row, column,
		    findex, column_node) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_matrix: parse the error term matrix
 *   @vlp: parser state
 *   @etsp: error terms set
 *   @findex: frequency index
 *   @node: error term matrix
 */
static int parse_matrix(vnacal_load_t *vlp, vnacal_etermset_t *etsp,
	int findex, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    int rows;
    yaml_node_item_t *item;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected sequence "
		"value", vcp->vc_filename, node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    rows = node->data.sequence.items.top - node->data.sequence.items.start;
    if (rows != etsp->ets_rows) {
	_vnacal_error(vcp, "%s (line %ld) error: expected %d rows "
		"but found %d",
		vcp->vc_filename, node->start_mark.line + 1,
		etsp->ets_rows, rows);
	errno = EBADMSG;
	return -1;
    }
    for (item = node->data.sequence.items.start;
         item < node->data.sequence.items.top; ++item) {
	int row = item - node->data.sequence.items.start;
	yaml_node_t *row_node;

	row_node = yaml_document_get_node(&vlp->vl_document, *item);
	if (parse_matrix_row(vlp, etsp, row, findex, row_node) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_data: parse the error term vector
 *   @vlp: parser state
 *   @etsp: error terms set
 *   @node: error term matrix
 */
static int parse_data(vnacal_load_t *vlp, vnacal_etermset_t *etsp,
	yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    int count;
    yaml_node_item_t *item;
    yaml_node_pair_t *pair;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected sequence "
		"value for \"data\"",
		vcp->vc_filename, node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    count = node->data.sequence.items.top - node->data.sequence.items.start;
    if (count != etsp->ets_frequencies) {
	_vnacal_error(vcp, "%s (line %ld) error: expected %d frequency entries "
		"but found %d",
		vcp->vc_filename, node->start_mark.line + 1,
		etsp->ets_frequencies, count);
	errno = EBADMSG;
	return -1;
    }
    for (item = node->data.sequence.items.start;
         item < node->data.sequence.items.top; ++item) {
	int findex = item - node->data.sequence.items.start;
	yaml_node_t *child = yaml_document_get_node(&vlp->vl_document, *item);
	double frequency = -1.0;
	yaml_node_t *matrix = NULL;

	for (pair = child->data.mapping.pairs.start;
	     pair < child->data.mapping.pairs.top; ++pair) {
	    yaml_node_t *key, *value;

	    /*
	     * Get key and value nodes.
	     */
	    key   = yaml_document_get_node(&vlp->vl_document, pair->key);
	    value = yaml_document_get_node(&vlp->vl_document, pair->value);

	    /*
	     * Ignore non-scalar keys.
	     */
	    if (key->type != YAML_SCALAR_NODE) {
		continue;
	    }

	    switch (key->data.scalar.value[0]) {
	    case 'e':
		if (strcmp((const char *)key->data.scalar.value, "e") == 0) {
		    matrix = value;
		    break;
		}
		break;

	    case 'f':
		if (strcmp((const char *)key->data.scalar.value, "f") == 0) {
		    if (parse_double(vlp, value, &frequency) == -1) {
			return -1;
		    }
		    break;
		}
		break;
	    }
	}
	if (frequency < 0.0 || matrix == NULL) {
	    _vnacal_error(vcp,
		    "%s (line %ld) error: missing required field "
		    "\"f\" or \"e\"",
		    vcp->vc_filename, child->start_mark.line + 1);
	    return -1;
	}
	if (findex > 1 && frequency <= etsp->ets_frequency_vector[findex - 1]) {
	    _vnacal_error(vcp,
		    "%s (line %ld) error: frequencies are not in "
		    "ascending order",
		    vcp->vc_filename, child->start_mark.line + 1);
	    return -1;
	}
	etsp->ets_frequency_vector[findex] = frequency;
	if (parse_matrix(vlp, etsp, findex, matrix) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_set: parse a set
 *   @vlp: parser state
 *)   @node: the set
 *   @set: set number
 */
static int parse_set(vnacal_load_t *vlp, int set, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    yaml_node_pair_t *pair;
    const char *name = NULL;
    int rows = -1, columns = -1, frequencies = -1;
    double complex z0 = VNADATA_DEFAULT_Z0;
    yaml_node_t *properties = NULL;
    yaml_node_t *data = NULL;
    vnacal_etermset_t *etsp = NULL;

    if (node->type != YAML_MAPPING_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected mapping "
		"value for \"set\"",
		vcp->vc_filename, node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    for (pair = node->data.mapping.pairs.start;
         pair < node->data.mapping.pairs.top; ++pair) {
	yaml_node_t *key, *value;

	/*
	 * Get key and value nodes.
	 */
	key   = yaml_document_get_node(&vlp->vl_document, pair->key);
	value = yaml_document_get_node(&vlp->vl_document, pair->value);

	/*
	 * Ignore non-scalar keys.
	 */
	if (key->type != YAML_SCALAR_NODE) {
	    continue;
	}

	switch (key->data.scalar.value[0]) {
	case 'c':
	    if (strcmp((const char *)key->data.scalar.value, "columns") == 0) {
		if (parse_int(vlp, value, &columns) == -1) {
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
		if (parse_int(vlp, value, &frequencies) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	case 'n':
	    if (strcmp((const char *)key->data.scalar.value, "name") == 0) {
		if (value->type != YAML_SCALAR_NODE) {
		    _vnacal_error(vcp, "%s (line %ld) error: expected scalar "
			    "value for name",
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
		if (parse_int(vlp, value, &rows) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	case 'z':
	    if (strcmp((const char *)key->data.scalar.value, "z0") == 0) {
		if (parse_complex(vlp, value, &z0) == -1) {
		    return -1;
		}
		break;
	    }
	    break;

	default:
	    break;
	}
    }
    if (name == NULL || rows < 0 || columns < 0 || frequencies < 0) {
	_vnacal_error(vcp,
		"%s (line %ld) error: missing required field \"name\", "
		"\"rows\", \"columns\" or \"frequencies\"",
		vcp->vc_filename, node->start_mark.line + 1);
	return -1;
    }
    if ((vcp->vc_set_vector[set] = etsp = _vnacal_etermset_alloc(vcp, name,
		    rows, columns, frequencies)) == NULL) {
	return -1;
    }
    etsp->ets_z0 = z0;
    if (properties != NULL) {
	etsp->ets_properties = parse_properties(vlp, properties);
	if (etsp->ets_properties == NULL) {
	    return -1;
	}
    }
    if (data != NULL) {
	if (parse_data(vlp, etsp, data) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * parse_sets: parse a sequence of sets
 *   @vlp: parser state
 *   @node: the set
 */
static int parse_sets(vnacal_load_t *vlp, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    yaml_node_item_t *item;
    int sets;

    if (node->type != YAML_SEQUENCE_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected sequence "
		"value for \"sets\"",
		vcp->vc_filename, node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    sets = node->data.sequence.items.top - node->data.sequence.items.start;
    if ((vcp->vc_set_vector = calloc(sets,
		    sizeof(vnacal_etermset_t *))) == NULL) {
	_vnacal_error(vcp, "%s (line %ld) error: "
		"calloc: %s",
		vcp->vc_filename, (long)node->start_mark.line + 1,
		strerror(errno));
	return -1;
    }
    vcp->vc_sets = sets;
    for (item = node->data.sequence.items.start;
         item < node->data.sequence.items.top; ++item) {
	yaml_node_t *child = yaml_document_get_node(&vlp->vl_document, *item);

	parse_set(vlp, item - node->data.sequence.items.start, child);
    }
    return 0;
}

/*
 * parse_document: parse the YAML document
 *   @vlp: parser state
 *   @node: top-level YAML node
 */
static int parse_document(vnacal_load_t *vlp, yaml_node_t *node)
{
    vnacal_t *vcp = vlp->vl_vcp;
    yaml_node_pair_t *pair;

    if (node->type != YAML_MAPPING_NODE) {
	_vnacal_error(vcp, "%s (line %ld) error: expected map at top level",
		vcp->vc_filename, node->start_mark.line + 1);
	errno = EBADMSG;
	return -1;
    }
    for (pair = node->data.mapping.pairs.start;
         pair < node->data.mapping.pairs.top; ++pair) {
	yaml_node_t *key, *value;

	/*
	 * Get key and value.  Ignore non-scalar keys.
	 */
	key = yaml_document_get_node(&vlp->vl_document, pair->key);
	value = yaml_document_get_node(&vlp->vl_document, pair->value);
	if (key->type != YAML_SCALAR_NODE) {
	    continue;
	}

	/*
	 * Process global properties.
	 */
	if (strcmp((const char *)key->data.scalar.value, "properties") == 0) {
	    if ((vcp->vc_properties = parse_properties(vlp, value)) == NULL) {
		return -1;
	    }
	}

	/*
	 * Process sets.
	 */
	if (strcmp((const char *)key->data.scalar.value, "sets") == 0) {
	    if (parse_sets(vlp, value) == -1) {
		return -1;
	    }
	}
    }
    return 0;
}

/*
 * vnacal_load: load the calibration from a file
 *   @pathname: calibration file name
 *   @dotdir: directory under $HOME or NULL
 *   @error_fn: error reporting callback or NULL
 *   @error_arg: arbitrary argument passed through to error_fn or NULL
 *
 *   If error_fn is non-NULL, then vnacal_load and subsequent functions report
 *   error messages using error_fn before returning failure to the caller.
 */
vnacal_t *vnacal_load(const char *pathname, const char *dotdir,
	vnacal_error_fn_t *error_fn, void *error_arg)
{
    vnacal_t *vcp;
    FILE *fp;
    vnacal_load_t vl;
    yaml_parser_t parser;
    yaml_node_t *root;
    double vnacal_version;
    bool delete_document = false;

    /*
     * Allocate the vnaset_t structure and initialize the error reporting function.
     */
    if ((vcp = (vnacal_t *)malloc(sizeof(vnacal_t))) == NULL) {
	if (error_fn != NULL) {
	    char buf[80];

	    (void)snprintf(buf, sizeof(buf), "vnacal_create: %s",
		    strerror(errno));
	    buf[sizeof(buf)-1] = '\000';
	    (*error_fn)(error_arg, buf);
	}
	return NULL;
    }
    (void)memset((void *)vcp, 0, sizeof(vnacal_t));
    vcp->vc_fprecision = VNACAL_DEFAULT_DATA_PRECISION;
    vcp->vc_dprecision = VNACAL_DEFAULT_FREQUENCY_PRECISION;
    vcp->vc_error_fn   = error_fn;
    vcp->vc_error_arg  = error_arg;

    /*
     * Init the parser state.
     */
    (void)memset((void *)&vl, 0, sizeof(vl));
    vl.vl_vcp = vcp;

    /*
     * Open and parse the file.
     */
    if ((fp = _vnacal_open(vcp, pathname, dotdir, "r")) == NULL) {
	vnacal_free(vcp);
	return NULL;
    }
    if (fscanf(fp, "#VNACAL %lf", &vnacal_version) != 1) {
	_vnacal_error(vcp, "%s (line 1) error: expected #VNACAL",
		vcp->vc_filename);
	errno = EBADMSG;
	goto error;
    }
    if (vnacal_version < 2.0 || vnacal_version >= 3.0) {
	_vnacal_error(vcp, "%s (line 1) error: unsupported version %.1f",
		vcp->vc_filename, vnacal_version);
	errno = ENOPROTOOPT;
	goto error;
    }
    yaml_parser_initialize(&parser);
    yaml_parser_set_input_file(&parser, fp);
    if (!yaml_parser_load(&parser, &vl.vl_document)) {
	_vnacal_error(vcp, "%s (line %ld) error: %s",
		vcp->vc_filename, (long)parser.problem_mark.line + 1,
		parser.problem);
	errno = EBADMSG;
	goto error;
    }
    delete_document = true;
    (void)fclose(fp);
    fp = NULL;
    if ((root = yaml_document_get_root_node(&vl.vl_document)) == NULL) {
	_vnacal_error(vcp, "%s error: empty YAML document",
		vcp->vc_filename);
	errno = EBADMSG;
	goto error;
    }
    if (parse_document(&vl, root) == -1) {
	goto error;
    }
    yaml_document_delete(&vl.vl_document);
    return vcp;

error:
    if (delete_document) {
	yaml_document_delete(&vl.vl_document);
    }
    if (fp != NULL) {
	(void)fclose(fp);
    }
    vnacal_free(vcp);
    return NULL;
}
