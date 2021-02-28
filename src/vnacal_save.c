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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "vnacal_internal.h"


/*
 * add_integer: add an integer scalar to the yaml_document_t
 *   @document: yaml document
 *   @value: integer value
 */
static int add_integer(yaml_document_t *document, int value)
{
    char buf[2 + 3 * sizeof(int)];
    int tag;

    (void)sprintf(buf, "%d", value);
    if ((tag = yaml_document_add_scalar(document, NULL,
		    (yaml_char_t *)buf, strlen(buf),
		    YAML_ANY_SCALAR_STYLE)) == 0) {
	return -1;
    }
    return tag;
}

/*
 * add_double: add double scalar to the yaml_document_t
 *   @document: yaml document
 *   @value: real value
 *   @precision: number of significant figures to show
 */
static int add_double(yaml_document_t *document, double value, int precision)
{
    char buf[3 * sizeof(double) + 10];
    int tag;

    assert(precision >= 1);
    (void)sprintf(buf, "%.*e", precision - 1, value);
    if ((tag = yaml_document_add_scalar(document, NULL,
		    (yaml_char_t *)buf, strlen(buf),
		    YAML_ANY_SCALAR_STYLE)) == 0) {
	return -1;
    }
    return tag;
}

/*
 * add_complex: add a complex scalar to the yaml_document_t
 *   @document: yaml document
 *   @value: complex value
 *   @precision: number of significant figures to show
 */
static int add_complex(yaml_document_t *document, double complex value,
	int precision)
{
    double real = creal(value);
    double imag = cimag(value);
    char buf[3 * sizeof(double complex) + 20];
    int tag;

    assert(precision >= 1);
    if (precision == VNACAL_MAX_PRECISION) {
	(void)sprintf(buf, "%+a %+aj", real, imag);
    } else {
	(void)sprintf(buf, "%+.*e %+.*ej",
	    precision - 1, real,
	    precision - 1, imag);
    }
    if ((tag = yaml_document_add_scalar(document, NULL,
		    (yaml_char_t *)buf, strlen(buf),
		    YAML_ANY_SCALAR_STYLE)) == 0) {
	return -1;
    }
    return tag;
}

/*
 * add_mapping_entry: add a simple scalar tag to value mapping entry
 *   @vcp:      a pointer to the object returned by vnacal_create or vnacal_load
 *   @document: yaml document
 *   @t_map:    map into which we're adding
 *   @key:      string valued key
 *   @value:    tag of value to insert
 */
static int add_mapping_entry(vnacal_t *vcp, yaml_document_t *document,
	int t_map, const char *key, int t_value)
{
    int t_key;

    if ((t_key = yaml_document_add_scalar(document, NULL,
		    (yaml_char_t *)key, strlen(key),
		    YAML_ANY_SCALAR_STYLE)) == 0) {
	_vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    if (yaml_document_append_mapping_pair(document, t_map,
		t_key, t_value) == 0) {
	_vnacal_error(vcp, "yaml_document_add_mapping_pair: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    return 0;
}

/*
 * add_properties: add a property list to the document
 *   @vcp:      a pointer to the object returned by vnacal_create or vnacal_load
 *   @document: yaml document
 *   @root:     property list root
 */
static int add_properties(vnacal_t *vcp, yaml_document_t *document,
	vnaproperty_t *root)
{
    if (root == NULL) {
	return 0;
    }
    switch (vnaproperty_type(root)) {
    case VNAPROPERTY_SCALAR:
	{
	    int item;
	    const char *value;
	    yaml_scalar_style_t style = YAML_ANY_SCALAR_STYLE;

	    if ((value = vnaproperty_scalar_get(root)) == NULL) {
		_vnacal_error(vcp, "vnaproperty_scalar_get: %s: %s",
			vcp->vc_filename, strerror(errno));
		return -1;
	    }
	    if (strchr(value, '\n') != NULL) {
		style = YAML_LITERAL_SCALAR_STYLE;
	    }
	    if ((item = yaml_document_add_scalar(document, NULL,
		(yaml_char_t *)value, strlen(value), style)) == 0) {
		_vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
			vcp->vc_filename, strerror(errno));
		return -1;
	    }
	    return item;
	}

    case VNAPROPERTY_MAP:
	{
	    int map;
	    const vnaproperty_map_pair_t *vmprp;

	    if ((map = yaml_document_add_mapping(document, NULL,
		    YAML_ANY_MAPPING_STYLE)) == 0) {
		_vnacal_error(vcp, "yaml_document_add_mapping: %s: %s",
			vcp->vc_filename, strerror(errno));
		return -1;
	    }
	    for (vmprp = vnaproperty_map_begin(root); vmprp != NULL;
		 vmprp = vnaproperty_map_next(vmprp)) {
		int value;

		if ((value = add_properties(vcp, document,
				vmprp->vmpr_value)) == -1) {
		    return -1;
		}
		if (add_mapping_entry(vcp, document, map,
			    vmprp->vmpr_key, value) == -1) {
		    return -1;
		}
	    }
	    return map;
	}

    case VNAPROPERTY_LIST:
	{
	    int sequence;
	    int count = vnaproperty_list_count(root);

	    if ((sequence = yaml_document_add_sequence(document, NULL,
			    YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
		_vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
			vcp->vc_filename, strerror(errno));
		return -1;
	    }
	    for (int i = 0; i < count; ++i) {
		vnaproperty_t *property;
		int value;

		if ((property = vnaproperty_list_get(root, i)) == NULL) {
		    _vnacal_error(vcp, "vnaproperty_list_get: %s: %s",
			    vcp->vc_filename, strerror(errno));
		    return -1;
		}
		if ((value = add_properties(vcp, document, property)) == -1) {
		    return -1;
		}
		if (yaml_document_append_sequence_item(document, sequence,
			    value) == 0) {
		    _vnacal_error(vcp,
			    "yaml_document_append_sequence_item: "
			    "%s: %s", vcp->vc_filename, strerror(errno));
		    return -1;
		}
	    }
	    return sequence;
	}

    default:
	abort();
    }
}

/*
 * vnacal_save: create or overwrite a calibration file with new data
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @pathname: calibration file name
 *   @dotdir: directory under $HOME or NULL
 *
 *   The pathname and basename parameters work as in vnacal_load except
 *   that the $HOME/{pathname} directory is created if necessary.
 */
int vnacal_save(vnacal_t *vcp, const char *pathname, const char *dotdir)
{
    FILE *fp;
    yaml_document_t document;
    yaml_version_directive_t version = { 1, 1 };
    yaml_tag_directive_t tags[1];
    yaml_emitter_t emitter;
    bool delete_document = false;
    int t_root, t_properties, t_sets;

    if ((fp = _vnacal_open(vcp, pathname, dotdir, "w")) == NULL)
	return -1;

    if (!yaml_document_initialize(&document, &version, &tags[0], &tags[0],
		0, 0)) {
	_vnacal_error(vcp, "yaml_document_initialize: %s: %s",
		vcp->vc_filename, strerror(errno));
	goto error;
    }
    delete_document = true;
    t_root = yaml_document_add_mapping(&document, NULL,
	    YAML_BLOCK_MAPPING_STYLE);
    if (t_root == 0) {
	_vnacal_error(vcp, "yaml_document_add_mapping: %s: %s",
		vcp->vc_filename, strerror(errno));
	goto error;
    }
    if ((t_properties = add_properties(vcp, &document,
		    vcp->vc_properties)) == -1) {
	goto error;
    }
    if (t_properties != 0) {
	if (add_mapping_entry(vcp, &document, t_root, "properties",
		    t_properties) == -1) {
	    goto error;
	}
    }
    if ((t_sets = yaml_document_add_sequence(&document, NULL,
		    YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
	_vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
		vcp->vc_filename, strerror(errno));
	goto error;
    }
    if (add_mapping_entry(vcp, &document, t_root, "sets", t_sets) == -1) {
	goto error;
    }
    for (int set = 0; set < vcp->vc_sets; ++set) {
	vnacal_etermset_t *etsp = vcp->vc_set_vector[set];
	int t_set, t_value, t_matrix, t_row, t_data;

	t_set = yaml_document_add_mapping(&document, NULL,
		YAML_ANY_MAPPING_STYLE);
	if (t_set == 0) {
	    _vnacal_error(vcp, "yaml_document_add_mapping: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (yaml_document_append_sequence_item(&document, t_sets, t_set) == 0) {
	    _vnacal_error(vcp, "yaml_document_append_sequence_item: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if ((t_value = yaml_document_add_scalar(&document, NULL,
			(yaml_char_t *)etsp->ets_setname,
			strlen(etsp->ets_setname),
			YAML_ANY_SCALAR_STYLE)) == 0) {
	    _vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_set, "name", t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_integer(&document, etsp->ets_rows)) == -1) {
	    _vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_set, "rows", t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_integer(&document, etsp->ets_columns)) == -1) {
	    _vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_set, "columns",
		    t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_integer(&document, etsp->ets_frequencies)) == -1) {
	    _vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_set, "frequencies",
		    t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_complex(&document, etsp->ets_z0,
			vcp->vc_dprecision)) == -1) {
	    _vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_set, "z0", t_value) == -1) {
	    goto error;
	}
	if ((t_matrix = yaml_document_add_sequence(&document, NULL,
			YAML_ANY_SEQUENCE_STYLE)) == 0) {
	    _vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	{
	    int t_set_properties;

	    if ((t_set_properties = add_properties(vcp, &document,
			    etsp->ets_properties)) == -1) {
		goto error;
	    }
	    if (t_set_properties != 0) {
		if (add_mapping_entry(vcp, &document, t_set, "properties",
			    t_set_properties) == -1) {
		    goto error;
		}
	    }
	}
	if ((t_data = yaml_document_add_sequence(&document, NULL,
			YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
	    _vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_set, "data", t_data) == -1) {
	    goto error;
	}
	for (int findex = 0; findex < etsp->ets_frequencies; ++findex) {
	    int t_fitem;

	    if ((t_fitem = yaml_document_add_mapping(&document, NULL,
			    YAML_ANY_MAPPING_STYLE)) == 0) {
		_vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if (!yaml_document_append_sequence_item(&document, t_data,
			t_fitem)) {
		_vnacal_error(vcp, "yaml_document_append_sequence_item: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if ((t_value = add_double(&document,
			    etsp->ets_frequency_vector[findex],
			    vcp->vc_fprecision)) == -1) {
		_vnacal_error(vcp, "yaml_document_add_scalar: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if (add_mapping_entry(vcp, &document, t_fitem, "f",
			t_value) == -1) {
		goto error;
	    }
	    if ((t_matrix = yaml_document_add_sequence(&document, NULL,
			    YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
		_vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if (add_mapping_entry(vcp, &document, t_fitem, "e",
			t_matrix) == -1) {
		goto error;
	    }
	    for (int row = 0; row < etsp->ets_rows; ++row) {
		if ((t_row = yaml_document_add_sequence(&document, NULL,
				YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
		    _vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
			    vcp->vc_filename, strerror(errno));
		    goto error;
		}
		if (yaml_document_append_sequence_item(&document,
			    t_matrix, t_row) == 0) {
		    _vnacal_error(vcp,
			    "yaml_document_append_sequence_item: "
			    "%s: %s",
			    vcp->vc_filename, strerror(errno));
		    goto error;
		}
		for (int column = 0; column < etsp->ets_columns; ++column) {
		    int index = row * etsp->ets_columns + column;
		    vnacal_error_terms_t *etp =
			&etsp->ets_error_term_matrix[index];
		    int t_triple;

		    if ((t_triple = yaml_document_add_sequence(&document, NULL,
				    vcp->vc_dprecision <= 4 ?
					YAML_FLOW_SEQUENCE_STYLE :
					YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
			_vnacal_error(vcp, "yaml_document_add_sequence: %s: %s",
				vcp->vc_filename, strerror(errno));
			goto error;
		    }
		    if (yaml_document_append_sequence_item(&document,
				t_row, t_triple) == 0) {
			_vnacal_error(vcp,
				"yaml_document_append_sequence_item: "
				"%s: %s",
				vcp->vc_filename, strerror(errno));
			goto error;
		    }
		    for (int item = 0; item < 3; ++item) {
			if ((t_value = add_complex(&document,
					etp->et_data_vectors[item][findex],
					vcp->vc_dprecision)) == -1) {
			    _vnacal_error(vcp,
				    "yaml_document_add_scalar: %s: %s",
				    vcp->vc_filename, strerror(errno));
			    goto error;
			}
			if (yaml_document_append_sequence_item(&document,
				    t_triple, t_value) == 0) {
			    _vnacal_error(vcp,
				    "yaml_document_append_sequence_item: "
				    "%s: %s",
				    vcp->vc_filename, strerror(errno));
			    goto error;
			}
		    }
		}
	    }
	}
    }

    /*
     * Write the output file.
     */
    (void)fprintf(fp, "#VNACAL 2.0\n");
    if (!yaml_emitter_initialize(&emitter)) {
	_vnacal_error(vcp, "yaml_emitter_initialize: %s: error",
		vcp->vc_filename);
	goto error;
    }
    yaml_emitter_set_output_file(&emitter, fp);
    yaml_emitter_set_encoding(&emitter, YAML_UTF8_ENCODING);
    yaml_emitter_set_canonical(&emitter, 0);
    //yaml_emitter_set_indent(&emitter, 2);
    yaml_emitter_set_width(&emitter, 80);
    yaml_emitter_set_unicode(&emitter, 1);
    yaml_emitter_set_break(&emitter, YAML_ANY_BREAK);
    if (!yaml_emitter_open(&emitter)) {
	_vnacal_error(vcp, "yaml_emitter_open: %s: error: %s",
		vcp->vc_filename, emitter.problem);
	goto error;
    }
    if (!yaml_emitter_dump(&emitter, &document)) {
	_vnacal_error(vcp, "yaml_emitter_dump: %s: error: %s",
		vcp->vc_filename, emitter.problem);
	goto error;
    }
    if (!yaml_emitter_close(&emitter)) {
	_vnacal_error(vcp, "yaml_emitter_close: %s: error: %s",
		vcp->vc_filename, emitter.problem);
	goto error;
    }
    (void)yaml_emitter_delete(&emitter);
    //(void)yaml_document_delete(&document);
    delete_document = false;

    if (fclose(fp) == -1) {
	_vnacal_error(vcp, "fclose: %s: %s", vcp->vc_filename,
	    strerror(errno));
	goto error;
    }
    return 0;

error:
    if (delete_document) {
	(void)yaml_document_delete(&document);
    }
    if (fp != NULL) {
	(void)fclose(fp);
    }
    return -1;
}
