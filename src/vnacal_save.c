/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "vnacal_internal.h"
#include "vnaproperty_internal.h"


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
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @document: yaml document
 *   @t_map:    map into which we're adding
 *   @key:      string valued key
 *   @value:    tag of value to insert
 */
static int add_mapping_entry(vnacal_t *vcp, yaml_document_t *document,
	int t_map, const char *key, int t_value)
{
    int t_key;

    errno = 0;
    if ((t_key = yaml_document_add_scalar(document, NULL,
		    (yaml_char_t *)key, strlen(key),
		    YAML_ANY_SCALAR_STYLE)) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_add_scalar: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    if (yaml_document_append_mapping_pair(document, t_map,
		t_key, t_value) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_add_mapping_pair: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    return 0;
}

/*
 * error_terms_common_arguments_t: common error term matrix arguments
 */
typedef struct error_terms_common_arguments {
    vnacal_t	       *etca_vcp;
    yaml_document_t    *etca_document;
    int			etca_mapping;
    int			etca_findex;
} error_terms_common_arguments_t;

/*
 * add_vector: add a complex vector to an existing map
 *   @etca:     arguments passed through from vnacal_save
 *   @name:     name of vector
 *   @vector:   data
 *   @length:   length of vector
 */
static int add_vector(error_terms_common_arguments_t etca,
	const char *name, double complex *const *vector, int length)
{
    vnacal_t *vcp = etca.etca_vcp;
    yaml_document_t *document = etca.etca_document;
    int mapping = etca.etca_mapping;
    int findex = etca.etca_findex;
    int t_vector, t_value;

    yaml_sequence_style_t sequence_style =
	(2 * length * (vcp->vc_dprecision + 8) <= 70) ?
	YAML_FLOW_SEQUENCE_STYLE : YAML_BLOCK_SEQUENCE_STYLE;

    if ((t_vector = yaml_document_add_sequence(document, NULL,
		    sequence_style)) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_add_sequence: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    if (add_mapping_entry(vcp, document, mapping, name, t_vector) == -1) {
	return -1;
    }
    for (int i = 0; i < length; ++i) {
	if ((t_value = add_complex(document, vector[i][findex],
			vcp->vc_dprecision)) == -1) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    return -1;
	}
	if (yaml_document_append_sequence_item(document,
		    t_vector, t_value) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_append_sequence_item: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    return -1;
	}
    }
    return 0;
}

/*
 * add_matrix: add a complex matrix to an existing map
 *   @etca:     arguments passed through from vnacal_save
 *   @name:        name of vector
 *   @matrix:      data
 *   @rows:        rows in matrix
 *   @columns:     columns in matrix
 *   @no_diagonal: matrix is missing the major diagonal
 */
static int add_matrix(error_terms_common_arguments_t etca,
	const char *name, double complex *const *matrix,
	int rows, int columns, bool no_diagonal)
{
    vnacal_t *vcp = etca.etca_vcp;
    yaml_document_t *document = etca.etca_document;
    int mapping = etca.etca_mapping;
    int findex = etca.etca_findex;
    int t_matrix, t_value;
    yaml_sequence_style_t row_style;

    row_style = (2 * columns * (vcp->vc_dprecision + 8) <= 70) ?
	YAML_FLOW_SEQUENCE_STYLE : YAML_BLOCK_SEQUENCE_STYLE;
    if ((t_matrix = yaml_document_add_sequence(document, NULL,
		    YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_add_sequence: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    if (add_mapping_entry(vcp, document, mapping, name, t_matrix) == -1) {
	return -1;
    }
    for (int row = 0; row < rows; ++row) {
	int t_row;

	if ((t_row = yaml_document_add_sequence(document, NULL,
			row_style)) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_sequence: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    return -1;
	}
	if (yaml_document_append_sequence_item(document,
		    t_matrix, t_row) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_append_sequence_item: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    return -1;
	}
	for (int column = 0; column < columns; ++column) {
	    if (no_diagonal && row == column) {
		if ((t_value = yaml_document_add_scalar(document, NULL,
				(yaml_char_t *)"~", 1,
				YAML_ANY_SCALAR_STYLE)) == 0) {
		    if (errno == 0) {
			errno = EINVAL;
		    }
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "yaml_document_add_scalar: %s: %s",
			    vcp->vc_filename, strerror(errno));
		    return -1;
		}
	    } else {
		if ((t_value = add_complex(document, (*matrix++)[findex],
				vcp->vc_dprecision)) == -1) {
		    if (errno == 0) {
			errno = EINVAL;
		    }
		    _vnacal_error(vcp, VNAERR_SYSTEM,
			    "yaml_document_add_scalar: %s: %s",
			    vcp->vc_filename, strerror(errno));
		    return -1;
		}
	    }
	    if (yaml_document_append_sequence_item(document,
			t_row, t_value) == 0) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"yaml_document_append_sequence_item: %s: %s",
			vcp->vc_filename, strerror(errno));
		return -1;
	    }
	}
    }
    return 0;
}

/*
 * add_error_parameters: add error parameters matrices for one frequency
 *   @etca:     arguments passed through from vnacal_save
 *   @calp:     pointer to calibration structure
 *   @vlp:      pointer to vnacal_layout_t structure
 */
static int add_error_parameters(error_terms_common_arguments_t etca,
	const vnacal_calibration_t *calp, const vnacal_layout_t *vlp)
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
	    if (add_vector(etca, "ts", ts, ts_terms) == -1) {
		return -1;
	    }

	    if (add_vector(etca, "ti", ti, ti_terms) == -1) {
		return -1;
	    }

	    if (add_vector(etca, "tx", tx, tx_terms) == -1) {
		return -1;
	    }

	    if (add_vector(etca, "tm", tm, tm_terms) == -1) {
		return -1;
	    }
	}
	if (calp->cal_type == VNACAL_TE10) {
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex **el = &e[el_offset];

	    if (add_matrix(etca, "el",
			el, el_rows, el_columns, true) == -1) {
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
	    if (add_vector(etca, "um", um, um_terms) == -1) {
		return -1;
	    }

	    if (add_vector(etca, "ui", ui, ui_terms) == -1) {
		return -1;
	    }

	    if (add_vector(etca, "ux", ux, ux_terms) == -1) {
		return -1;
	    }

	    if (add_vector(etca, "us", us, us_terms) == -1) {
		return -1;
	    }
	}
	if (calp->cal_type == VNACAL_UE10) {
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex **el = &e[el_offset];

	    if (add_matrix(etca, "el",
			el, el_rows, el_columns, true) == -1) {
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
	    if (add_matrix(etca, "ts", ts,
			ts_rows, ts_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "ti", ti,
			ti_rows, ti_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "tx", tx,
			tx_rows, tx_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "tm", tm,
			tm_rows, tm_columns, false) == -1) {
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
	    if (add_matrix(etca, "um", um,
			um_rows, um_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "ui", ui,
			ui_rows, ui_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "ux", ux,
			ux_rows, ux_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "us", us,
			us_rows, us_columns, false) == -1) {
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
	    if (add_matrix(etca, "um", &packed_um[0][0],
			um_terms, m_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "ui", &packed_ui[0][0],
			ui_terms, m_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "ux", &packed_ux[0][0],
			ux_terms, m_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "us", &packed_us[0][0],
			us_terms, m_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "el",
			el, el_rows, el_columns, true) == -1) {
		return -1;
	    }
	}
	return 0;

    case VNACAL_E12:
	{
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
	    if (add_matrix(etca, "el", &packed_el[0][0],
			el_terms, m_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "er", &packed_er[0][0],
			er_terms, m_columns, false) == -1) {
		return -1;
	    }
	    if (add_matrix(etca, "em", &packed_em[0][0],
			em_terms, m_columns, false) == -1) {
		return -1;
	    }
	}
	return 0;

    default:
	abort();
    }
    return -1;
}

/*
 * add_properties: add a property list to the document
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @document: yaml document
 *   @root:     property list root
 */
static int add_properties(vnacal_t *vcp, yaml_document_t *document,
	const vnaproperty_t *root)
{
    vnaproperty_yaml_t vyml;

    (void)memset((void *)&vyml, 0, sizeof(vyml));
    vyml.vyml_document = (void *)document;
    vyml.vyml_filename = vcp->vc_filename;
    vyml.vyml_error_fn = vcp->vc_error_fn;
    vyml.vyml_error_arg = vcp->vc_error_arg;

    return _vnaproperty_yaml_export(&vyml, root);
}

/*
 * vnacal_save: create or overwrite a calibration file with new data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @pathname: calibration file name
 *
 *   The pathname and basename parameters work as in vnacal_load except
 *   that the $HOME/{pathname} directory is created if necessary.
 */
int vnacal_save(vnacal_t *vcp, const char *pathname)
{
    FILE *fp;
    yaml_document_t document;
    yaml_version_directive_t version = { 1, 1 };
    yaml_tag_directive_t tags[1];
    yaml_emitter_t emitter;
    bool delete_document = false;
    int t_root, t_properties, t_calibrations;

    if ((fp = fopen(pathname, "w")) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "fopen: %s: %s",
		vcp->vc_filename, strerror(errno));
	return -1;
    }
    free((void *)vcp->vc_filename);
    if ((vcp->vc_filename = strdup(pathname)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"strdup: %s", strerror(errno));
	goto error;
    }
    errno = 0;
    if (!yaml_document_initialize(&document, &version, &tags[0], &tags[0],
		0, 0)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_initialize: %s: %s",
		vcp->vc_filename, strerror(errno));
	goto error;
    }
    delete_document = true;
    t_root = yaml_document_add_mapping(&document, NULL,
	    YAML_BLOCK_MAPPING_STYLE);
    if (t_root == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_add_mapping: %s: %s",
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
    if ((t_calibrations = yaml_document_add_sequence(&document, NULL,
		    YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"yaml_document_add_sequence: %s: %s",
		vcp->vc_filename, strerror(errno));
	goto error;
    }
    if (add_mapping_entry(vcp, &document, t_root, "calibrations",
		t_calibrations) == -1) {
	goto error;
    }
    for (int ci = 0; ci < vcp->vc_calibration_allocation; ++ci) {
	vnacal_calibration_t *calp = vcp->vc_calibration_vector[ci];
	int t_calibration, t_value, t_matrix, t_data;
	vnacal_layout_t vl;
	const char *cp;

	if (calp == NULL) {
	    continue;
	}
	_vnacal_layout(&vl, calp->cal_type, calp->cal_rows, calp->cal_columns);
	t_calibration = yaml_document_add_mapping(&document, NULL,
		YAML_ANY_MAPPING_STYLE);
	if (t_calibration == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_mapping: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (yaml_document_append_sequence_item(&document, t_calibrations,
		    t_calibration) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_append_sequence_item: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if ((t_value = yaml_document_add_scalar(&document, NULL,
			(yaml_char_t *)calp->cal_name,
			strlen(calp->cal_name),
			YAML_ANY_SCALAR_STYLE)) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration,
		    "name", t_value) == -1) {
	    goto error;
	}
	cp = vnacal_type_to_name(calp->cal_type);
	if ((t_value = yaml_document_add_scalar(&document, NULL,
			(yaml_char_t *)cp, strlen(cp),
			YAML_ANY_SCALAR_STYLE)) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration,
		    "type", t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_integer(&document, calp->cal_rows)) == -1) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration,
		    "rows", t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_integer(&document, calp->cal_columns)) == -1) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration, "columns",
		    t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_integer(&document, calp->cal_frequencies)) == -1) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration, "frequencies",
		    t_value) == -1) {
	    goto error;
	}
	if ((t_value = add_complex(&document, calp->cal_z0,
			vcp->vc_dprecision)) == -1) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_scalar: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration,
		    "z0", t_value) == -1) {
	    goto error;
	}
	if ((t_matrix = yaml_document_add_sequence(&document, NULL,
			YAML_ANY_SEQUENCE_STYLE)) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_sequence: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	{
	    int t_calibration_properties;

	    if ((t_calibration_properties = add_properties(vcp, &document,
			    calp->cal_properties)) == -1) {
		goto error;
	    }
	    if (t_calibration_properties != 0) {
		if (add_mapping_entry(vcp, &document, t_calibration,
			    "properties", t_calibration_properties) == -1) {
		    goto error;
		}
	    }
	}
	if ((t_data = yaml_document_add_sequence(&document, NULL,
			YAML_BLOCK_SEQUENCE_STYLE)) == 0) {
	    if (errno == 0) {
		errno = EINVAL;
	    }
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "yaml_document_add_sequence: %s: %s",
		    vcp->vc_filename, strerror(errno));
	    goto error;
	}
	if (add_mapping_entry(vcp, &document, t_calibration, "data",
		    t_data) == -1) {
	    goto error;
	}
	for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
	    int t_per_frequency;
	    error_terms_common_arguments_t etca;

	    if ((t_per_frequency = yaml_document_add_mapping(&document, NULL,
			    YAML_ANY_MAPPING_STYLE)) == 0) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"yaml_document_add_sequence: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if (!yaml_document_append_sequence_item(&document, t_data,
			t_per_frequency)) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"yaml_document_append_sequence_item: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if ((t_value = add_double(&document,
			    calp->cal_frequency_vector[findex],
			    vcp->vc_fprecision)) == -1) {
		if (errno == 0) {
		    errno = EINVAL;
		}
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"yaml_document_add_scalar: %s: %s",
			vcp->vc_filename, strerror(errno));
		goto error;
	    }
	    if (add_mapping_entry(vcp, &document, t_per_frequency, "f",
			t_value) == -1) {
		goto error;
	    }
	    etca.etca_vcp = vcp;
	    etca.etca_document = &document;
	    etca.etca_mapping = t_per_frequency;
	    etca.etca_findex = findex;
	    if (add_error_parameters(etca, calp, &vl) == -1) {
		goto error;
	    }
	}
    }

    /*
     * Write the output file.
     */
    (void)fprintf(fp, "#VNACal 1.0\n");
    if (!yaml_emitter_initialize(&emitter)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM, "yaml_emitter_initialize: %s: %s",
		vcp->vc_filename, strerror(errno));
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
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM, "yaml_emitter_open: %s: %s",
		vcp->vc_filename, emitter.problem);
	goto error;
    }
    delete_document = false;	/* deleted by yaml_emitter_dump */
    if (!yaml_emitter_dump(&emitter, &document)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM, "yaml_emitter_dump: %s: %s",
		vcp->vc_filename, emitter.problem);
	goto error;
    }
    if (!yaml_emitter_close(&emitter)) {
	if (errno == 0) {
	    errno = EINVAL;
	}
	_vnacal_error(vcp, VNAERR_SYSTEM, "yaml_emitter_close: %s: %s",
		vcp->vc_filename, emitter.problem);
	goto error;
    }
    (void)yaml_emitter_delete(&emitter);
    if (fclose(fp) == -1) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "fclose: %s: %s",
		vcp->vc_filename, strerror(errno));
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
    free((void *)vcp->vc_filename);
    vcp->vc_filename = NULL;
    return -1;
}
