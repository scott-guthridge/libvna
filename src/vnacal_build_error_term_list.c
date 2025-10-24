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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <yaml.h>
#include "vnacal_internal.h"


/*
 * make_error_term_matrix: create a new error term matrix and append to list
 *   @vcp: vnacal structure
 *   @anchor: address of address telling where to insert the next element
 *   @type: error term matrix type
 *   @name: name of error term matrix (must remain valid)
 *   @matrix: error term matrix/vector of vectors (one entry per findex)
 *   @rows: rows in matrix
 *   @columns: columns in matrix
 */
static int make_error_term_matrix(vnacal_calibration_t *calp,
	vnacal_error_term_matrix_t ***anchor,
	vnacal_error_term_matrix_type_t type, const char *name,
	double complex **matrix, int rows, int columns)
{
    vnacal_t *vcp = calp->cal_vcp;
    const int cells = rows * columns;
    vnacal_error_term_matrix_t *vetmp = NULL;

    assert(type != VETM_VECTOR || rows == 1);
    if ((vetmp = malloc(sizeof(vnacal_error_term_matrix_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	goto error;
    }
    (void)memset((void *)vetmp, 0, sizeof(*vetmp));
    vetmp->vetm_calp = calp;
    vetmp->vetm_type = type;
    vetmp->vetm_name = name;
    if ((vetmp->vetm_matrix = calloc(cells,
		    sizeof(double complex *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto error;
    }
    vetmp->vetm_rows = rows;
    vetmp->vetm_columns = columns;
    vetmp->vetm_next = NULL;
    (void)memcpy((void *)vetmp->vetm_matrix, (void *)matrix,
	    cells * sizeof(double complex *));
    **anchor = vetmp;
    *anchor = &vetmp->vetm_next;
    return 0;

error:
    _vnacal_free_error_term_matrices(&vetmp);
    return -1;
}

/*
 * _vnacal_build_error_term_list: make list of ET matrices for load/save
 *   @calp: calibration structure
 *   @head: address to hold head of generated list
 */
int _vnacal_build_error_term_list(vnacal_calibration_t *calp,
		const vnacal_layout_t *vlp, vnacal_error_term_matrix_t **head)
{
    double complex **e = calp->cal_error_term_vector;
    vnacal_error_term_matrix_t **anchor = head;
    int rc = -1;

    *head = NULL;
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
	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "ts",
			ts, 1, ts_terms) == -1) {
		goto out;
	    }

	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "ti",
			ti, 1, ti_terms) == -1) {
		goto out;
	    }

	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "tx",
			tx, 1, tx_terms) == -1) {
		goto out;
	    }

	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "tm",
			tm, 1, tm_terms) == -1) {
		goto out;
	    }
	}
	if (calp->cal_type == VNACAL_TE10) {
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex **el = &e[el_offset];

	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX_ND, "el",
			el, el_rows, el_columns) == -1) {
		goto out;
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
	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "um",
			um, 1, um_terms) == -1) {
		goto out;
	    }

	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "ui",
			ui, 1, ui_terms) == -1) {
		goto out;
	    }

	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "ux",
			ux, 1, ux_terms) == -1) {
		goto out;
	    }

	    if (make_error_term_matrix(calp, &anchor, VETM_VECTOR, "us",
			us, 1, us_terms) == -1) {
		goto out;
	    }
	}
	if (calp->cal_type == VNACAL_UE10) {
	    const int el_rows    = VL_EL_ROWS(vlp);
	    const int el_columns = VL_EL_COLUMNS(vlp);
	    const int el_offset  = VL_EL_OFFSET(vlp);
	    double complex **el = &e[el_offset];

	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX_ND, "el",
			el, el_rows, el_columns) == -1) {
		goto out;
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
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "ts",
			ts, ts_rows, ts_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "ti",
			ti, ti_rows, ti_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "tx",
			tx, tx_rows, tx_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "tm",
			tm, tm_rows, tm_columns) == -1) {
		goto out;
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
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "um",
			um, um_rows, um_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "ui",
			ui, ui_rows, ui_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "ux",
			ux, ux_rows, ux_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "us",
			us, us_rows, us_columns) == -1) {
		goto out;
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
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "um",
			&packed_um[0][0], um_terms, m_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "ui",
			&packed_ui[0][0], ui_terms, m_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "ux",
			&packed_ux[0][0], ux_terms, m_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "us",
			&packed_us[0][0], us_terms, m_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX_ND, "el",
			el, el_rows, el_columns) == -1) {
		goto out;
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
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "el",
			&packed_el[0][0], el_terms, m_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "er",
			&packed_er[0][0], er_terms, m_columns) == -1) {
		goto out;
	    }
	    if (make_error_term_matrix(calp, &anchor, VETM_MATRIX, "em",
			&packed_em[0][0], em_terms, m_columns) == -1) {
		goto out;
	    }
	}
	return 0;

    default:
	abort();
    }
    rc = 0;

out:
    if (rc == -1) {
	_vnacal_free_error_term_matrices(head);
    }
    return rc;
}

/*
 * _vnacal_free_error_term_matrices: free the memory for an error term matrix
 *   @vetmp: error term matrix structure
 */
void _vnacal_free_error_term_matrices(vnacal_error_term_matrix_t **vetmpp)
{
    vnacal_error_term_matrix_t *vetmp;

    while ((vetmp = *vetmpp) != NULL) {
	*vetmpp = vetmp->vetm_next;
	free((void *)vetmp->vetm_matrix);
	free((void *)vetmp);
    }
}
