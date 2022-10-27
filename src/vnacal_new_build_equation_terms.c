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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_new_internal.h"

//#define DEBUG

#ifdef DEBUG
/*
 * print_term: print an expanded term of an equation
 *   @vntp:   structure representing an expanded term of an equation
 *   @vlp:    layout structure
 *   @sindex: system index
 *   @with_v: include V terms if true
 */
static void print_term(const vnacal_new_term_t *vntp,
	const vnacal_layout_t *vlp, int sindex, bool with_v)
{
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    const int v_columns = VL_V_COLUMNS(vlp);
    const int unity     = _vl_unity_offset(vlp, sindex);
    int term = vntp->vnt_xindex;
    bool negative = vntp->vnt_negative;
    bool first_factor = true;

    /*
     * Move the unity term back to the left-hand side.
     */
    term -= sindex * (vlp->vl_t_terms - 1);
    if (term == -1) {
	term = unity;
	negative = !negative;
    } else if (term >= unity) {
	++term;
    }

    /*
     * Print
     */
    if (negative) {
	(void)printf(" -");
    } else {
	(void)printf(" +");
    }
    if (vntp->vnt_m_cell >= 0) {
	const int m_cell   = vntp->vnt_m_cell;
	const int m_row    = m_cell / m_columns;
	const int m_column = m_cell % m_columns;

	(void)printf("m%d%s%d",
		m_row + 1, m_row >= 9 ? "," : "", m_column + 1);
	first_factor = false;
    }
    if (vntp->vnt_s_cell >= 0) {
	const int s_cell  = vntp->vnt_s_cell;
	const int s_row    = s_cell / s_columns;
	const int s_column = s_cell % s_columns;

	if (!first_factor) {
	    (void)printf("*");
	}
	(void)printf("s%d%s%d",
		s_row + 1, s_row >= 9 ? "," : "", s_column + 1);
	first_factor = false;
    }
    if (with_v && vntp->vnt_v_cell >= 0) {
	const int v_cell   = vntp->vnt_v_cell;
	const int v_row    = v_cell / v_columns;
	const int v_column = v_cell % v_columns;

	if (!first_factor) {
	    (void)printf("*");
	}
	(void)printf("v%d%s%d",
		v_row + 1, v_row >= 9 ? "," : "", v_column + 1);
	first_factor = false;
    }
    if (!first_factor) {
	(void)printf("*");
    }
    if (VL_IS_T(vlp)) {	/* T terms */
	if (term < VL_TI_OFFSET(vlp)) {		/* ts */
	    const int ts_cell    = term;
	    const int ts_row     = ts_cell / VL_TS_COLUMNS(vlp);
	    const int ts_column  = ts_cell % VL_TS_COLUMNS(vlp);

	    (void)printf("ts%d%s%d",
		ts_row + 1, ts_row >= 9 ? "," : "", ts_column + 1);

	} else if (term < VL_TX_OFFSET(vlp)) {	/* ti */
	    const int ti_cell    = term    - VL_TI_OFFSET(vlp);
	    const int ti_row     = ti_cell / VL_TI_COLUMNS(vlp);
	    const int ti_column  = ti_cell % VL_TI_COLUMNS(vlp);

	    (void)printf("ti%d%s%d",
		ti_row + 1, ti_row >= 9 ? "," : "", ti_column + 1);

	} else if (term < VL_TM_OFFSET(vlp)) {	/* tx */
	    const int tx_cell    = term    - VL_TX_OFFSET(vlp);
	    const int tx_row     = tx_cell / VL_TX_COLUMNS(vlp);
	    const int tx_column  = tx_cell % VL_TX_COLUMNS(vlp);

	    (void)printf("tx%d%s%d",
		tx_row + 1, tx_row >= 9 ? "," : "", tx_column + 1);

	} else {				/* tm */
	    const int tm_cell    = term    - VL_TM_OFFSET(vlp);
	    const int tm_row     = tm_cell / VL_TM_COLUMNS(vlp);
	    const int tm_column  = tm_cell % VL_TM_COLUMNS(vlp);

	    (void)printf("tm%d%s%d",
		tm_row + 1, tm_row >= 9 ? "," : "", tm_column + 1);
	}

    } else {		/* U terms */
	if (term < VL_UI_OFFSET(vlp)) {		/* um */
	    const int um_cell    = term;
	    const int um_row     = um_cell / VL_UM_COLUMNS(vlp);
	    const int um_column  = um_cell % VL_UM_COLUMNS(vlp);

	    (void)printf("um%d%s%d",
		um_row + 1, um_row >= 9 ? "," : "", um_column + 1);

	} else if (term < VL_UX_OFFSET(vlp)) {	/* ui */
	    const int ui_cell    = term    - VL_UI_OFFSET(vlp);
	    const int ui_row     = ui_cell / VL_UI_COLUMNS(vlp);
	    const int ui_column  = ui_cell % VL_UI_COLUMNS(vlp);

	    (void)printf("ui%d%s%d",
		ui_row + 1, ui_row >= 9 ? "," : "", ui_column + 1);

	} else if (term < VL_US_OFFSET(vlp)) {	/* ux */
	    const int ux_cell    = term    - VL_UX_OFFSET(vlp);
	    const int ux_row     = ux_cell / VL_UX_COLUMNS(vlp);
	    const int ux_column  = ux_cell % VL_UX_COLUMNS(vlp);

	    (void)printf("ux%d%s%d",
		ux_row + 1, ux_row >= 9 ? "," : "", ux_column + 1);

	} else {				/* us */
	    const int us_cell    = term    - VL_US_OFFSET(vlp);
	    const int us_row     = us_cell / VL_US_COLUMNS(vlp);
	    const int us_column  = us_cell % VL_US_COLUMNS(vlp);

	    (void)printf("us%d%s%d",
		us_row + 1, us_row >= 9 ? "," : "", us_column + 1);
	}
    }
}

/*
 * print_equation: print the generated no-V and V equations
 *   @vnep: equation structure
 */
static void print_equation(const vnacal_new_equation_t *vnep)
{
    const vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    const vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const vnacal_new_term_t *vntp;
    int sindex = 0;

    (void)printf("eq%d%s%d\n",
	vnep->vne_row + 1, vnep->vne_row >= 9 ? "," : "",
	vnep->vne_column + 1);
    if (VL_HAS_COLUMN_SYSTEMS(vlp)) {
	sindex = vnep->vne_column;
	(void)printf("sindex %2d\n", sindex + 1);
    }
    (void)printf("no-v:");
    for (vntp = vnep->vne_term_list_no_v; vntp != NULL;
	    vntp = vntp->vnt_next_no_v) {
	print_term(vntp, vlp, sindex, false);
    }
    (void)printf(" == 0\n");
    (void)printf("v:   ");
    for (vntp = vnep->vne_term_list; vntp != NULL;
	    vntp = vntp->vnt_next) {
	print_term(vntp, vlp, sindex, true);
    }
    (void)printf(" == 0\n");
    (void)printf("\n");
}
#endif /* DEBUG */

/*
 * add_term: add a term to the current equation
 *   @vnmp: represents a measured calibration standard
 *   @vncppp_anchor: address of pointer where next term should be linked
 *   @xindex: column in expanded coefficient matrix or -1 for right hand side
 *   @v_columns: number of columns in v_matrix
 *   @negative: true if term has a minus sign
 *   @m: index of measurement in vnsm_m_matrix or -1
 *   @s: index of parameter   in vnsm_s_matrix or -1
 *   @v: index of value       in vnsp_v_matrix or -1
 */
static int add_term(vnacal_new_measurement_t *vnmp,
	vnacal_new_term_t ***vncppp_anchor, int xindex,
	int v_columns, bool negative, int m, int s, int v)
{
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    vnacal_new_term_t *vntp;

    if ((vntp = malloc(sizeof(vnacal_new_term_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return -1;
    }
    (void)memset((void *)vntp, 0, sizeof(vnacal_new_term_t));
    vntp->vnt_xindex = xindex;
    vntp->vnt_negative = negative;
    vntp->vnt_m_cell = m;
    vntp->vnt_s_cell = s;
    vntp->vnt_v_cell = v;
    if (v % (v_columns + 1) == 0) {
	*vncppp_anchor[0] = vntp;
	vncppp_anchor[0] = &vntp->vnt_next_no_v;
    }
    *vncppp_anchor[1] = vntp;
    vncppp_anchor[1] = &vntp->vnt_next;
    return 0;
}

/*
 * build_terms_t8: build coefficients for T8/TE10 error terms
 *   @vnep: vnacal_new_equation_t structure
 *
 * Build the coefficients of: -Ts S V - Ti V + M Tx S V + M Tm V == 0
 * for a single equation.
 *
 * Dimensions (m_rows <= m_columns)
 *   ts: m_rows    x m_columns (diagonal)
 *   ti: m_rows    x m_columns (diagonal)
 *   tx: m_columns x m_columns (diagonal)
 *   tm: m_columns x m_columns (diagonal)
 *   eq: m_rows    x m_columns
 *   m:  m_rows    x m_columns
 *   s:  m_columns x m_columns
 *   v:  m_columns x m_columns
 */
static int build_terms_t8(vnacal_new_equation_t *vnep)
{
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int eq_row    = vnep->vne_row;
    const int eq_column = vnep->vne_column;
    vnacal_new_term_t **anchors[2] = {
	&vnep->vne_term_list_no_v,
	&vnep->vne_term_list
    };
    int base_coefficient = 0;

    /*
     * Add the non-zero Ts terms.
     */
    for (int v_row = 0; v_row < m_columns; ++v_row) {
	const int s_cell = eq_row * m_columns + v_row;
	const int v_cell = v_row  * m_columns + eq_column;
	vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	if (!vnmp->vnm_connectivity_matrix[v_cell]) {
	    continue;
	}
	assert(vnprp != NULL);
	if (vnprp != vnp->vn_zero) {
	    if (add_term(vnmp, anchors, base_coefficient + eq_row,
			/*v_columns*/m_columns, /*negative*/true,
			/*m*/-1, s_cell, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_rows;

    /*
     * Add the Ti term.
     */
    {
	const int v_cell = eq_row * m_columns + eq_column;

	if (vnmp->vnm_connectivity_matrix[v_cell]) {
	    if (add_term(vnmp, anchors, base_coefficient + eq_row,
			/*v_columns*/m_columns, /*negative*/true,
			/*m*/-1, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_rows;

    /*
     * Add the non-zero Tx terms.
     */
    for (int tx_d = 0; tx_d < m_columns; ++tx_d) {
	const int m_cell = eq_row * m_columns + tx_d;

	for (int v_row = 0; v_row < m_columns; ++v_row) {
	    const int s_cell = tx_d  * m_columns + v_row;
	    const int v_cell = v_row * m_columns + eq_column;
	    vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	    if (!vnmp->vnm_connectivity_matrix[v_cell]) {
		continue;
	    }
	    assert(vnprp != NULL);
	    if (vnprp != vnp->vn_zero) {
		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (add_term(vnmp, anchors, base_coefficient + tx_d,
			    /*v_columns*/m_columns, /*negative*/false,
			    m_cell, s_cell, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_columns;

    /*
     * Add the non-zero Tm terms.
     */
    for (int tm_d = 0; tm_d < m_columns; ++tm_d) {
	const int m_cell = eq_row * m_columns + tm_d;
	const int v_cell = tm_d * m_columns + eq_column;

	if (!vnmp->vnm_connectivity_matrix[v_cell]) {
	    continue;
	}
	assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	if (tm_d == 0) {	/* tm11 == 1.0 */
	    if (add_term(vnmp, anchors, -1,
			/*v_columns*/m_columns, /*negative*/true,
			m_cell, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	} else {
	    if (add_term(vnmp, anchors, base_coefficient + tm_d - 1,
			/*v_columns*/m_columns, /*negative*/false,
			m_cell, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_columns - 1;

    return 0;
}

/*
 * build_terms_u8: build coefficients for U8/UE10 error terms
 *   @vnep: vnacal_new_equation_t structure
 *
 * Build the coefficients of: V Um M + V Ui - V S Ux M - V S Us == 0
 * for a single equation.
 *
 * Dimensions (m_rows >= m_columns)
 *   um: m_rows    x m_rows    (diagonal)
 *   ui: m_rows    x m_columns (diagonal)
 *   ux: m_rows    x m_rows    (diagonal)
 *   us: m_rows    x m_columns (diagonal)
 *   eq: m_rows    x m_columns
 *   m:  m_rows    x m_columns
 *   s:  m_rows    x m_rows
 *   v:  m_rows    x m_rows
 */
static int build_terms_u8(vnacal_new_equation_t *vnep)
{
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int eq_row    = vnep->vne_row;
    const int eq_column = vnep->vne_column;
    vnacal_new_term_t **anchors[2] = {
	&vnep->vne_term_list_no_v,
	&vnep->vne_term_list
    };
    int base_coefficient = 0;

    /*
     * Add the Um terms.
     */
    for (int um_d = 0; um_d < m_rows; ++um_d) {
	const int v_cell = eq_row * m_rows + um_d;
	const int m_cell = um_d * m_columns + eq_column;

	if (!vnmp->vnm_connectivity_matrix[v_cell]) {
	    continue;
	}
	assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	if (um_d == 0) { /* um11 == 1.0 */
	    if (add_term(vnmp, anchors, -1,
			/*v_columns*/m_rows, /*negative*/true,
			m_cell, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	} else {
	    if (add_term(vnmp, anchors, base_coefficient + um_d - 1,
			/*v_columns*/m_rows, /*negative*/false,
			m_cell, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_rows - 1;

    /*
     * Add the Ui term.
     */
    {
	const int v_cell = eq_row * m_rows + eq_column;

	if (vnmp->vnm_connectivity_matrix[v_cell]) {
	    if (add_term(vnmp, anchors, base_coefficient + eq_column,
			/*v_columns*/m_rows, /*negative*/false,
			/*m*/-1, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_columns;

    /*
     * Add the non-zero Ux terms.
     */
    for (int ux_d = 0; ux_d < m_rows; ++ux_d) {
	const int m_cell = ux_d * m_columns + eq_column;

	for (int v_column = 0; v_column < m_rows; ++v_column) {
	    const int v_cell = eq_row   * m_rows + v_column;
	    const int s_cell = v_column * m_rows + ux_d;
	    vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	    if (!vnmp->vnm_connectivity_matrix[v_cell]) {
		continue;
	    }
	    assert(vnprp != NULL);
	    if (vnprp != vnp->vn_zero) {
		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (add_term(vnmp, anchors, base_coefficient + ux_d,
			    /*v_columns*/m_rows, /*negative*/true,
			    m_cell, s_cell, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_rows;

    /*
     * Add the non-zero Us terms.
     */
    for (int v_column = 0; v_column < m_rows; ++v_column) {
	const int v_cell = eq_row * m_rows + v_column;
	const int s_cell = v_column * m_rows + eq_column;
	vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	if (!vnmp->vnm_connectivity_matrix[v_cell]) {
	    continue;
	}
	assert(vnprp != NULL);
	if (vnprp != vnp->vn_zero) {
	    if (add_term(vnmp, anchors, base_coefficient + eq_column,
			/*v_columns*/m_rows, /*negative*/true,
			/*m*/-1, s_cell, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_columns;

    return 0;
}

/*
 * build_terms_t16: build coefficients for T16 error terms
 *   @vnep: vnacal_new_equation_t structure
 *
 * Build the coefficients of: -Ts S V - Ti V + M Tx S V + M Tm V == 0
 * for a single equation.
 *
 * Dimensions (m_rows <= m_columns)
 *   ts: m_rows    x m_columns
 *   ti: m_rows    x m_columns
 *   tx: m_columns x m_columns
 *   tm: m_columns x m_columns
 *   eq: m_rows    x m_columns
 *   m:  m_rows    x m_columns
 *   s:  m_columns x m_columns
 *   v:  m_columns x m_columns
 */
static int build_terms_t16(vnacal_new_equation_t *vnep)
{
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int eq_row    = vnep->vne_row;
    const int eq_column = vnep->vne_column;
    vnacal_new_term_t **anchors[2] = {
	&vnep->vne_term_list_no_v,
	&vnep->vne_term_list
    };
    int base_coefficient = 0;

    /*
     * Add the non-zero Ts terms.
     */
    for (int ts_column = 0; ts_column < m_columns; ++ts_column) {
	const int ts_cell = eq_row * m_columns + ts_column;

	for (int v_row = 0; v_row < m_columns; ++v_row) {
	    const int s_cell = ts_column * m_columns + v_row;
	    const int v_cell = v_row * m_columns + eq_column;
	    vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	    if (vnprp != NULL && vnprp == vnp->vn_zero) {
		continue;
	    }
	    if (add_term(vnmp, anchors, base_coefficient + ts_cell,
			/*v_columns*/m_columns, /*negative*/true,
			/*m*/-1, s_cell, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_rows * m_columns;

    /*
     * Add the Ti terms.
     */
    for (int ti_column = 0; ti_column < m_columns; ++ti_column) {
	const int ti_cell = eq_row * m_columns + ti_column;
	const int v_cell = ti_column * m_columns + eq_column;

	if (add_term(vnmp, anchors, base_coefficient + ti_cell,
		    /*v_columns*/m_columns, /*negative*/true,
		    /*m*/-1, /*s*/-1, v_cell) == -1) {
	    return -1;
	}
    }
    base_coefficient += m_rows * m_columns;

    /*
     * Add the non-zero Tx terms.
     */
    for (int tx_row = 0; tx_row < m_columns; ++tx_row) {
	for (int tx_column = 0; tx_column < m_columns; ++tx_column) {
	    const int tx_cell = tx_row * m_columns + tx_column;
	    const int m_cell = eq_row * m_columns + tx_row;

	    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	    for (int v_row = 0; v_row < m_columns; ++v_row) {
		const int v_cell = v_row * m_columns + eq_column;
		const int s_cell = tx_column * m_columns + v_row;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		if (vnprp != NULL && vnprp == vnp->vn_zero) {
		    continue;
		}
		if (add_term(vnmp, anchors, base_coefficient + tx_cell,
			    /*v_columns*/m_columns, /*negative*/false,
			    m_cell, s_cell, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_columns * m_columns;

    /*
     * Add the Tm terms.
     */
    for (int tm_row = 0; tm_row < m_columns; ++tm_row) {
	for (int tm_column = 0; tm_column < m_columns; ++tm_column) {
	    const int tm_cell = tm_row * m_columns + tm_column;
	    const int m_cell = eq_row * m_columns + tm_row;
	    const int v_cell = tm_column * m_columns + eq_column;

	    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	    if (tm_cell == 0) {	/* tm11 == 1.0 */
		if (add_term(vnmp, anchors, -1,
			    /*v_columns*/m_columns, /*negative*/true,
			    m_cell, /*s*/-1, v_cell) == -1) {
		    return -1;
		}
	    } else {
		if (add_term(vnmp, anchors, base_coefficient + tm_cell - 1,
			    /*v_columns*/m_columns, /*negative*/false,
			    m_cell, /*s*/-1, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_columns * m_columns - 1;

    return 0;
}

/*
 * build_terms_u16: build coefficients for U16 error terms
 *   @vnep: vnacal_new_equation_t structure
 *
 * Build the coefficients of: V Um M + V Ui - V S Ux M - V S Us == 0
 * for a single equation.
 *
 * Dimensions (m_rows >= m_columns)
 *   um: m_rows    x m_rows
 *   ui: m_rows    x m_columns
 *   ux: m_rows    x m_rows
 *   us: m_rows    x m_columns
 *   eq: m_rows    x m_columns
 *   m:  m_rows    x m_columns
 *   s:  m_rows    x m_rows
 *   v:  m_rows    x m_rows
 */
static int build_terms_u16(vnacal_new_equation_t *vnep)
{
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int eq_row    = vnep->vne_row;
    const int eq_column = vnep->vne_column;
    vnacal_new_term_t **anchors[2] = {
	&vnep->vne_term_list_no_v,
	&vnep->vne_term_list
    };
    int base_coefficient = 0;

    /*
     * Add the Um terms.
     */
    for (int um_row = 0; um_row < m_rows; ++um_row) {
	for (int um_column = 0; um_column < m_rows; ++um_column) {
	    const int um_cell = um_row * m_rows + um_column;
	    const int v_cell = eq_row * m_rows + um_row;
	    const int m_cell = um_column * m_columns + eq_column;

	    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	    if (um_cell == 0) {	/* um11 == 1.0 */
		if (add_term(vnmp, anchors, -1,
			    /*v_columns*/m_rows, /*negative*/true,
			    m_cell, /*s*/-1, v_cell) == -1) {
		    return -1;
		}
	    } else {
		if (add_term(vnmp, anchors, base_coefficient + um_cell - 1,
			    /*v_columns*/m_rows, /*negative*/false,
			    m_cell, /*s*/-1, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_rows * m_rows - 1;

    /*
     * Add the Ui terms.
     */
    for (int ui_row = 0; ui_row < m_rows; ++ui_row) {
	const int ui_cell = ui_row * m_columns + eq_column;
	const int v_cell  = eq_row * m_rows + ui_row;

	if (add_term(vnmp, anchors, base_coefficient + ui_cell,
		    /*v_columns*/m_rows, /*negative*/false,
		    /*m*/-1, /*s*/-1, v_cell) == -1) {
	    return -1;
	}
    }
    base_coefficient += m_rows * m_columns;

    /*
     * Add the non-zero Ux terms.
     */
    for (int ux_row = 0; ux_row < m_rows; ++ux_row) {
	for (int ux_column = 0; ux_column < m_rows; ++ux_column) {
	    const int ux_cell = ux_row * m_rows + ux_column;
	    const int m_cell = ux_column * m_columns + eq_column;

	    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	    for (int v_column = 0; v_column < m_rows; ++v_column) {
		const int v_cell = eq_row * m_rows + v_column;
		const int s_cell = v_column * m_rows + ux_row;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		if (vnprp != NULL && vnprp == vnp->vn_zero) {
		    continue;
		}
		if (add_term(vnmp, anchors, base_coefficient + ux_cell,
			    m_rows, /*negative*/true,
			    m_cell, s_cell, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_rows * m_rows;

    /*
     * Add the non-zero Us terms.
     */
    for (int us_row = 0; us_row < m_rows; ++us_row) {
	const int us_cell = us_row * m_columns + eq_column;

	for (int v_column = 0; v_column < m_rows; ++v_column) {
	    const int v_cell = eq_row * m_rows + v_column;
	    const int s_cell = v_column * m_rows + us_row;
	    vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	    if (vnprp != NULL && vnprp == vnp->vn_zero) {
		continue;
	    }
	    if (add_term(vnmp, anchors, base_coefficient + us_cell,
			/*v_columns*/m_rows, /*negative*/true,
			/*m*/-1, s_cell, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_rows * m_columns;

    return 0;
}

/*
 * build_terms_ue14: build coefficients for UE14 error terms
 *   @vnep: vnacal_new_equation_t structure
 *
 * Build the coefficients of: V Um M + V Ui - V S Ux M - V S Us == 0
 * for a single equation.
 *
 * Dimensions (m_rows >= m_columns)
 *   um: m_rows    x m_rows    (diagonal)
 *   ui: m_rows    x 1         (diagonal)
 *   ux: m_rows    x m_rows    (diagonal)
 *   us: m_rows    x 1         (diagonal)
 *   eq: m_rows    x m_columns (each column belongs to an independent system)
 *   m:  m_rows    x m_columns (each column belongs to an independent system)
 *   s:  m_rows    x m_rows
 *   v:  m_rows    x m_rows
 */
static int build_terms_ue14(vnacal_new_equation_t *vnep)
{
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int eq_row    = vnep->vne_row;
    const int eq_column = vnep->vne_column;
    vnacal_new_term_t **anchors[2] = {
	&vnep->vne_term_list_no_v,
	&vnep->vne_term_list
    };
    int base_coefficient = 0;

    /*
     * Add the Um terms.
     */
    for (int um_d = 0; um_d < m_rows; ++um_d) {
	const int m_cell = um_d * m_columns + eq_column;
	const int v_cell = eq_row * m_rows + um_d;

	if (!vnmp->vnm_connectivity_matrix[v_cell]) {
	    continue;
	}
	assert(vnmp->vnm_m_matrix[m_cell] != NULL);
	if (um_d == eq_column) {
	    if (add_term(vnmp, anchors, -1,
			/*v_columns*/m_rows, /*negative*/true,
			m_cell, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	    --base_coefficient;
	} else {
	    if (add_term(vnmp, anchors, base_coefficient + um_d,
			/*v_columns*/m_rows, /*negative*/false,
			m_cell, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += m_rows;

    /*
     * Add the Ui term.
     */
    {
	const int v_cell = eq_row * m_rows + eq_column;

	if (vnmp->vnm_connectivity_matrix[v_cell]) {
	    if (add_term(vnmp, anchors, base_coefficient,
			/*v_columns*/m_rows, /*negative*/false,
			/*m*/-1, /*s*/-1, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += 1;

    /*
     * Add the non-zero Ux terms.
     */
    for (int ux_d = 0; ux_d < m_rows; ++ux_d) {
	const int m_cell = ux_d * m_columns + eq_column;

	for (int v_column = 0; v_column < m_rows; ++v_column) {
	    const int v_cell = eq_row * m_rows + v_column;
	    const int s_cell = v_column * m_rows + ux_d;
	    vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

	    if (!vnmp->vnm_connectivity_matrix[v_cell]) {
		continue;
	    }
	    assert(vnprp != NULL);
	    if (vnprp != vnp->vn_zero) {
		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (add_term(vnmp, anchors, base_coefficient + ux_d,
			    /*v_columns*/m_rows, /*negative*/true,
			    m_cell, s_cell, v_cell) == -1) {
		    return -1;
		}
	    }
	}
    }
    base_coefficient += m_rows;

    /*
     * Add the non-zero Us terms.
     */
    for (int v_column = 0; v_column < m_rows; ++v_column) {
	const int v_cell = eq_row * m_rows + v_column;
	const int s_cell = v_column * m_rows + eq_column;
	vnacal_new_parameter_t *vnprp =
	    vnmp->vnm_s_matrix[s_cell];

	if (!vnmp->vnm_connectivity_matrix[v_cell]) {
	    continue;
	}
	assert(vnprp != NULL);
	if (vnprp != vnp->vn_zero) {
	    if (add_term(vnmp, anchors, base_coefficient,
			/*v_columns*/m_rows, /*negative*/true,
			/*m*/-1, s_cell, v_cell) == -1) {
		return -1;
	    }
	}
    }
    base_coefficient += 1;

    return 0;
}

/*
 * _vnacal_new_build_equation_terms: build vnacal_new_term_t lists
 *   @vnep: equation generated from a measurement of a standard
 *
 * Given a partially filled vnacal_new_equation_t structure, build lists
 * of expanded algebraic terms making the equation.  For example with T
 * error terms, we have the following matrix equation representing the
 * equations for a measured standard:
 *
 *     -Ts S V - Ti V + M Tx S V + M Tm V == 0
 *
 * Where the Ts, Ti, Tx, and Tm matrices are the t11, t12, t21, and
 * t22, error term elements we need to solve for, respectively, S is
 * the s-parameters matrix of the calibration standard, M is the matrix
 * of measurements as seen by the VNA, and V is a weighting matrix that
 * transforms the residuals of the equations to errors in the M matrix.
 *
 * For example take the case of 2x2 T8 parameters.  If we multiply
 * out the matrix equation above, we get four equations:
 *
 *     M11 equation:
 *      -ts11 s11 v11 - ts11 s12 v21 - ti11 v11 +
 *          m11 tx11 s11 v11 + m11 tx11 s12 v21 +
 *          m12 tx22 s21 v11 + m12 tx22 s22 v21 +
 *          m11 tm11 v11 + m12 tm22 v21        == 0
 *
 *     M12 equation:
 *      -s11 ts11 v12 - s12 ts11 v22 - ti11 v12 +
 *          m11 tx11 s11 v12 + m11 tx11 s12 v22 +
 *          m12 tx22 s21 v12 + m12 tx22 s22 v22 +
 *          m11 tm11 v12 + m12 tm22 v22        == 0
 *
 *     M21 equation:
 *      -s21 ts22 v11 - s22 ts22 v21 - ti22 v21 +
 *          m21 tx11 s11 v11 + m21 tx11 s12 v21 +
 *          m22 tx22 s21 v11 + m22 tx22 s22 v21 +
 *          m21 tm11 v11 + m22 tm22 v21        == 0
 *
 *     M22 equation:
 *      -s21 ts22 v12 - s22 ts22 v22 - ti22 v22 +
 *          m21 tx11 s11 v12 + m21 tx11 s12 v22 +
 *          m22 tx22 s21 v12 + m22 tx22 s22 v22 +
 *          m21 tm11 v12 + m22 tm22 v22        == 0
 *
 * Which of these we generate depends on vne_row and vne_column in the
 * vnacal_new_equation_t structure.  In addition to these coefficients,
 * we also build a thread through the list elements representing the
 * subset along the major diagonal of the V matrix.  This is useful when
 * we're not using the V matrices and assume V is the identity matrix.
 *
 * For T error terms, we set tm11=1 and move the associated terms of
 * the equation to the right hand side.  For U erorr terms, we similarly
 * set um=11 and move associated terms to the right.
 */
int _vnacal_new_build_equation_terms(vnacal_new_equation_t *vnep)
{
    vnacal_new_measurement_t *vnmp = vnep->vne_vnmp;
    vnacal_new_t *vnp = vnmp->vnm_vnp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;

    /*
     * Add terms based on error term type.
     */
    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	if (build_terms_t8(vnep) == -1) {
	    return -1;
	}
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
	if (build_terms_u8(vnep) == -1) {
	    return -1;
	}
	break;

    case VNACAL_T16:
	if (build_terms_t16(vnep) == -1) {
	    return -1;
	}
	break;

    case VNACAL_U16:
	if (build_terms_u16(vnep) == -1) {
	    return -1;
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	if (build_terms_ue14(vnep) == -1) {
	    return -1;
	}
	break;

    case VNACAL_E12:
    default:
	abort();
    }
#ifdef DEBUG
    print_equation(vnep);
#endif /* DEBUG */
    return 0;
}
