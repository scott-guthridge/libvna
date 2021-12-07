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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/*
 * int_cmp: compare integers for qsort
 *   @vp1: first element
 *   @vp2: second element
 */
static int int_cmp(const void *vp1, const void *vp2)
{
    return *(const int *)vp1 - *(const int *)vp2;
}

/*
 * add_term: add a term to the current equation
 *   @function: name of user-called function
 *   @vnmp: represents a measured calibration standard
 *   @vncppp_anchor: address of pointer where next term should be linked
 *   @coefficient: column of this term in the coefficient matrix
 *   @negative: true if coefficient has a minus sign
 *   @m: index of measurement in vnm_m_matrix, or -1
 *   @s: index of parameter in vnm_s_matrix, or -1
 */
static int add_term(const char *function, vnacal_new_measurement_t *vnmp,
	vnacal_new_coefficient_t ***vncppp_anchor, int coefficient,
	bool negative, int m, int s)
{
    vnacal_new_t *vnp = vnmp->vnm_ncp;
    vnacal_t *vcp = vnp->vn_vcp;
    vnacal_new_coefficient_t *vncp;

    if ((vncp = malloc(sizeof(vnacal_new_coefficient_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return -1;
    }
    (void)memset((void *)vncp, 0, sizeof(vnacal_new_coefficient_t));
    vncp->vnc_coefficient = coefficient;
    vncp->vnc_negative = negative;
    vncp->vnc_m_cell = m;
    vncp->vnc_s_cell = s;
    **vncppp_anchor = vncp;
    *vncppp_anchor = &vncp->vnc_next;

    return 0;
}

/*
 * add_equation: add an equation to an vnacal_new_measurement_t structure
 *   @function: name of user-called function
 *   @vnmp: represents a measured calibration standard
 *   @vneppp_anchor: address of pointer where next equation should be linked
 *   @eq_row: row in matrix of equations
 *   @eq_column: column in matrix of equations
 */
static int add_equation(const char *function, vnacal_new_measurement_t *vnmp,
	vnacal_new_equation_t ***vneppp_anchor, int eq_row, int eq_column)
{
    vnacal_new_t *vnp = vnmp->vnm_ncp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &vnp->vn_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    vnacal_new_equation_t *vnep = NULL;
    vnacal_new_coefficient_t **vncpp_anchor = NULL;
    int base_coefficient = 0;

    /*
     * Construct the equation structure and link it onto the
     * vnacal_new_measurement_t structure.
     */
    vnep = malloc(sizeof(vnacal_new_equation_t));
    if (vnep == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return -1;
    }
    (void)memset((void *)vnep, 0, sizeof(vnacal_new_equation_t));
    vnep->vne_vnmp = vnmp;
    vnep->vne_row = eq_row;
    vnep->vne_column = eq_column;
    vncpp_anchor = &vnep->vne_coefficient_list;
    **vneppp_anchor = vnep;
    *vneppp_anchor = &vnep->vne_next;

    /*
     * Add terms based on error term type.
     */
    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	{
	    const int ts_diagonals = MIN(m_rows, s_rows);
	    const int ti_diagonals = MIN(m_rows, s_columns);
	    const int tx_diagonals = MIN(m_columns, s_rows);
	    const int tm_diagonals = MIN(m_columns, s_columns);

	    /*
	     * Add the non-zero Ts term.
	     */
	    if (eq_row < ts_diagonals) {
		const int s_cell = eq_row * s_columns + eq_column;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp != vnp->vn_zero) {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + eq_row,
				/*negative*/false, /*m*/-1, s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += ts_diagonals;

	    /*
	     * Add the Ti term.
	     */
	    if (eq_row < ti_diagonals && eq_row == eq_column) {
		if (add_term(function, vnmp, &vncpp_anchor,
			    base_coefficient + eq_row,
			    /*negative*/false, /*m*/-1, /*s*/-1) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += ti_diagonals;

	    /*
	     * Add the non-zero Tx terms.
	     */
	    for (int tx_d = 0; tx_d < tx_diagonals; ++tx_d) {
		const int m_cell = eq_row * m_columns + tx_d;
		const int s_cell = tx_d * s_columns + eq_column;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp != vnp->vn_zero) {
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + tx_d,
				/*negative*/true,
				/*m*/m_cell, /*s*/s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += tx_diagonals;

	    /*
	     * Add the Tm term.
	     */
	    if (eq_column < tm_diagonals) {
		const int m_cell = eq_row * m_columns + eq_column;

		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (eq_column == 0) {	/* let tm11 = 1.0 */
		    if (add_term(function, vnmp, &vncpp_anchor, -1,
				/*negative*/false,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		} else {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + eq_column - 1,
				/*negative*/true,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += tm_diagonals - 1;
	}
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
	{
	    const int um_diagonals = MIN(s_rows, m_rows);
	    const int ui_diagonals = MIN(s_rows, m_columns);
	    const int ux_diagonals = MIN(s_columns, m_rows);
	    const int us_diagonals = MIN(s_columns, m_columns);

	    /*
	     * Add the Um term.
	     */
	    if (eq_row < um_diagonals) {
		const int m_cell = eq_row * m_columns + eq_column;

		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (eq_row == 0) { /* let um11 = 1.0 */
		    if (add_term(function, vnmp, &vncpp_anchor, -1,
				/*negative*/true,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		} else {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + eq_row - 1,
				/*negative*/false,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += um_diagonals - 1;

	    /*
	     * Add the Ui term.
	     */
	    if (eq_row < ui_diagonals && eq_row == eq_column) {
		if (add_term(function, vnmp, &vncpp_anchor,
			    base_coefficient + eq_row,
			    /*negative*/false, /*m*/-1, /*s*/-1) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += ui_diagonals;

	    /*
	     * Add the non-zero Ux terms.
	     */
	    for (int ux_d = 0; ux_d < ux_diagonals; ++ux_d) {
		const int m_cell = ux_d * m_columns + eq_column;
		const int s_cell = eq_row * s_columns + ux_d;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp != vnp->vn_zero) {
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + ux_d,
				/*negative*/true,
				/*m*/m_cell, /*s*/s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += ux_diagonals;

	    /*
	     * Add the non-zero Us term.
	     */
	    if (eq_column < us_diagonals) {
		const int s_cell = eq_row * s_columns + eq_column;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp != vnp->vn_zero) {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + eq_column,
				/*negative*/true, /*m*/-1, s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += us_diagonals;
	}
	break;

    case VNACAL_T16:
	{
	    const int ts_rows	 = m_rows;
	    const int ts_columns = s_rows;
	    const int ti_rows	 = m_rows;
	    const int ti_columns = s_columns;
	    const int tx_rows	 = m_columns;
	    const int tx_columns = s_rows;
	    const int tm_rows	 = m_columns;
	    const int tm_columns = s_columns;

	    /*
	     * Add the non-zero Ts terms.
	     */
	    for (int ts_column = 0; ts_column < ts_columns; ++ts_column) {
		const int ts_row = eq_row;
		const int ts_cell = ts_row * ts_columns + ts_column;
		const int s_cell = ts_column * s_columns + eq_column;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp == vnp->vn_zero) {
		    continue;
		}
		if (add_term(function, vnmp, &vncpp_anchor,
			    base_coefficient + ts_cell,
			    /*negative*/false, /*m*/-1, s_cell) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += ts_rows * ts_columns;

	    /*
	     * Add the Ti term.
	     */
	    if (eq_row < ti_rows && eq_column < ti_columns) {
		const int ti_row     = eq_row;
		const int ti_column  = eq_column;
		const int ti_cell = ti_row * ti_columns + ti_column;

		if (add_term(function, vnmp, &vncpp_anchor,
			    base_coefficient + ti_cell,
			    /*negative*/false, /*m*/-1, /*s*/-1) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += ti_rows * ti_columns;

	    /*
	     * Add the non-zero Tx terms.
	     */
	    for (int tx_row = 0; tx_row < tx_rows; ++tx_row) {
		for (int tx_column = 0; tx_column < tx_columns; ++tx_column) {
		    const int tx_cell = tx_row * tx_columns + tx_column;
		    const int m_cell = eq_row * m_columns + tx_row;
		    const int s_cell = tx_column * s_columns + eq_column;
		    vnacal_new_parameter_t *vnprp =
			vnmp->vnm_s_matrix[s_cell];

		    assert(vnprp != NULL);
		    if (vnprp == vnp->vn_zero) {
			continue;
		    }
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + tx_cell,
				/*negative*/true,
				/*m*/m_cell, /*s*/s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += tx_rows * tx_columns;

	    /*
	     * Add the Tm terms.
	     */
	    for (int tm_row = 0; tm_row < tm_rows; ++tm_row) {
		const int tm_column = eq_column;
		const int tm_cell = tm_row * tm_columns + tm_column;
		const int m_cell = eq_row * m_columns + tm_row;

		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (tm_cell == 0) {	/* let tm11 = 1.0 */
		    if (add_term(function, vnmp, &vncpp_anchor, -1,
				/*negative*/false,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		} else {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + tm_cell - 1,
				/*negative*/true, /*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += tm_rows * tm_columns - 1;
	}
	break;

    case VNACAL_U16:
	{
	    const int um_rows	 = s_rows;
	    const int um_columns = m_rows;
	    const int ui_rows	 = s_rows;
	    const int ui_columns = m_columns;
	    const int ux_rows	 = s_columns;
	    const int ux_columns = m_rows;
	    const int us_rows	 = s_columns;
	    const int us_columns = m_columns;

	    /*
	     * Add the Um terms.
	     */
	    for (int um_column = 0; um_column < um_columns; ++um_column) {
		const int um_row = eq_row;
		const int um_cell = um_row * um_columns + um_column;
		const int m_cell = um_column * m_columns + eq_column;

		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (um_cell == 0) {	/* let um11 = 1.0 */
		    if (add_term(function, vnmp, &vncpp_anchor, -1,
				/*negative*/true, m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		} else {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + um_cell - 1,
				/*negative*/false, m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += um_rows * um_columns - 1;

	    /*
	     * Add the Ui term.
	     */
	    if (eq_row < ui_rows && eq_column < ui_columns) {
		const int ui_row     = eq_row;
		const int ui_column  = eq_column;
		const int ui_cell = ui_row * ui_columns + ui_column;

		if (add_term(function, vnmp, &vncpp_anchor,
			    base_coefficient + ui_cell,
			    /*negative*/false, /*m*/-1, /*s*/-1) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += ui_rows * ui_columns;

	    /*
	     * Add the non-zero Ux terms.
	     */
	    for (int ux_row = 0; ux_row < ux_rows; ++ux_row) {
		for (int ux_column = 0; ux_column < ux_columns; ++ux_column) {
		    const int ux_cell = ux_row * ux_columns + ux_column;
		    const int m_cell = ux_column * m_columns + eq_column;
		    const int s_cell = eq_row * s_columns + ux_row;
		    vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		    assert(vnprp != NULL);
		    if (vnprp == vnp->vn_zero) {
			continue;
		    }
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + ux_cell,
				/*negative*/true,
				/*m*/m_cell, /*s*/s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += ux_rows * ux_columns;

	    /*
	     * Add the non-zero Us terms.
	     */
	    for (int us_row = 0; us_row < us_rows; ++us_row) {
		const int us_column = eq_column;
		const int us_cell = us_row * us_columns + us_column;
		const int s_cell = eq_row * s_columns + us_row;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		if (vnprp == vnp->vn_zero) {
		    continue;
		}
		if (add_term(function, vnmp, &vncpp_anchor,
			    base_coefficient + us_cell,
			    /*negative*/true, /*m*/-1, /*s*/s_cell) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += us_rows * us_columns;
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	{
	    const int um_diagonals = MIN(s_rows, m_rows);
	    const int ui_diagonals = 1;
	    const int ux_diagonals = MIN(s_columns, m_rows);
	    const int us_diagonals = 1;

	    /*
	     * Add the Um term.
	     */
	    if (eq_row < um_diagonals) {
		const int m_cell = eq_row * m_columns + eq_column;

		assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		if (eq_row == eq_column) {
		    if (add_term(function, vnmp, &vncpp_anchor, -1,
				/*negative*/true,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		} else {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + eq_row -
				(eq_row > eq_column),
				/*negative*/false,
				/*m*/m_cell, /*s*/-1) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += um_diagonals - 1;

	    /*
	     * Add the Ui term.
	     */
	    if (eq_row == eq_column) {
		if (add_term(function, vnmp, &vncpp_anchor, base_coefficient,
			    /*negative*/false, /*m*/-1, /*s*/-1) == -1) {
		    return -1;
		}
	    }
	    base_coefficient += ui_diagonals;

	    /*
	     * Add the non-zero Ux terms.
	     */
	    for (int ux_d = 0; ux_d < ux_diagonals; ++ux_d) {
		const int m_cell = ux_d * m_columns + eq_column;
		const int s_cell = eq_row * s_columns + ux_d;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp != vnp->vn_zero) {
		    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient + ux_d, /*negative*/true,
				/*m*/m_cell, /*s*/s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += ux_diagonals;

	    /*
	     * Add the non-zero Us term.
	     */
	    if (eq_column < s_columns) {
		const int s_cell = eq_row * s_columns + eq_column;
		vnacal_new_parameter_t *vnprp = vnmp->vnm_s_matrix[s_cell];

		assert(vnprp != NULL);
		if (vnprp != vnp->vn_zero) {
		    if (add_term(function, vnmp, &vncpp_anchor,
				base_coefficient,
				/*negative*/true, /*m*/-1, s_cell) == -1) {
			return -1;
		    }
		}
	    }
	    base_coefficient += us_diagonals;
	}
	break;

    case VNACAL_E12:
    default:
	abort();
    }
    return 0;
}

/*
 * _vnacal_new_add_common: common function to add calibration equations
 *   @vnaa: common argument structure
 */
int _vnacal_new_add_common(vnacal_new_add_arguments_t vnaa)
{
    /* make short aliases for the commonly used arguments */
    const char *const function = vnaa.vnaa_function;
    vnacal_new_t *const vnp = vnaa.vnaa_cmp;
    const double complex *const *const a_matrix = vnaa.vnaa_a_matrix;
    const int a_rows = vnaa.vnaa_a_rows;
    const int a_columns = vnaa.vnaa_a_columns;
    const double complex *const *const b_matrix = vnaa.vnaa_b_matrix;
    const int b_rows = vnaa.vnaa_b_rows;
    const int b_columns = vnaa.vnaa_b_columns;
    const int b_diagonals = MIN(b_rows, b_columns);
    const int b_cells = vnaa.vnaa_m_is_diagonal ?
	b_diagonals : b_rows * b_columns;
    const int *const s_matrix = vnaa.vnaa_s_matrix;
    const int s_rows = vnaa.vnaa_s_rows;
    const int s_columns = vnaa.vnaa_s_columns;
    const int s_diagonals = MIN(s_rows, s_columns);
    const int s_ports = MAX(s_rows, s_columns);
    const int s_cells = vnaa.vnaa_s_is_diagonal ?
	s_diagonals : s_rows * s_columns;
    const int *const s_port_map = vnaa.vnaa_s_port_map;
    vnacal_t *const vcp = vnp->vn_vcp;
    const vnacal_layout_t *const vlp = &vnp->vn_layout;
    const int full_m_rows = VL_M_ROWS(vlp);
    const int full_m_columns = VL_M_COLUMNS(vlp);
    const int full_s_rows = VL_S_ROWS(vlp);
    const int full_s_columns = VL_S_COLUMNS(vlp);
    const int full_s_ports = MAX(full_s_rows, full_s_columns);
    const int frequencies = vnp->vn_frequencies;

    /* parameter type: 'T' or 'U' */
    char ptype = '\000';

    /* mimimum allowed size of the b matrix */
    int min_b_rows, min_b_columns;

    /* map from b_matrix index to vnm_m_matrix index */
    int m_cell_map[MAX(1, MIN(b_cells, full_m_rows * full_m_columns))];

    /* map from s_matrix index to vnm_s_matrix index */
    int s_cell_map[s_cells];

    /* which VNA ports are connected to the standard */
    bool port_connected[full_s_ports];

    /* which rows and columns in vnm_m_matrix and vnm_s_matrix were given */
    bool m_row_given[full_m_rows];
    bool m_column_given[full_m_columns];
    bool s_row_given[full_s_rows];
    bool s_column_given[full_s_columns];

    /* new measured standard we're adding */
    vnacal_new_measurement_t *vnmp = NULL;

    /* alias for vnm_m_matrix */
    double complex **full_m_matrix;

    /* alias for vnm_s_matrix */
    vnacal_new_parameter_t **full_s_matrix;

    /* list of new equations */
    vnacal_new_equation_t *ncep_head = NULL;

    /* address where next equation should be linked */
    vnacal_new_equation_t **vnepp_anchor = &ncep_head;

    /* return code */
    int rc = -1;

    /*
     * Validate that b is not NULL.  If it is, make the error messages
     * reflect the name of the actual parameter.
     */
    if (b_matrix == NULL) {
	if (vnaa.vnaa_m_type == 'a') {
	    _vnacal_error(vcp, VNAERR_USAGE, "%s: NULL m matrix", function);
	} else {
	    _vnacal_error(vcp, VNAERR_USAGE, "%s: NULL b matrix", function);
	}
	goto out;
    }

    /*
     * Collect per parameter type information used below.
     */
    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	ptype = 'T';
	min_b_rows    = s_ports;
	min_b_columns = s_ports;
	break;

    case VNACAL_T16:
	ptype = 'T';
	min_b_rows    = s_rows;
	min_b_columns = full_m_columns;
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	ptype = 'U';
	min_b_rows    = s_ports;
	min_b_columns = s_ports;
	break;

    case VNACAL_U16:
	ptype = 'U';
	min_b_rows    = full_m_rows;
	min_b_columns = s_columns;
	break;

    case VNACAL_E12:
    default:
	abort();
    }
    assert(ptype != '\000');

    /*
     * Check the S matrix size.
     */
    if (s_rows < 1 || s_rows > full_s_rows) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid s_rows value: %d",
		function, s_rows);
	goto out;
    }
    if (s_columns < 1 || s_columns > full_s_columns) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid s_columns value: %d",
		function, s_columns);
	goto out;
    }

    /*
     * When a rectangular S matrix is given, it means that we don't
     * fully know the S parameters of the standard.  In T parameters,
     * the entire S column must be known and in U parameters, the entire
     * S row must be known.  Make sure the shape of the S parameter
     * matrix is consistent with the error term type so that no required
     * parameters are unknown.
     */
    if (s_rows < s_columns && s_rows != full_s_rows && ptype == 'T') {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: s_rows cannot be less than %d",
		function, MIN(s_columns, full_s_rows));
	goto out;
    }
    if (s_rows > s_columns && s_columns != full_s_columns && ptype == 'U') {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: s_columns cannot be less than %d",
		function, MIN(s_rows, full_s_columns));
	goto out;
    }
    assert(!vnaa.vnaa_s_is_diagonal || s_rows == s_columns);

    /*
     * Make sure a port map was provided if the S matrix is smaller than
     * the calibration matrix.
     */
    if (s_port_map == NULL &&
	    (s_rows != full_s_rows || s_columns != full_s_ports)) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: port map is required when the given S "
		"matrix is smaller than that of the calibration", function);
	goto out;
    }

    /*
     * Validate the dimensions of the B matrix.
     */
    if (b_rows != min_b_rows && b_rows != full_m_rows) {
	if (min_b_rows == full_m_rows) {
	    if (vnaa.vnaa_m_type == 'a') {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: b_rows must be %d",
			function, full_m_rows);
	    } else {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: m_rows must be %d",
			function, full_m_rows);
	    }
	} else {
	    if (vnaa.vnaa_m_type == 'a') {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: b_rows must be %d or %d",
			function, min_b_rows, full_m_rows);
	    } else {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: m_rows must be %d or %d",
			function, min_b_rows, full_m_rows);
	    }
	}
	goto out;
    }
    if (b_columns != min_b_columns && b_columns != full_m_columns) {
	if (min_b_columns == full_m_columns) {
	    if (vnaa.vnaa_m_type == 'a') {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: b_columns must be %d",
			function, full_m_columns);
	    } else {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: m_columns must be %d",
			function, full_m_columns);
	    }
	} else {
	    if (vnaa.vnaa_m_type == 'a') {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: b_columns must be %d or %d",
			function, min_b_columns, full_m_columns);
	    } else {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: m_columns must be %d or %d",
			function, min_b_columns, full_m_columns);
	    }
	}
	goto out;
    }

    /*
     * If an A matrix was given, validate its dimensions.  Normally, it
     * must be square and of the same dimensions as b_column.  In UE14,
     * however, it's a row-vector of 1x1 a matrices.
     */
    if (a_matrix != NULL) {
	const int rows = VL_IS_UE14(vlp) ? 1 : b_columns;

	if (a_rows != rows || a_columns != b_columns) {
	    _vnacal_error(vcp, VNAERR_USAGE, "%s: 'a' matrix must be %d x %d",
		    function, rows, b_columns);
	    goto out;
	}
    }

    /*
     * If a port map was given, check for out-of-bounds and duplicate
     * port indices.  At the same time, initialize the port_connected array
     * to indicate which s-ports appear in s_port_map.  Set all elements
     * to true if no port map was given.
     */
    if (s_port_map != NULL) {
	int max_port = 0;

	(void)memset((void *)&port_connected, 0, sizeof(port_connected));
	for (int s_port_index = 0; s_port_index < s_ports; ++s_port_index) {
	    int port = s_port_map[s_port_index];

	    if (port < 1) {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: %d: invalid port index",
			function, port);
		goto out;
	    }
	    if (port > max_port) {
		max_port = port;
	    }
	    if (s_port_index < s_rows && max_port > full_s_rows) {
		_vnacal_error(vcp, VNAERR_USAGE, "%s: port index %d exceeds "
			"calibration matrix row bound", function, max_port);
		goto out;
	    }
	    if (s_port_index < s_columns && max_port > full_s_columns) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: s_port_index index %d exceeds "
			"calibration matrix column bound", function, max_port);
		goto out;
	    }
	    if (port_connected[port - 1]) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: s_port_index index %d arg "
			"appears more than once", function, port);
		goto out;
	    }
	    port_connected[port - 1] = true;
	}
    } else {
	for (int port_index = 0; port_index < full_s_ports; ++port_index) {
	    port_connected[port_index] = true;
	}
    }

    /*
     * Create maps between the cells of the B and S matrices given in
     * the argument structure to the cells of the M and S matrices in the
     * vnacal_new_measurement_t structure, taking the diagonal cases and port
     * maps into account.  These maps are used to populate vnm_m_matrix
     * and vnm_s_matrix below.	At the same time, construct vectors
     * of flags indicating exactly which rows and columns were given.
     * These are used below to decide which equations to generate.
     */
    (void)memset((void *)m_row_given, 0, full_m_rows * sizeof(bool));
    (void)memset((void *)m_column_given, 0, full_m_columns * sizeof(bool));
    (void)memset((void *)s_row_given, 0, full_s_rows * sizeof(bool));
    (void)memset((void *)s_column_given, 0, full_s_columns * sizeof(bool));
    if (s_port_map != NULL) {
	int m_port_map[s_ports];

	/*
	 * Make a sorted version of the port map for M.	 The rows and
	 * columns of the M matrix remain in relative order even when the
	 * port map is used to reorder rows and columns of the S matrix.
	 */
	(void)memcpy((void *)m_port_map, (void *)s_port_map,
		s_ports * sizeof(int));
	qsort((void *)m_port_map, s_ports, sizeof(int), int_cmp);

	/*
	 * Make a map from the cells of the argument B matrix to the cells
	 * of the vnacal_new_measurement_t M matrix in (with port map).
	 */
	if (vnaa.vnaa_m_is_diagonal) {
	    for (int b_diagonal = 0; b_diagonal < b_diagonals; ++b_diagonal) {
		int full_m_row = (b_rows < full_m_rows) ?
		    m_port_map[b_diagonal] - 1: b_diagonal;
		int full_m_column = (b_columns < full_m_columns) ?
		    m_port_map[b_diagonal] - 1: b_diagonal;

		m_cell_map[b_diagonal] =
		    full_m_row * full_s_columns + full_m_column;
		m_row_given[full_m_row] = true;
		m_column_given[full_m_column] = true;
	    }
	} else {
	    for (int b_row = 0; b_row < b_rows; ++b_row) {
		int full_m_row = (b_rows < full_m_rows) ?
		    m_port_map[b_row] - 1: b_row;

		m_row_given[full_m_row] = true;
		for (int b_column = 0; b_column < b_columns; ++b_column) {
		    int full_m_column = (b_columns < full_m_columns) ?
			m_port_map[b_column] - 1: b_column;

		    m_cell_map[b_row * b_columns + b_column] =
			full_m_row * full_m_columns + full_m_column;
		    m_column_given[full_m_column] = true;
		}
	    }
	}

	/*
	 * Make a map from the cells of the argument strucure S matrix
	 * to the cells of the vnacal_new_measurement_t S matrix (with
	 * port map).
	 */
	if (vnaa.vnaa_s_is_diagonal) {
	    for (int s_diagonal = 0; s_diagonal < s_diagonals; ++s_diagonal) {
		const int full_diagonal = s_port_map[s_diagonal] - 1;

		s_cell_map[s_diagonal] = full_diagonal * (full_s_columns + 1);
		s_row_given[full_diagonal] = true;
		s_column_given[full_diagonal] = true;
	    }
	} else {
	    for (int s_row = 0; s_row < s_rows; ++s_row) {
		int full_s_row = s_port_map[s_row] - 1;

		s_row_given[full_s_row] = true;
		for (int s_column = 0; s_column < s_columns; ++s_column) {
		    int full_s_column = s_port_map[s_column] - 1;

		    s_cell_map[s_row * s_columns + s_column] =
			full_s_row * full_s_columns + full_s_column;
		    s_column_given[full_s_column] = true;
		}
	    }
	}
    } else {
	/*
	 * Make a map from the cells of the argument B matrix to the
	 * cells of the vnacal_new_measurement_t M matrix in (no port map).
	 * This is not a simple 1:1 map if the b matrix contains only the
	 * diagonal vector.
	 */
	if (vnaa.vnaa_m_is_diagonal) {
	    for (int b_diagonal = 0; b_diagonal < b_diagonals; ++b_diagonal) {
		m_cell_map[b_diagonal] = b_diagonal * (full_m_columns + 1);
		m_row_given[b_diagonal] = true;
		m_column_given[b_diagonal] = true;
	    }
	} else {
	    for (int b_cell = 0; b_cell < b_cells; ++b_cell) {
		m_cell_map[b_cell] = b_cell;
	    }
	    for (int b_row = 0; b_row < b_rows; ++b_row) {
		m_row_given[b_row] = true;
	    }
	    for (int b_column = 0; b_column < b_rows; ++b_column) {
		m_column_given[b_column] = true;
	    }
	}

	/*
	 * Make a map from the cells of the argument strucure S matrix
	 * to the cells of the vnacal_new_measurement_t S matrix (no
	 * port map).
	 */
	if (vnaa.vnaa_s_is_diagonal) {
	    for (int s_diagonal = 0; s_diagonal < s_diagonals; ++s_diagonal) {
		s_cell_map[s_diagonal] = s_diagonal * (full_s_columns + 1);
		s_row_given[s_diagonal] = true;
		s_column_given[s_diagonal] = true;
	    }
	} else {
	    for (int s_cell = 0; s_cell < s_cells; ++s_cell) {
		s_cell_map[s_cell] = s_cell;
	    }
	    for (int s_row = 0; s_row < s_rows; ++s_row) {
		s_row_given[s_row] = true;
	    }
	    for (int s_column = 0; s_column < s_columns; ++s_column) {
		s_column_given[s_column] = true;
	    }
	}
    }

    /*
     * Allocate and init the vnacal_new_measurement_t structure and its
     * vectors of per-frequency M values.
     */
    if ((vnmp = malloc(sizeof(vnacal_new_measurement_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	goto out;
    }
    (void)memset((void *)vnmp, 0, sizeof(vnacal_new_measurement_t));
    if ((vnmp->vnm_m_matrix = full_m_matrix =
		calloc(full_m_rows * full_m_columns,
		    sizeof(double complex *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto out;
    }
    for (int b_cell = 0; b_cell < b_cells; ++b_cell) {
	int full_m_cell = m_cell_map[b_cell];

	if ((full_m_matrix[full_m_cell] = calloc(frequencies,
			sizeof(double complex))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    goto out;
	}
    }
    if ((vnmp->vnm_s_matrix = full_s_matrix =
		calloc(full_s_rows * full_s_columns,
		    sizeof(vnacal_new_parameter_t *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc %s", strerror(errno));
	goto out;
    }
    vnmp->vnm_ncp = vnp;

    /*
     * If no 'a' matrix was given, just copy the m vectors.
     */
    if (a_matrix == NULL) {
	for (int b_cell = 0; b_cell < b_cells; ++b_cell) {
	    int full_m_cell = m_cell_map[b_cell];

	    (void)memcpy((void *)&full_m_matrix[full_m_cell][0],
		    (void *)&b_matrix[b_cell][0],
		    frequencies * sizeof(double complex));
	}
    }
    /*
     * Else if UE14, the 'a' matrix is a row vector of 1x1 matrices.
     * Divide each column by its corresponding 'a' entry.
     */
    else if (VL_IS_UE14(vlp)) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int b_column = 0; b_column < b_columns; ++b_column) {
		double complex a = a_matrix[b_column][findex];

		if (a == 0.0) {
		    _vnacal_error(vcp, VNAERR_MATH,
			    "%s: 'a' matrix is singular at frequency index %d",
			    function, findex);
		    goto out;
		}
		for (int b_row = 0; b_row < b_rows; ++b_row) {
		    int b_cell = b_row * b_columns + b_column;
		    int full_m_cell = m_cell_map[b_cell];

		    full_m_matrix[full_m_cell][findex] =
			b_matrix[b_cell][findex] / a;
		}
	    }
	}
    }
    /*
     * Else, find M = B * A^-1 for each frequency.
     */
    else {
	assert(!vnaa.vnaa_m_is_diagonal);
	for (int findex = 0; findex < frequencies; ++findex) {
	    double complex a[a_rows][a_columns];
	    double complex b[b_rows][b_columns];
	    double complex m[b_rows][b_columns];
	    double complex determinant;

	    assert(a_rows == a_columns);
	    assert(a_rows == b_columns);
	    for (int a_row = 0; a_row < a_rows; ++a_row) {
		for (int a_column = 0; a_column < a_columns; ++a_column) {
		    int a_cell = a_row * a_columns + a_column;

		    a[a_row][a_column] = a_matrix[a_cell][findex];
		}
	    }
	    for (int b_row = 0; b_row < b_rows; ++b_row) {
		for (int b_column = 0; b_column < b_columns; ++b_column) {
		    int b_cell = b_row * b_columns + b_column;

		    b[b_row][b_column] = b_matrix[b_cell][findex];
		}
	    }
	    determinant = _vnacommon_mrdivide(&m[0][0], &b[0][0], &a[0][0],
		    b_rows, b_columns);
	    if (determinant == 0.0) {
		_vnacal_error(vcp, VNAERR_MATH,
			"%s: 'a' matrix is singular at frequency index %d",
			function, findex);
		goto out;
	    }
	    for (int m_row = 0; m_row < b_rows; ++m_row) {
		for (int m_column = 0; m_column < b_columns; ++m_column) {
		    int m_cell = m_row * b_columns + m_column;
		    int full_m_cell = m_cell_map[m_cell];

		    full_m_matrix[full_m_cell][findex] = m[m_row][m_column];
		}
	    }
	}
    }

    /*
     * Construct the vnacal_new_measurement_t S matrix.
     */
    if ((vnmp->vnm_s_matrix = full_s_matrix =
		calloc(full_s_rows * full_s_columns,
		    sizeof(vnacal_new_parameter_t *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto out;
    }
    for (int s_cell = 0; s_cell < s_cells; ++s_cell) {
	if ((full_s_matrix[s_cell_map[s_cell]] =
		    _vnacal_new_get_parameter(function, vnp,
			s_matrix[s_cell])) == NULL) {
	    goto out;
	}
    }

    /*
     * If the given S matrix is diagonal, fill in the off-diagonal
     * entries with zeros.
     */
    if (vnaa.vnaa_s_is_diagonal) {
	for (int r = 0; r < full_s_rows; ++r) {
	    for (int c = 0; c < full_s_columns; ++c) {
		if (r != c && port_connected[r] && port_connected[c]) {
		    const int cell = r * full_s_columns + c;

		    assert(full_s_matrix[cell] == NULL);
		    full_s_matrix[cell] = vnp->vn_zero;
		}
	    }
	}
    }

    /*
     * When the calibration standard connects to only a subset of the VNA
     * ports, we don't presume to know anything about the S parameters
     * of the remaining ports except that they have no connection to
     * the connected ports.  For example, suppose vnm_s_matrix is 5x5
     * and the given s_matrix is 2x3 with a map containing { 1, 2, 3 }
     * zero based, i.e. offset one from vnm_s_matrix.  Using the notation
     * small s11, s12, etc.  for the elements of s_matrix, and capital
     * S11, S12, etc., both one-based, for the elements of the full
     * matrix, we have:
     *
     *		*   0	0   0	*
     *		0   s11 s12 ?	0
     *		0   s21 s22 ?	0
     *		0   s31 s32 ?	0
     *		*   0	0   0	*
     *
     * Where 0's indicate cells we know to be zero.  We know nothing
     * about S11, S15, S51 or S55.  Similarly, S24, S34 and S35 (the
     * ones with question marks) are unknown because they weren't given
     * in s_matrix.  But we know that there are no external connections
     * between the 1,5 group and the 2,3,4 group.  Reflect the cells
     * known to be zero in vnm_s_matrix.
     */
    if (s_port_map != NULL) {
	for (int r = 0; r < full_s_rows; ++r) {
	    for (int c = 0; c < full_s_columns; ++c) {
		if ((port_connected[r] && !port_connected[c]) ||
		    (!port_connected[r] && port_connected[c])) {
		    const int cell = r * full_s_columns + c;

		    assert(full_s_matrix[cell] == NULL);
		    full_s_matrix[cell] = vnp->vn_zero;
		}
	    }
	}
    }

    /*
     * For all calibration types except T16 and U16 that handle leakage
     * terms within the linear system, find the transitive closure
     * of the S parameter matrix to determine which port pairs have a
     * signal path through the calibration standard.  What we're really
     * interested in here is which pairs of ports -don't- have a signal
     * path between them.  If we can't prove that there isn't a signal
     * path, we assume there is one.  For example, in a 3x3 system with
     * a short standard on port 1, we know that S11 is -1, S12, S13,
     * S21, and S31 are 0, but we don't know S22, S23, S32 or S33, and
     * we assume they're non-zero.
     *
     * When a vector standard is given, we -could- check if the value at
     * a given frequency is exactly zero and include that case, but doing
     * so would require a separate reachability graph for each frequency,
     * and when one supplies a vector standard, it's unlikely that any
     * value is exactly zero, so we ignore that case.  For finding the
     * transitive closure, we use the Floyd Warshall algorithm.
     */
    switch (VL_TYPE(vlp)) {
	bool *matrix;

    case VNACAL_T8:
    case VNACAL_U8:
    case VNACAL_TE10:
    case VNACAL_UE10:
    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	if ((vnmp->vnm_reachability_matrix = matrix = calloc(full_s_rows *
			full_s_columns, sizeof(bool))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    goto out;
	}
	for (int r = 0; r < full_s_rows; ++r) {
	    for (int c = 0; c < full_s_columns; ++c) {
		int cell = r * full_s_columns + c;

		matrix[cell] = full_s_matrix[cell] != vnp->vn_zero;
	    }
	}
	for (int i = 0; i < MIN(full_s_rows, full_s_columns); ++i) {
	    for (int r = 0; r < full_s_rows; ++r) {
		const int ri_cell = r * full_s_columns + i;

		for (int c = 0; c < full_s_columns; ++c) {
		    const int rc_cell = r * full_s_columns + c;
		    const int ic_cell = i * full_s_columns + c;

		    matrix[rc_cell] |= matrix[ri_cell] && matrix[ic_cell];
		}
	    }
	}
	break;

    case VNACAL_T16:
    case VNACAL_U16:
	break;

    case VNACAL_E12:
    default:
	abort();
    }

    /*
     * Generate the equations and add them to a temporary list.  In T
     * parameters, there are at most m_rows x s_columns equations, and
     * in U parameters, there are at most s_rows x m_columns equations.
     * We can only create equations for the rows and columns that the
     * user actually gave us, however.
     */
    if (VL_IS_UE14(vlp)) {
	for (int eq_column = 0; eq_column < full_m_columns; ++eq_column) {
	    for (int eq_row = 0; eq_row < full_s_rows; ++eq_row) {
		if (s_row_given[eq_row] && m_column_given[eq_column]) {
		    rc = add_equation(function, vnmp, &vnepp_anchor,
			    eq_row, eq_column);
		    if (rc == -1) {
			goto out;
		    }
		}
	    }
	}
    } else if (ptype == 'T') {
	for (int eq_row = 0; eq_row < full_m_rows; ++eq_row) {
	    for (int eq_column = 0; eq_column < full_s_columns; ++eq_column) {
		if (m_row_given[eq_row] && s_column_given[eq_column]) {
		    rc = add_equation(function, vnmp, &vnepp_anchor,
			    eq_row, eq_column);
		    if (rc == -1) {
			goto out;
		    }
		}
	    }
	}
    } else {
	for (int eq_row = 0; eq_row < full_s_rows; ++eq_row) {
	    for (int eq_column = 0; eq_column < full_m_columns; ++eq_column) {
		if (s_row_given[eq_row] && m_column_given[eq_column]) {
		    rc = add_equation(function, vnmp, &vnepp_anchor,
			    eq_row, eq_column);
		    if (rc == -1) {
			goto out;
		    }
		}
	    }
	}
    }

    /*
     * Link the new standard onto the vnacal_new_t.
     */
    *vnp->vn_measurement_anchor = vnmp;
    vnp->vn_measurement_anchor = &vnmp->vnm_next;
    vnmp = NULL; /* no longer ours to free */

    /*
     * Link the equations onto their respective systems.
     */
    while (ncep_head != NULL) {
	vnacal_new_equation_t *vnep = ncep_head;
	const int system = VL_IS_UE14(vlp) ? vnep->vne_column : 0;
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[system];

	ncep_head = vnep->vne_next;
	vnep->vne_next = NULL;
	*vnsp->vns_equation_anchor = vnep;
	vnsp->vns_equation_anchor = &vnep->vne_next;
	if (++vnsp->vns_equation_count > vnp->vn_max_equations) {
	    vnp->vn_max_equations = vnsp->vns_equation_count;
	}
	++vnp->vn_equations;
    }
    rc = 0;

out:
    /*
     * Clean up.
     */
    while (ncep_head != NULL) {
	vnacal_new_equation_t *vnep = ncep_head;

	ncep_head = vnep->vne_next;
	while (vnep->vne_coefficient_list != NULL) {
	    vnacal_new_coefficient_t *vncp = vnep->vne_coefficient_list;

	    vnep->vne_coefficient_list = vncp->vnc_next;
	    free((void *)vncp);
	}
	free((void *)vnep);
    }
    _vnacal_new_free_measurement(vnmp);

    return rc;
}

/*
 * vnacal_new_add_single_reflect: add a single reflect on the given port
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s11: s11 parameter index, e.g. VNACAL_MATCH
 *   @port: VNA port on which the measurement is made
 */
int vnacal_new_add_single_reflect(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	int s11, int port)
{
    vnacal_new_add_arguments_t vnaa;

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= (const double complex *const *)a;
    vnaa.vnaa_a_rows		= a_rows;
    vnaa.vnaa_a_columns		= a_columns;
    vnaa.vnaa_b_matrix		= (const double complex *const *)b;
    vnaa.vnaa_b_rows		= b_rows;
    vnaa.vnaa_b_columns		= b_columns;
    vnaa.vnaa_s_matrix		= &s11;
    vnaa.vnaa_s_rows		= 1;
    vnaa.vnaa_s_columns		= 1;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= true;
    vnaa.vnaa_m_type		= 'a';
    vnaa.vnaa_s_port_map	= &port;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_single_reflect_m: add a single reflect on the given port
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s11: s11 parameter index, e.g. VNACAL_MATCH
 *   @port: VNA port on which the measurement is made
 */
int vnacal_new_add_single_reflect_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	int s11, int port)
{
    vnacal_new_add_arguments_t vnaa;

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= NULL;
    vnaa.vnaa_a_rows		= 0;
    vnaa.vnaa_a_columns		= 0;
    vnaa.vnaa_b_matrix		= (const double complex *const *)m;
    vnaa.vnaa_b_rows		= m_rows;
    vnaa.vnaa_b_columns		= m_columns;
    vnaa.vnaa_s_matrix		= &s11;
    vnaa.vnaa_s_rows		= 1;
    vnaa.vnaa_s_columns		= 1;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= true;
    vnaa.vnaa_m_type		= 'm';
    vnaa.vnaa_s_port_map	= &port;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_double_reflect: add a pair of reflects
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s11: s11 parameter index
 *   @s22: s22 parameter index
 *   @port1: first reflect is on this VNA port
 *   @port2: second reflect is on this VNA port
 */
int vnacal_new_add_double_reflect(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	int s11, int s22, int port1, int port2)
{
    vnacal_new_add_arguments_t vnaa;
    int s_vector[2] = { s11, s22 };
    int port_map[2] = { port1, port2 };

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= (const double complex *const *)a;
    vnaa.vnaa_a_rows		= a_rows;
    vnaa.vnaa_a_columns		= a_columns;
    vnaa.vnaa_b_matrix		= (const double complex *const *)b;
    vnaa.vnaa_b_rows		= b_rows;
    vnaa.vnaa_b_columns		= b_columns;
    vnaa.vnaa_s_matrix		= s_vector;
    vnaa.vnaa_s_rows		= 2;
    vnaa.vnaa_s_columns		= 2;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= true;
    vnaa.vnaa_m_type		= 'a';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_double_reflect_m: add a pair of reflects
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s11: s11 parameter index, e.g. VNACAL_MATCH
 *   @port1: first reflect is on this VNA port
 *   @port2: second reflect is on this VNA port
 */
int vnacal_new_add_double_reflect_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	int s11, int s22, int port1, int port2)
{
    vnacal_new_add_arguments_t vnaa;
    int s_vector[2] = { s11, s22 };
    int port_map[2] = { port1, port2 };

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= NULL;
    vnaa.vnaa_a_rows		= 0;
    vnaa.vnaa_a_columns		= 0;
    vnaa.vnaa_b_matrix		= (const double complex *const *)m;
    vnaa.vnaa_b_rows		= m_rows;
    vnaa.vnaa_b_columns		= m_columns;
    vnaa.vnaa_s_matrix		= s_vector;
    vnaa.vnaa_s_rows		= 2;
    vnaa.vnaa_s_columns		= 2;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= true;
    vnaa.vnaa_m_type		= 'm';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * TODO: In the T8 and U8 reflect cases, we need only the M diagonal.
 * Consider adding an alternative add_double_reflect to the API that takes
 * just m11 and m22 vectors intead of the full M matrix.  Note that in
 * AB, we still need the full A, so this shortcut can be used only in
 * The T8/U8 M cases.
 *
 * int vnacal_new_add_single_reflect_m8(vnacal_new_t *vnp,
 *   const double complex *m11, int s11, int port)
 *
 * int vnacal_new_add_double_reflect_m8(vnacal_new_t *vnp,
 *   const double complex *m11, const double complex *m22,
 *   int s11, int s22, int port1, int port2).
 */

/*
 * vnacal_new_add_line: add an arbitrary two-port-standard
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s_2x2: 2x2 matrix of parameter indices (by rows) describing the standard
 *   @port1: first VNA port attached to standard
 *   @port2: second VNA port attached to standard
 */
int vnacal_new_add_line(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	const int *s_2x2, int port1, int port2)
{
    vnacal_new_add_arguments_t vnaa;
    int port_map[2] = { port1, port2 };

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= (const double complex *const *)a;
    vnaa.vnaa_a_rows		= a_rows;
    vnaa.vnaa_a_columns		= a_columns;
    vnaa.vnaa_b_matrix		= (const double complex *const *)b;
    vnaa.vnaa_b_rows		= b_rows;
    vnaa.vnaa_b_columns		= b_columns;
    vnaa.vnaa_s_matrix		= s_2x2;
    vnaa.vnaa_s_rows		= 2;
    vnaa.vnaa_s_columns		= 2;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= false;
    vnaa.vnaa_m_type		= 'a';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_line_m: add an arbitrary two-port-standard
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s_2x2: 2x2 matrix of parameter indices (by rows) describing the standard
 *   @port1: first VNA port attached to standard
 *   @port2: second VNA port attached to standard
 */
int vnacal_new_add_line_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	const int *s_2x2, int port1, int port2)
{
    vnacal_new_add_arguments_t vnaa;
    int port_map[2] = { port1, port2 };

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= NULL;
    vnaa.vnaa_a_rows		= 0;
    vnaa.vnaa_a_columns		= 0;
    vnaa.vnaa_b_matrix		= (const double complex *const *)m;
    vnaa.vnaa_b_rows		= m_rows;
    vnaa.vnaa_b_columns		= m_columns;
    vnaa.vnaa_s_matrix		= s_2x2;
    vnaa.vnaa_s_rows		= 2;
    vnaa.vnaa_s_columns		= 2;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= false;
    vnaa.vnaa_m_type		= 'm';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_through: add a perfect through between two VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @port1: first VNA port attached to through
 *   @port2: second VNA port attached to through
 */
int vnacal_new_add_through(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	int port1, int port2)
{
    vnacal_new_add_arguments_t vnaa;
    int s_2x2[2][2] = {{ VNACAL_ZERO, VNACAL_ONE },
		       { VNACAL_ONE, VNACAL_ZERO }};
    int port_map[2] = { port1, port2 };

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= (const double complex *const *)a;
    vnaa.vnaa_a_rows		= a_rows;
    vnaa.vnaa_a_columns		= a_columns;
    vnaa.vnaa_b_matrix		= (const double complex *const *)b;
    vnaa.vnaa_b_rows		= b_rows;
    vnaa.vnaa_b_columns		= b_columns;
    vnaa.vnaa_s_matrix		= &s_2x2[0][0];
    vnaa.vnaa_s_rows		= 2;
    vnaa.vnaa_s_columns		= 2;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= false;
    vnaa.vnaa_m_type		= 'a';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_through_m: add a perfect through between two VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @port1: first VNA port attached to through
 *   @port2: second VNA port attached to through
 */
int vnacal_new_add_through_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	int port1, int port2)
{
    vnacal_new_add_arguments_t vnaa;
    int s_2x2[2][2] = {{ VNACAL_ZERO, VNACAL_ONE },
		       { VNACAL_ONE, VNACAL_ZERO }};
    int port_map[2] = { port1, port2 };

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= NULL;
    vnaa.vnaa_a_rows		= 0;
    vnaa.vnaa_a_columns		= 0;
    vnaa.vnaa_b_matrix		= (const double complex *const *)m;
    vnaa.vnaa_b_rows		= m_rows;
    vnaa.vnaa_b_columns		= m_columns;
    vnaa.vnaa_s_matrix		= &s_2x2[0][0];
    vnaa.vnaa_s_rows		= 2;
    vnaa.vnaa_s_columns		= 2;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= false;
    vnaa.vnaa_m_type		= 'm';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_mapped_matrix: add a matrix of measurements with port map
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s: matrix of parameter indices (by rows)
 *   @s_rows: number of rows in s
 *   @s_columns: number of columns in s
 *   @port_map: vector of VNA port numbers corresponding to these ports
 */
int vnacal_new_add_mapped_matrix(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	const int *s, int s_rows, int s_columns,
	const int *port_map)
{
    vnacal_new_add_arguments_t vnaa;

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= (const double complex *const *)a;
    vnaa.vnaa_a_rows		= a_rows;
    vnaa.vnaa_a_columns		= a_columns;
    vnaa.vnaa_b_matrix		= (const double complex *const *)b;
    vnaa.vnaa_b_rows		= b_rows;
    vnaa.vnaa_b_columns		= b_columns;
    vnaa.vnaa_s_matrix		= s;
    vnaa.vnaa_s_rows		= s_rows;
    vnaa.vnaa_s_columns		= s_columns;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= false;
    vnaa.vnaa_m_type		= 'a';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}

/*
 * vnacal_new_add_mapped_matrix_m: add a matrix of measurements with port map
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s: matrix of parameter indices (by rows)
 *   @s_rows: number of rows in s
 *   @s_columns: number of columns in s
 *   @port_map: vector of VNA port numbers corresponding to these ports
 */
int vnacal_new_add_mapped_matrix_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	const int *s, int s_rows, int s_columns,
	const int *port_map)
{
    vnacal_new_add_arguments_t vnaa;

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    (void)memset((void *)&vnaa, 0, sizeof(vnaa));
    vnaa.vnaa_function		= __func__;
    vnaa.vnaa_cmp		= vnp;
    vnaa.vnaa_a_matrix		= NULL;
    vnaa.vnaa_a_rows		= 0;
    vnaa.vnaa_a_columns		= 0;
    vnaa.vnaa_b_matrix		= (const double complex *const *)m;
    vnaa.vnaa_b_rows		= m_rows;
    vnaa.vnaa_b_columns		= m_columns;
    vnaa.vnaa_s_matrix		= s;
    vnaa.vnaa_s_rows		= s_rows;
    vnaa.vnaa_s_columns		= s_columns;
    vnaa.vnaa_m_is_diagonal	= false;
    vnaa.vnaa_s_is_diagonal	= false;
    vnaa.vnaa_m_type		= 'm';
    vnaa.vnaa_s_port_map	= port_map;

    return _vnacal_new_add_common(vnaa);
}
