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
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _vnacal_type_to_name: convert type to type name
 *   @type: error term type
 */
const char *_vnacal_type_to_name(vnacal_type_t type)
{
    switch (type) {
    case VNACAL_T8:
	return "T8";
    case VNACAL_U8:
	return "U8";
    case VNACAL_TE10:
	return "TE10";
    case VNACAL_UE10:
	return "UE10";
    case VNACAL_T16:
	return "T16";
    case VNACAL_U16:
	return "U16";
    case VNACAL_UE14:
	return "UE14";
    case _VNACAL_E12_UE14:
	return "E12_UE14";
    case VNACAL_E12:
	return "E12";
    }
    abort();
}

/*
 * _vnacal_type_to_name: convert type name to type
 *   @name: error parameter type name
 */
vnacal_type_t _vnacal_name_to_type(const char *name)
{
    switch (name[0]) {
    case 'E':
	if (strcmp(name, "E12") == 0) {
	    return VNACAL_E12;
	}
	if (strcmp(name, "E12_UE14") == 0) {
	    return _VNACAL_E12_UE14;
	}
	break;

    case 'T':
	switch (name[1]) {
	case '1':
	    if (strcmp(name, "T16") == 0) {
		return VNACAL_T16;
	    }
	    break;

	case '8':
	    if (strcmp(name, "T8") == 0) {
		return VNACAL_T8;
	    }
	    break;

	case 'E':
	    if (strcmp(name, "TE10") == 0) {
		return VNACAL_TE10;
	    }
	    break;

	default:
	    break;
	}
	break;

    case 'U':
	switch (name[1]) {
	case '1':
	    if (strcmp(name, "U16") == 0) {
		return VNACAL_U16;
	    }
	    break;

	case '8':
	    if (strcmp(name, "U8") == 0) {
		return VNACAL_U8;
	    }
	    break;

	case 'E':
	    if (strcmp(name, "UE10") == 0) {
		return VNACAL_UE10;
	    }
	    if (strcmp(name, "UE14") == 0) {
		return VNACAL_UE14;
	    }
	    break;

	default:
	    break;
	}
	break;

    default:
	break;
    }
    return (vnacal_type_t)-1;
}

/*
 * _vnacal_layout: init the error term layout structure
 *  @vlp: pointer to vnacal_layout_t structure
 *  @type: error term type
 *  @m_rows: number of VNA ports that detect signal
 *  @m_columns: number of VNA ports that generate signal
 *
 * Notes:
 *  T16, TE10 and T8 are scattering-transfer, "T", parameters, while U16,
 *  UE10, UE14 and U8 are inverse scattering-transfer, "U", parameters.
 *  Both T and U matrices are, by definition, always 2x2, but in this
 *  context, the elements of the matrices are themselves matrices.
 *
 *  We refer to the four sub-matrices of T as Ts, Ti, Tx and Tm, and
 *  the four sub-matrices of U as Um, Ui, Ux and Us, indicating their
 *  coefficients in the following matrix equations:
 *
 *    Ts S + Ti = M Tx S + M Tm
 *    Um M + Ui = S Ux M + S Us
 *
 *    with dimensions:
 *      Ts: m_rows    x s_rows         Um: s_rows    x m_rows
 *      Ti: m_rows    x s_columns      Ui: s_rows    x m_columns
 *      Tx: m_columns x s_rows         Ux: s_columns x m_rows
 *      Tm: m_columns x s_columns      Us: s_columns x m_columns
 *      S: s_rows     x s_columns
 *      M: m_rows     x m_columns
 *
 *    where:
 *      Ts, Us: error terms with coefficients involving only S
 *      Ti, Ui: error terms with coefficients that are only 1 or 0
 *      Tx, Ux: error terms with coefficients involving both M and S
 *      Tm, Um: error terms with coefficients involving on M
 *      S:  s-parameter matrix of a calibration standard
 *      M:  measured values of the calibration standard
 *
 *  While the matrix equations above are valid for any values of m_rows,
 *  m_columns, s_rows and s_columns, we apply some practical constraints.
 *  In T, we constrain m_rows <= m_columns.  In U, we constrain m_rows
 *  >= m_columns.  In both systems, we set s_rows and s_columns to
 *  MAX(m_rows, m_columns).  These avoid systems with more equations
 *  than measurements and make it possible to solve both matrix equations
 *  for M, ensuring that the matrix to be inverted is square.
 *
 *    M = (Ts S + Ti) (Tx S + Tm)^-1
 *      = (Um - S Ux)^-1 (S Us - Ui)
 *
 *  In T16 and U16, the four sub-matrices are complete.  In T8, U8, TE10,
 *  UE10 and UE14, the sub-matrices are diagonal matrices and only the
 *  diagonal elements are stored.
 *
 *  In TE10, UE10 and UE14, we also include a scattering parameter matrix,
 *  El, containing the off-diagonal leakage terms.  This matrix contains
 *  only of the off diagonal elements.
 *
 *  UE14 is a generalization of the classic SOLT 12-term calibration where
 *  each column of the measurement matrix has its own independent error
 *  parameters, i.e. it comprises m_columns separate m_rows x 1 systems.
 *  UE14 has the advantage that it can compensate for a switch placed
 *  on the DUT side of the reflection bridges.
 */
void _vnacal_layout(vnacal_layout_t *vlp, vnacal_type_t type,
	int m_rows, int m_columns)
{
    const int diagonals = MIN(m_rows, m_columns);
    const int ports     = MAX(m_rows, m_columns);
    const int s_rows    = ports;
    const int s_columns = ports;

    /*
     * Set universal members.
     */
    (void)memset((void *)vlp, 0, sizeof(*vlp));
    vlp->vl_type      = type;
    vlp->vl_m_rows    = m_rows;
    vlp->vl_m_columns = m_columns;

    /*
     * Set per-type members.
     */
    switch (type) {
    case VNACAL_T16:
	{
	    const int ts_offset = 0;
	    const int ti_offset = ts_offset + m_rows    * s_rows;
	    const int tx_offset = ti_offset + m_rows    * s_columns;
	    const int tm_offset = tx_offset + m_columns * s_rows;
	    const int t_terms	= tm_offset + m_columns * s_columns;

	    vlp->vl_ti_offset	 = ti_offset;
	    vlp->vl_tx_offset	 = tx_offset;
	    vlp->vl_tm_offset	 = tm_offset;
	    vlp->vl_t_terms	 = t_terms;
	    vlp->vl_el_offset	 = t_terms;
	    vlp->vl_el_terms	 = 0;
	    vlp->vl_error_terms  = t_terms;
	}
	break;

    case VNACAL_TE10:
    case VNACAL_T8:
	{
	    const int ts_offset = 0;
	    const int ti_offset = ts_offset + MIN(m_rows, s_rows);
	    const int tx_offset = ti_offset + MIN(m_rows, s_columns);
	    const int tm_offset = tx_offset + MIN(m_columns, s_rows);
	    const int t_terms	= tm_offset + MIN(m_columns, s_columns);
	    const int el_terms  = type == VNACAL_TE10 ?
		m_rows * m_columns - diagonals : 0;

	    vlp->vl_ti_offset	 = ti_offset;
	    vlp->vl_tx_offset	 = tx_offset;
	    vlp->vl_tm_offset	 = tm_offset;
	    vlp->vl_t_terms	 = t_terms;
	    vlp->vl_el_offset	 = t_terms;
	    vlp->vl_el_terms	 = el_terms;
	    vlp->vl_error_terms  = t_terms + el_terms;
	}
	break;

    case VNACAL_U16:
	{
	    const int um_offset = 0;
	    const int ui_offset = um_offset + s_rows    * m_rows;
	    const int ux_offset = ui_offset + s_rows    * m_columns;
	    const int us_offset = ux_offset + s_columns * m_rows;
	    const int u_terms	= us_offset + s_columns * m_columns;

	    vlp->vl_ui_offset	 = ui_offset;
	    vlp->vl_ux_offset	 = ux_offset;
	    vlp->vl_us_offset	 = us_offset;
	    vlp->vl_u_terms	 = u_terms;
	    vlp->vl_el_offset	 = u_terms;
	    vlp->vl_el_terms	 = 0;
	    vlp->vl_error_terms  = u_terms;
	}
	break;

    case VNACAL_UE10:
    case VNACAL_U8:
	{
	    const int um_offset = 0;
	    const int ui_offset = um_offset + MIN(s_rows, m_rows);
	    const int ux_offset = ui_offset + MIN(s_rows, m_columns);
	    const int us_offset = ux_offset + MIN(s_columns, m_rows);
	    const int u_terms	= us_offset + MIN(s_columns, m_columns);
	    const int el_terms  = type == VNACAL_UE10 ?
		m_rows * m_columns - diagonals : 0;

	    vlp->vl_ui_offset	 = ui_offset;
	    vlp->vl_ux_offset	 = ux_offset;
	    vlp->vl_us_offset	 = us_offset;
	    vlp->vl_u_terms	 = u_terms;
	    vlp->vl_el_offset    = u_terms;
	    vlp->vl_el_terms     = el_terms;
	    vlp->vl_error_terms  = u_terms + el_terms;
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	{
	    const int um_offset = 0;
	    const int ui_offset = um_offset + MIN(s_rows, m_rows);
	    const int ux_offset = ui_offset + 1;
	    const int us_offset = ux_offset + MIN(s_columns, m_rows);
	    const int u_terms	= us_offset + 1;
	    const int el_terms  = m_rows * m_columns - diagonals;

	    vlp->vl_ui_offset	 = ui_offset;
	    vlp->vl_ux_offset	 = ux_offset;
	    vlp->vl_us_offset	 = us_offset;
	    vlp->vl_u_terms	 = u_terms;
	    vlp->vl_el_offset    = m_columns * u_terms;
	    vlp->vl_el_terms     = el_terms;
	    vlp->vl_error_terms  = m_columns * u_terms + el_terms;
	}
	break;

    case VNACAL_E12:
	{
	    const int el_terms  = m_rows;
	    const int er_terms  = m_rows;
	    const int et_terms  = 0;			/* not stored */
	    const int em_terms  = m_rows;
	    const int el_offset = 0;
	    const int er_offset = el_offset + el_terms;
	    const int et_offset = er_offset + er_terms;
	    const int em_offset = et_offset + et_terms;
	    const int e_terms   = et_offset + em_terms;

	    vlp->vl_el_offset	= el_offset;
	    vlp->vl_er_offset	= er_offset;
	    vlp->vl_et_offset	= et_offset;
	    vlp->vl_em_offset	= em_offset;
	    vlp->vl_e_terms	= e_terms;
	    vlp->vl_el_terms    = el_terms;
	    vlp->vl_error_terms	= m_columns * e_terms;
	}
	break;
    }
}
