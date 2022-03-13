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

#ifndef _VNACAL_LAYOUT_H
#define _VNACAL_LAYOUT_H

#include "vnacal.h"

/*
 * vnacal_layout_t: type and layout of the error terms
 */
typedef struct vnacal_layout {
    /* error term type */
    vnacal_type_t vl_type;

    /* number of rows in the measurement matrix */
    int vl_m_rows;

    /* number of columns in the measurement matrix */
    int vl_m_columns;

    /* offset of the Ti, Ui and Er matrices */
    int vl_ti_offset;

    /* offset of the Tx, Ux and Et matrices */
    int vl_tx_offset;

    /* offset of the Tm, Us and Em matrices */
    int vl_tm_offset;

    /* number of terms in the full T, U or E matrix */
    int vl_t_terms;

    /* offset of the El matrix */
    int vl_el_offset;

    /* number of terms in the El matrix */
    int vl_el_terms;

    /* total number of error terms */
    int vl_error_terms;

} vnacal_layout_t;

/*
 * Alias members for U
 */
#define vl_ui_offset	vl_ti_offset
#define vl_ux_offset	vl_tx_offset
#define vl_us_offset	vl_tm_offset
#define vl_u_terms	vl_t_terms

/*
 * Alias members for E
 */
#define vl_er_offset	vl_ti_offset
#define vl_et_offset	vl_tx_offset
#define vl_em_offset	vl_tm_offset
#define vl_e_terms	vl_t_terms

/*
 * VL_TYPE: return the error parameter type
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_TYPE(vlp)		((vlp)->vl_type)

/*
 * VL_M_ROWS, VL_M_COLUMNS: return the dimensions of the measurement matrix
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_M_ROWS(vlp)		((vlp)->vl_m_rows)
#define VL_M_COLUMNS(vlp)	((vlp)->vl_m_columns)

/*
 * VL_M_PORTS: return the number of VNA ports
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_M_PORTS(vlp)		MAX((vlp)->vl_m_rows, (vlp)->vl_m_columns)

/*
 * VL_S_ROWS, VL_S_COLUMNS: return the dimensions of the s-parameter matrix
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_S_ROWS(vlp)		VL_M_PORTS(vlp)
#define VL_S_COLUMNS(vlp)	VL_M_PORTS(vlp)

/*
 * VL_ERROR_TERMS: return the total number of error terms
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_ERROR_TERMS(vlp)	((vlp)->vl_error_terms)

/*
 * VNACAL_IS_T: test if we're using T parameters
 *   @type: error term type
 */
#define VNACAL_IS_T(type) \
    ((type) == VNACAL_T8 || (type) == VNACAL_TE10 || (type) == VNACAL_T16)

/*
 * VL_IS_T: test if we're using T parameters
 *   @type: error term type
 */
#define	VL_IS_T(vlp) \
	VNACAL_IS_T((vlp)->vl_type)

/*
 * VNACAL_IS_UE14: test if type is VNACAL_UE14 or _VNACAL_E12_UE14
 *   @type: error term type
 */
#define VNACAL_IS_UE14(type) \
	((type) == VNACAL_UE14 || (type) == _VNACAL_E12_UE14)

/*
 * VNACAL_IS_UE14: test if the type is VNACAL_UE14 or _VNACAL_E12_UE14
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_IS_UE14(vlp) \
	VNACAL_IS_UE14((vlp)->vl_type)

/*
 * VNACAL_HAS_COLUMN_SYSTEMS: return true if each column is a separate system
 *   @type: error term type
 */
#define VNACAL_HAS_COLUMN_SYSTEMS(type) \
	((type) == VNACAL_E12 || VNACAL_IS_UE14(type))

/*
 * VNACAL_HAS_COLUMN_SYSTEMS: return true if each column is a separate system
 *   @vlp: pointer to vnacal_layout_t structure
 */
#define VL_HAS_COLUMN_SYSTEMS(vlp) \
	VNACAL_HAS_COLUMN_SYSTEMS((vlp)->vl_type)


/***********************************************************************
 * T terms: 2x2 T matrix of matrices:
 *
 *	[ Ts Ti ]
 *	[ Tx Tm ]
 *
 * The sub-matrices are related to the measurement matrix, M, and the
 * scattering parameter matrix, S by the following matrix equation:
 *
 *	Ts S + Ti = M Tx S + M Tm
 *
 * The dimensions follow directly from the equation:
 *	Ts: m_rows x s_rows
 *	Ti: m_rows x s_columns
 *	Tx: m_columns x s_rows
 *	Tm: m_columns x s_columns
 *
 * The sub-matrices are named based on their coefficients:
 *	Ts terms have only S coefficients.
 *	Ti terms have only 1 or 0 (identity matrix) coefficients
 *	Tx terms have both M and S coefficients
 *	Tm terms have only M coefficients.
 *
 * In T16 error terms, the four sub-matrices are complete.  In T8 and
 * TE10, the four sub-matrices are diagonal and only the diagonal elements
 * are stored.
 *
 * Solved in terms of M and S:
 *	M = (Ts S + Ti) (Tx S + Tm)^-1		if m_columns == s_columns
 *	S = (Ts - M Tx)^-1 (M Tm - Ti)		if m_rows    == s_rows
 *
 * From E terms:				if m_columns == s_columns
 *	Ts = Er - El Et^-1 Em
 *	Ti = El Et^-1
 *	Tx = - Et^-1 Em
 *	Tm = Et^-1
 *
 * From U terms:				if M and S square
 *	Ts = (Um - Ui Us^-1 Ux)^-1
 *	Ti = -Um^-1 Ui (Us - Ux Um^-1 Ui)^-1
 *	Tx = -Us^-1 Ux (Um - Ui Us^-1 Ux)^-1
 *	Tm = (Us - Ux Um^-1 Ui)^-1
 **********************************************************************/

#define VL_TS_ROWS(vlp)		VL_M_ROWS(vlp)
#define VL_TS_COLUMNS(vlp)	VL_S_ROWS(vlp)
#define VL_TS_TERMS(vlp)	((vlp)->vl_ti_offset - 0)
#define VL_TS_OFFSET(vlp)	(0)

#define VL_TI_ROWS(vlp)		VL_M_ROWS(vlp)
#define VL_TI_COLUMNS(vlp)	VL_S_COLUMNS(vlp)
#define VL_TI_TERMS(vlp)	((vlp)->vl_tx_offset - (vlp)->vl_ti_offset)
#define VL_TI_OFFSET(vlp)	((vlp)->vl_ti_offset)

#define VL_TX_ROWS(vlp)		VL_M_COLUMNS(vlp)
#define VL_TX_COLUMNS(vlp)	VL_S_ROWS(vlp)
#define VL_TX_TERMS(vlp)	((vlp)->vl_tm_offset - (vlp)->vl_tx_offset)
#define VL_TX_OFFSET(vlp)	((vlp)->vl_tx_offset)

#define VL_TM_ROWS(vlp)		VL_M_COLUMNS(vlp)
#define VL_TM_COLUMNS(vlp)	VL_S_COLUMNS(vlp)
#define VL_TM_TERMS(vlp)	((vlp)->vl_t_terms - (vlp)->vl_tm_offset)
#define VL_TM_OFFSET(vlp)	((vlp)->vl_tm_offset)


/***********************************************************************
 * U terms: 2x2 U matrix (inverse of T) of matrices:
 *
 *	[ Um Ui ]
 *	[ Ux Us ]
 *
 * The sub-matrices are related to the measurement matrix, M, and the
 * scattering parameter matrix, S by the following matrix equation:
 *
 *	Um M + Ui = S (Ux M + Us)
 *
 * The dimensions follow directly from the equation:
 *	Um: s_rows x m_rows
 *	Ui: s_rows x m_columns
 *	Ux: s_columns x m_rows
 *	Us: s_columns x m_columns
 *
 * In U16 error terms, the four sub-matrices are complete.  In U8 and
 * UE10, the four sub-matrices are diagonal and only the diagonal elements
 * are stored.
 *
 * Solved in terms of M and S:
 *	M = (Um - S Ux)^-1 (S Us - Ui)		if m_rows    == s_rows
 *	S = (Um M + Ui) (Ux M + Us)^-1		if m_columns == s_columns
 *
 * From E terms:
 *	Um = Er^-1				if m_rows    == s_rows
 *	Ui = -Er^-1 El
 *	Ux = Em Er^-1
 *	Us = Et - Em Er^-1 El
 *
 * From T terms:				if M and S square
 *	Um = (Ts - Ti Tm^-1 Tx)^-1
 *	Ui = -Ts^-1 Ti (Tm - Tx Ts^-1 Ti)^-1
 *	Ux = -Tm^-1 Tx (Ts - Ti Tm^-1 Tx)^-1
 *	Us = (Tm - Tx Ts^-1 Ti)^-1
 **********************************************************************/

#define VL_UM_ROWS(vlp)		VL_S_ROWS(vlp)
#define VL_UM_COLUMNS(vlp)	VL_M_ROWS(vlp)
#define VL_UM_TERMS(vlp)	((vlp)->vl_ui_offset - 0)
#define VL_UM_OFFSET(vlp)	(0)

#define VL_UI_ROWS(vlp)		VL_S_ROWS(vlp)
#define VL_UI_COLUMNS(vlp)	VL_M_COLUMNS(vlp)
#define VL_UI_TERMS(vlp)	((vlp)->vl_ux_offset - (vlp)->vl_ui_offset)
#define VL_UI_OFFSET(vlp)	((vlp)->vl_ui_offset)

#define VL_UX_ROWS(vlp)		VL_S_COLUMNS(vlp)
#define VL_UX_COLUMNS(vlp)	VL_M_ROWS(vlp)
#define VL_UX_TERMS(vlp)	((vlp)->vl_us_offset - (vlp)->vl_ux_offset)
#define VL_UX_OFFSET(vlp)	((vlp)->vl_ux_offset)

#define VL_US_ROWS(vlp)		VL_S_COLUMNS(vlp)
#define VL_US_COLUMNS(vlp)	VL_M_COLUMNS(vlp)
#define VL_US_TERMS(vlp)	((vlp)->vl_u_terms - (vlp)->vl_us_offset)
#define VL_US_OFFSET(vlp)	((vlp)->vl_us_offset)


/***********************************************************************
 * UE14 is a sequence of independent mx1 systems, one for each m_column:
 *
 *      Um is an s_rows x m_rows    diagonal matrix
 *      Ui is an s_rows x 1         diagonal matrix
 *      Ux is an s_columns x m_rows diagonal matrix
 *      Us is an s_columns x 1      diagonal matrix
 *      El is m_rows x 1            rectangular matrix, off diagonal terms only
 *
 * We can solve for M column by column as follows:
 *
 *  M(:,1) = El(:,1) + A^-1 B with:
 *      A = c1_um11 - c1_ux11 s11               - c1_ux22 s12
 *                  - c1_ux11 s21       c1_um22 - c1_ux22 s22
 *
 *      B = c1_us11 + c1_us11 s11
 *                    c1_us11 s21
 *
 *  M(:,2) = El(:,2) + A^-1 B with:
 *      A = c2_um11 - c2_ux11 s11               - c2_ux22 s12
 *                  - c2_ux11 s21       c2_um22 - c2_ux22 s22
 *
 *      B = c2_us11 + c2_us11 s11
 *                    c2_us11 s21
 *
 *  ...
 *
 * We can solve for S column by column as follows:
 *
 *  S = B A^-1 with:
 *      A = c1_ux11 m11 + c1_us11        c2_ux11 m12			...
 *          c1_ux22 m21                  c2_ux22 m22 + c2_us11
 *
 *      B = c1_um11 m11 + c1_ui11        c2_um11 m12			...
 *          c1_um22 m21                  c2_um22 m22 + c2_ui11
 *
 *  where we've already subtracted El from M.
 *
 **********************************************************************/

#define VL_UM14_ROWS(vlp)	VL_S_ROWS(vlp)
#define VL_UM14_COLUMNS(vlp)	VL_M_ROWS(vlp)
#define VL_UM14_TERMS(vlp)	((vlp)->vl_ui_offset - 0)
#define VL_UM14_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_u_terms)

#define VL_UI14_ROWS(vlp)	VL_S_ROWS(vlp)
#define VL_UI14_COLUMNS(vlp)	(1)
#define VL_UI14_TERMS(vlp)	((vlp)->vl_ux_offset - (vlp)->vl_ui_offset)
#define VL_UI14_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_u_terms + (vlp)->vl_ui_offset)

#define VL_UX14_ROWS(vlp)	VL_S_COLUMNS(vlp)
#define VL_UX14_COLUMNS(vlp)	VL_M_ROWS(vlp)
#define VL_UX14_TERMS(vlp)	((vlp)->vl_us_offset - (vlp)->vl_ux_offset)
#define VL_UX14_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_u_terms + (vlp)->vl_ux_offset)

#define VL_US14_ROWS(vlp)	VL_S_COLUMNS(vlp)
#define VL_US14_COLUMNS(vlp)	(1)
#define VL_US14_TERMS(vlp)	((vlp)->vl_u_terms - (vlp)->vl_us_offset)
#define VL_US14_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_u_terms + (vlp)->vl_us_offset)


/***********************************************************************
 * E terms: 2x2 S matrix of matrices:
 *
 *	[ El Er ]
 *	[ Et Em ]
 *
 * The sub-matrices are related to the measurement matrix, M, and the
 * scattering parameter matrix, S by the following matrix equations:
 *
 *      (I - S Em) Er^-1 (M - El) = S Et	if m_rows    == s_rows
 *      (M - El) Et^-1 (I - Em S) = Er S	if m_columns == s_columns
 *
 * The dimensions follow from the equations:
 *	El: m_rows x m_columns		leakage
 *	Er: m_rows x s_rows		reflection tracking
 *	Et: s_columns x m_columns	transmission tracking
 *	Em: s_columns x s_rows		match
 *
 * The El sub-matrix matches the dimensions of M.  The other matrices
 * may be complete or diagonal, depending on the number of error terms.
 *
 * Solved for M:
 *      M = El + Er (I - S Em)^-1 S Et
 *        = El + Er S (I - Em S)^-1 Et
 *
 * Solved for S:
 *      S = Er^-1 (M - El) (Et + Em Er^-1 (M - El))^-1	if mr == sr && mc == sc
 *        = ((M - El) Et^-1 Em + Er)^-1 (M - El) Et^-1	if mr == sr && mc == sc
 *
 * From T terms:				if m_columns == s_columns
 *	El = Ti Tm^-1
 *	Er = Ts - Ti Tm^-1 Tx
 *	Et = Tm^-1
 *	Em = -Tm^-1 Tx
 *
 * From U terms:				if s_rows == m_rows
 *	El = -Um^-1 Ui
 *	Er =  Um^-1
 *	Et =  Us - Ux Um^-1 Ui
 *	Em =  Ux Um^-1
 *
 **********************************************************************/

#define VL_EL_ROWS(vlp)		VL_M_ROWS(vlp)
#define VL_EL_COLUMNS(vlp)	VL_M_COLUMNS(vlp)
#define VL_EL_TERMS(vlp)	((vlp)->vl_el_terms)
#define VL_EL_OFFSET(vlp)	((vlp)->vl_el_offset)

#define VL_ER_ROWS(vlp)		VL_M_ROWS(vlp)
#define VL_ER_COLUMNS(vlp)	VL_S_ROWS(vlp)
#define VL_ER_TERMS(vlp)	((vlp)->vl_et_offset - (vlp)->vl_er_offset)
#define VL_ER_OFFSET(vlp)	((vlp)->vl_er_offset)

#define VL_ET_ROWS(vlp)		VL_S_COLUMNS(vlp)
#define VL_ET_COLUMNS(vlp)	VL_M_COLUMNS(vlp)
#define VL_ET_TERMS(vlp)	((vlp)->vl_em_offset - (vlp)->vl_et_offset)
#define VL_ET_OFFSET(vlp)	((vlp)->vl_et_offset)

#define VL_EM_ROWS(vlp)		VL_S_COLUMNS(vlp)
#define VL_EM_COLUMNS(vlp)	VL_S_ROWS(vlp)
#define VL_EM_TERMS(vlp)	((vlp)->vl_e_terms - (vlp)->vl_em_offset)
#define VL_EM_OFFSET(vlp)	((vlp)->vl_em_offset)


/***********************************************************************
 * E12 is a sequence of independent systems, one for each m_column where:
 *
 *	El is an m_rows x 1         column vector
 *	Er is an m_rows x s_rows    diagonal matrix
 *	Et is an s_columns x 1      rectangular matrix, which contains the
 *	   m_column'th column of the identity matrix, and isn't stored
 *	Em is an s_columns x s_rows diagonal matrix
 *
 * We can solve for M column by column as follows:
 *
 *  M(:,1) = El(:,1) + B A^-1 Et(:,1) with:
 *	A = [ 1 - c1_em11 s11			    - c1_em11 s12 ]
 *	    [	- c1_em22 s21			  1 - c1_em22 s22 ]
 *
 *	B = [ c1_er11 s11			      c1_er11 s12 ]
 *	    [ c1_er22 s21			      c1_er22 s22 ]
 *
 *  M(:,2) = El(:,2) + B A^-1 Et(:,2) with:
 *	A = [ 1 - c2_em11 s11			    - c2_em11 s12 ]
 *	    [	- c2_em22 s21			  1 - c2_em22 s22 ]
 *
 *	B = [ c2_er11 s11			      c2_er11 s12 ]
 *	    [ c2_er22 s21			      c2_er22 s22 ]
 *
 *  ...
 *
 * We can solve for S column by column as follows:
 *
 *  S = B A^-1 with:
 *      B = (m11 - c1_el11) / c1_er11		(m12 - c2_el12) / c2_er11  ...
 *          (m21 - c1_el11) / c1_er22		(m22 - c2_el22) / c2_er22
 *
 *      A = 1 + c1_em11 b11			    c2_em11 b12            ...
 *              c1_em22 b21                     1 + c2_em22 b22
 *
 **********************************************************************/

#define VL_EL12_ROWS(vlp)	VL_M_ROWS(vlp)
#define VL_EL12_COLUMNS(vlp)	(1)
#define VL_EL12_TERMS(vlp)	((vlp)->vl_el_terms)
#define VL_EL12_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_e_terms + 0)

#define VL_ER12_ROWS(vlp)	VL_M_ROWS(vlp)
#define VL_ER12_COLUMNS(vlp)	VL_S_ROWS(vlp)
#define VL_ER12_TERMS(vlp)	((vlp)->vl_et_offset - (vlp)->vl_er_offset)
#define VL_ER12_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_e_terms + (vlp)->vl_er_offset)

#define VL_ET12_ROWS(vlp)	VL_S_COLUMNS(vlp)
#define VL_ET12_COLUMNS(vlp)	(1)
#define VL_ET12_TERMS(vlp)	(0)
#define VL_ET12_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_e_terms + (vlp)->vl_et_offset)

#define VL_EM12_ROWS(vlp)	VL_S_COLUMNS(vlp)
#define VL_EM12_COLUMNS(vlp)	VL_S_ROWS(vlp)
#define VL_EM12_TERMS(vlp)	((vlp)->vl_e_terms - (vlp)->vl_em_offset)
#define VL_EM12_OFFSET(vlp, m_column) \
	((m_column) * (vlp)->vl_e_terms + (vlp)->vl_em_offset)


/*
 * _vl_unity_offset: return the index of the term that always has value 1.0
 *   @vlp: pointer to vnacal_layout_t structure
 *   @column: the column in UE14
 */
static inline int _vl_unity_offset(const vnacal_layout_t *vlp, int system)
{
    switch (vlp->vl_type) {
    case VNACAL_T8:			/* always tm11 */
    case VNACAL_TE10:
    case VNACAL_T16:
	return VL_TM_OFFSET(vlp);

    case VNACAL_U8:			/* always um11 */
    case VNACAL_UE10:
    case VNACAL_U16:
	return VL_UM_OFFSET(vlp);

    case VNACAL_UE14:			/* unity term is column */
    case _VNACAL_E12_UE14:
	return system;

    case VNACAL_E12:			/* no unity term in E12 */
	break;

    default:
	break;
    }
    return -1;
}

#endif /* _VNACAL_LAYOUT_H */
