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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_apply_args_t: arguments to _vnacal_apply_common
 */
typedef struct vnacal_apply_args {
    /* name of user-called function */
    const char			       *vaa_function;

    /* calibration to use */
    vnacal_t			       *vaa_vcp;
    int					vaa_ci;

    /* vector of frequencies */
    const double		       *vaa_frequency_vector;
    int					vaa_frequencies;

    /* matrix of voltages leaving each VNA port */
    const double complex        *const *vaa_a_matrix;
    int					vaa_a_rows;
    int					vaa_a_columns;

    /* matrix of voltages entering each VNA port */
    const double complex        *const *vaa_b_matrix;
    int					vaa_b_rows;
    int					vaa_b_columns;

    /* 'a' for a/b or 'm' for m */
    char				vaa_m_type;

    /* result */
    vnadata_t			       *vaa_s_parameters;

} vnacal_apply_args_t;

/*
 * fill_t8: fill in the A & B matrices for VNACAL_T8, VNACAL_TE10
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error terms
 *   @m: measurement matrix
 *   @a: s-parameter coefficient matrix
 *   @b: right-hand side matrix
 */
static void fill_t8(const vnacal_layout_t *vlp,
	const double complex *e, double complex *m,
	double complex *a, double complex *b)
{
    const vnacal_type_t type = VL_TYPE(vlp);
    const int m_rows         = VL_M_ROWS(vlp);
    const int m_columns      = VL_M_COLUMNS(vlp);
    const int s_rows         = VL_S_ROWS(vlp);
    const int s_columns      = VL_S_COLUMNS(vlp);
    const double complex *const ts = &e[VL_TS_OFFSET(vlp)];
    const double complex *const ti = &e[VL_TI_OFFSET(vlp)];
    const double complex *const tx = &e[VL_TX_OFFSET(vlp)];
    const double complex *const tm = &e[VL_TM_OFFSET(vlp)];
    const double complex *const el = &e[VL_EL_OFFSET(vlp)];

    /*
     * Special-case 2x2 M with a 1x2 calibration.
     */
    if (m_rows == 1 && m_columns == 2) {
	if (type == VNACAL_TE10) {
	    m[1] -= el[0];		/* m12 -= el12 */
	    m[2] -= el[0];		/* m21 -= el12 */
	}
	a[0] =  ts[0] - m[0] * tx[0];	/* a11 =  ts11 - m11 tx11 */
	a[1] =        - m[1] * tx[1];	/* a12 =       - m12 tx22 */
	a[2] =        - m[2] * tx[1];   /* a21 =       - m21 tx22 */
	a[3] =  ts[0] - m[3] * tx[0];   /* a22 =  ts11 - m22 tx11 */
	b[0] = -ti[0] + m[0] * tm[0];	/* b11 = -ti11 + m11 tm11 */
	b[1] =          m[1] * tm[1];	/* b12 =         m12 tm22 */
	b[2] =          m[2] * tm[1];	/* b21 =         m21 tm22 */
	b[3] = -ti[0] + m[3] * tm[0];	/* b22 = -ti11 + m22 tm11 */
	return;
    }

    /*
     * If the calibration type has error terms handled outside of the
     * linear system, subtract those out of the M matrix.
     */
    assert(m_rows == m_columns);
    if (type == VNACAL_TE10) {
	const double complex *el_cur = el;

	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		if (m_row != m_column) {
		    int m_cell = m_row * m_columns + m_column;

		    m[m_cell] -= *el_cur++;
		}
	    }
	}
	assert(el_cur == &e[VL_EL_OFFSET(vlp) + VL_EL_TERMS(vlp)]);
    }

    /*
     * Build the A matrix.
     */
    for (int a_row = 0; a_row < m_rows; ++a_row) {
	for (int a_column = 0; a_column < s_rows; ++a_column) {
	    int a_cell = a_row * s_rows + a_column;

	    a[a_cell] = 0.0;
	    if (a_row == a_column) {
		a[a_cell] += ts[a_row];
	    }
	    a[a_cell] -= m[a_cell] * tx[a_column];
	}
    }

    /*
     * Build the B matrix.
     */
    for (int b_row = 0; b_row < m_rows; ++b_row) {
	for (int b_column = 0; b_column < s_columns; ++b_column) {
	    int b_cell = b_row * s_columns + b_column;

	    b[b_cell] = 0.0;
	    if (b_row == b_column) {
		b[b_cell] -= ti[b_row];
	    }
	    b[b_cell] += m[b_cell] * tm[b_column];
	}
    }
}

/*
 * fill_u8: fill in the A & B matrices for VNACAL_U8, VNACAL_UE10
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error terms
 *   @m: measurement matrix
 *   @a: s-parameter coefficient matrix
 *   @b: right-hand side matrix
 */
static void fill_u8(const vnacal_layout_t *vlp,
	const double complex *e, double complex *m,
	double complex *a, double complex *b)
{
    const vnacal_type_t type = VL_TYPE(vlp);
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    const double complex *const um = &e[VL_UM_OFFSET(vlp)];
    const double complex *const ui = &e[VL_UI_OFFSET(vlp)];
    const double complex *const ux = &e[VL_UX_OFFSET(vlp)];
    const double complex *const us = &e[VL_US_OFFSET(vlp)];
    const double complex *const el = &e[VL_EL_OFFSET(vlp)];

    /*
     * Special-case 2x2 M with a 2x1 calibration.
     */
    if (m_rows == 2 && m_columns == 1) {
	if (type == VNACAL_UE10) {
	    m[1] -= el[0];		/* m12 -= el12 */
	    m[2] -= el[0];		/* m21 -= el12 */
	}
	a[0] = us[0] + m[0] * ux[0];	/* a11 = us11 + m11 ux11 */
	a[1] =         m[1] * ux[1];	/* a12 =        m12 ux22 */
	a[2] =         m[2] * ux[1];    /* a21 =        m21 ux22 */
	a[3] = us[0] + m[3] * ux[0];    /* a22 = us11 + m22 ux11 */
	b[0] = ui[0] + m[0] * um[0];	/* b11 = ui11 + m11 um11 */
	b[1] =         m[1] * um[1];	/* b12 =        m12 um22 */
	b[2] =         m[2] * um[1];	/* b21 =        m21 um22 */
	b[3] = ui[0] + m[3] * um[0];	/* b22 = ui11 + m22 um11 */
	return;
    }

    /*
     * If the calibration type has error terms handled outside of the
     * linear system, subtract those out of the M matrix.
     */
    assert(m_rows == m_columns);
    if (type == VNACAL_UE10) {
	const double complex *el_cur = el;

	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		if (m_row != m_column) {
		    int m_cell = m_row * m_columns + m_column;

		    m[m_cell] -= *el_cur++;
		}
	    }
	}
	assert(el_cur == &e[VL_EL_OFFSET(vlp) + VL_EL_TERMS(vlp)]);
    }

    /*
     * Build the A matrix.
     */
    for (int a_row = 0; a_row < s_columns; ++a_row) {
	for (int a_column = 0; a_column < m_columns; ++a_column) {
	    int a_cell = a_row * m_columns + a_column;

	    a[a_cell] = m[a_cell] * ux[a_row];
	    if (a_row == a_column) {
		a[a_cell] += us[a_row];
	    }
	}
    }

    /*
     * Build the B matrix.
     */
    for (int b_row = 0; b_row < s_rows; ++b_row) {
	for (int b_column = 0; b_column < m_columns; ++b_column) {
	    int b_cell = b_row * m_columns + b_column;

	    b[b_cell] = m[b_cell] * um[b_row];
	    if (b_row == b_column) {
		b[b_cell] += ui[b_row];
	    }
	}
    }
}

/*
 * fill_t16: fill in the A & B matrices for VNACAL_T16
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error terms
 *   @m: measurement matrix
 *   @a: s-parameter coefficient matrix
 *   @b: right-hand side matrix
 */
static void fill_t16(const vnacal_layout_t *vlp,
	const double complex *e, double complex *m,
	double complex *a, double complex *b)
{
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    const double complex *const ts = &e[VL_TS_OFFSET(vlp)];
    const double complex *const ti = &e[VL_TI_OFFSET(vlp)];
    const double complex *const tx = &e[VL_TX_OFFSET(vlp)];
    const double complex *const tm = &e[VL_TM_OFFSET(vlp)];

    /*
     * Special-case 2x2 M with a 1x2 calibration.
     */
    if (m_rows == 1 && m_columns == 2) {
	a[0] =  ts[0] - m[0] * tx[0] - m[1] * tx[2];
	a[1] =  ts[1] - m[0] * tx[1] - m[1] * tx[3];
	a[2] =  ts[1] - m[2] * tx[3] - m[3] * tx[1];
	a[3] =  ts[0] - m[2] * tx[2] - m[3] * tx[0];
	b[0] = -ti[0] + m[0] * tm[0] + m[1] * tm[2];
	b[1] = -ti[1] + m[0] * tm[1] + m[1] * tm[3];
	b[2] = -ti[1] + m[2] * tm[3] + m[3] * tm[1];
	b[3] = -ti[0] + m[2] * tm[2] + m[3] * tm[0];
	return;
    }

    /*
     * Build the A matrix.
     */
    assert(m_rows == m_columns);
    for (int a_row = 0; a_row < m_rows; ++a_row) {
	for (int a_column = 0; a_column < s_rows; ++a_column) {
	    int a_cell = a_row * s_rows + a_column;

	    a[a_cell] = ts[a_cell];
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		int m_cell = a_row * m_columns + m_column;
		int tx_cell = m_column * s_rows + a_column;

		a[a_cell] -= m[m_cell] * tx[tx_cell];
	    }
	}
    }

    /*
     * Build the B matrix.
     */
    for (int b_row = 0; b_row < m_rows; ++b_row) {
	for (int b_column = 0; b_column < s_columns; ++b_column) {
	    int b_cell = b_row * s_columns + b_column;

	    b[b_cell] = -ti[b_cell];
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		int m_cell = b_row * m_columns + m_column;
		int tm_cell = m_column * s_columns + b_column;

		b[b_cell] += m[m_cell] * tm[tm_cell];
	    }
	}
    }
}

/*
 * fill_u16: fill in the A & B matrices for VNACAL_U16
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error terms
 *   @m: measurement matrix
 *   @a: s-parameter coefficient matrix
 *   @b: right-hand side matrix
 */
static void fill_u16(const vnacal_layout_t *vlp,
	const double complex *e, double complex *m,
	double complex *a, double complex *b)
{
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    const double complex *const um = &e[VL_UM_OFFSET(vlp)];
    const double complex *const ui = &e[VL_UI_OFFSET(vlp)];
    const double complex *const ux = &e[VL_UX_OFFSET(vlp)];
    const double complex *const us = &e[VL_US_OFFSET(vlp)];

    /*
     * Special-case 2x2 M with a 2x1 calibration.
     */
    if (m_rows == 2 && m_columns == 1) {
	a[0] = us[0] + m[0] * ux[0] + m[2] * ux[1];
	a[1] = us[1] + m[1] * ux[3] + m[3] * ux[2];
	a[2] = us[1] + m[0] * ux[2] + m[2] * ux[3];
	a[3] = us[0] + m[1] * ux[1] + m[3] * ux[0];
	b[0] = ui[0] + m[0] * um[0] + m[2] * um[1];
	b[1] = ui[1] + m[1] * um[3] + m[3] * um[2];
	b[2] = ui[1] + m[0] * um[2] + m[2] * um[3];
	b[3] = ui[0] + m[1] * um[1] + m[3] * um[0];
	return;
    }

    /*
     * Build the A matrix.
     */
    assert(m_rows == m_columns);
    for (int a_row = 0; a_row < s_columns; ++a_row) {
	for (int a_column = 0; a_column < m_columns; ++a_column) {
	    int a_cell = a_row * m_columns + a_column;

	    a[a_cell] = 0.0;
	    for (int m_row = 0; m_row < m_rows; ++m_row) {
		int m_cell = m_row * m_columns + a_column;
		int ux_cell = a_row * m_rows + m_row;

		a[a_cell] += m[m_cell] * ux[ux_cell];
	    }
	    a[a_cell] += us[a_cell];
	}
    }

    /*
     * Build the B matrix.
     */
    for (int b_row = 0; b_row < s_rows; ++b_row) {
	for (int b_column = 0; b_column < m_columns; ++b_column) {
	    int b_cell = b_row * m_columns + b_column;

	    b[b_cell] = 0.0;
	    for (int m_row = 0; m_row < m_rows; ++m_row) {
		int m_cell = m_row * m_columns + b_column;
		int um_cell = b_row * m_rows + m_row;

		b[b_cell] += m[m_cell] * um[um_cell];
	    }
	    b[b_cell] += ui[b_cell];
	}
    }
}

/*
 * fill_ue14: fill in the A & B matrices for VNACAL_UE14, _VNACAL_E12_UE14
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error terms
 *   @m: measurement matrix
 *   @a: s-parameter coefficient matrix
 *   @b: right-hand side matrix
 */
static void fill_ue14(const vnacal_layout_t *vlp,
	const double complex *e, double complex *m,
	double complex *a, double complex *b)
{
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);

    /*
     * Special-case 2x2 M with a 2x1 calibration.
     */
    if (m_rows == 2 && m_columns == 1) {
	const double complex *const um = &e[VL_UM_OFFSET(vlp)];
	const double complex *const ui = &e[VL_UI_OFFSET(vlp)];
	const double complex *const ux = &e[VL_UX_OFFSET(vlp)];
	const double complex *const us = &e[VL_US_OFFSET(vlp)];
	const double complex *const el = &e[VL_EL_OFFSET(vlp)];

	m[1] -= el[0];			/* m12 -= el12 */
	m[2] -= el[0];			/* m21 -= el12 */
	a[0] = us[0] + m[0] * ux[0];	/* a11 = 1_us11 + m11 1_ux11 */
	a[1] =         m[1] * ux[1];	/* a12 =          m12 1_ux22 */
	a[2] =         m[2] * ux[1];    /* a21 =          m21 1_ux22 */
	a[3] = us[0] + m[3] * ux[0];    /* a22 = 1_us11 + m22 1_ux11 */
	b[0] = ui[0] + m[0] * um[0];	/* b11 = 1_ui11 + m11 1_um11 */
	b[1] =         m[1] * um[1];	/* b12 =          m12 1_um22 */
	b[2] =         m[2] * um[1];	/* b21 =          m21 1_um22 */
	b[3] = ui[0] + m[3] * um[0];	/* b22 = 1_ui11 + m22 1_um11 */
	return;
    }

    /*
     * Subtract leakage terms handled outside of the linear system from M.
     */
    assert(m_rows == m_columns);
    {
	const double complex *el_cur = &e[VL_EL_OFFSET(vlp)];

	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		if (m_row != m_column) {
		    int m_cell = m_row * m_columns + m_column;

		    m[m_cell] -= *el_cur++;
		}
	    }
	}
	assert(el_cur == &e[VL_EL_OFFSET(vlp) + VL_EL_TERMS(vlp)]);
    }

    /*
     * In UE14, each column in the calibration represents an independent
     * c_rows x 1 system with its own separate error parameters.  Each
     * system contributes a single column to the A and B matrices.
     */
    for (int m_column = 0; m_column < m_columns; ++m_column) {
	const double complex *const um = &e[VL_UM14_OFFSET(vlp, m_column)];
	const double complex *const ui = &e[VL_UI14_OFFSET(vlp, m_column)];
	const double complex *const ux = &e[VL_UX14_OFFSET(vlp, m_column)];
	const double complex *const us = &e[VL_US14_OFFSET(vlp, m_column)];

	/*
	 * Add a column to the A matrix.
	 */
	for (int a_row = 0; a_row < s_columns; ++a_row) {
	    int a_cell = a_row * m_columns + m_column;

	    a[a_cell] = m[a_cell] * ux[a_row];
	    if (a_row == m_column) {
		a[a_cell] += us[0];
	    }
	}

	/*
	 * Add a column to the B matrix.
	 */
	for (int b_row = 0; b_row < s_rows; ++b_row) {
	    int b_cell = b_row * m_columns + m_column;

	    b[b_cell] = m[b_cell] * um[b_row];
	    if (b_row == m_column) {
		b[b_cell] += ui[0];
	    }
	}
    }
}

/*
 * fill_e12: fill in the A & B matrices for VNACAL_E12
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error terms
 *   @m: measurement matrix
 *   @a: s-parameter coefficient matrix
 *   @b: right-hand side matrix
 */
static void fill_e12(const vnacal_layout_t *vlp,
	const double complex *e, double complex *m,
	double complex *a, double complex *b)
{
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);

    /*
     * Special-case 2x2 M with a 2x1 calibration.
     */
    if (m_rows == 2 && m_columns == 1) {
	const double complex *const el = &e[VL_EL12_OFFSET(vlp, 0)];
	const double complex *const er = &e[VL_ER12_OFFSET(vlp, 0)];
	const double complex *const em = &e[VL_EM12_OFFSET(vlp, 0)];

	m[0] -= el[0];
	m[1] -= el[1];
	m[2] -= el[1];
	m[3] -= el[0];
	b[0] = m[0] / er[0];
	b[1] = m[1] / er[1];
	b[2] = m[2] / er[1];
	b[3] = m[3] / er[0];
	a[0] = 1.0 + em[0] * b[0];
	a[1] =       em[1] * b[1];
	a[2] =       em[1] * b[2];
	a[3] = 1.0 + em[0] * b[3];
	return;
    }

    /*
     * In E12, each column in the calibration represents an independent
     * m_rows x 1 system with its own separate error parameters.  Each
     * system contributes a single column to the A and B matrices.
     */
    assert(m_rows == m_columns);
    for (int m_column = 0; m_column < m_columns; ++m_column) {
	const double complex *const el = &e[VL_EL12_OFFSET(vlp, m_column)];
	const double complex *const er = &e[VL_ER12_OFFSET(vlp, m_column)];
	const double complex *const em = &e[VL_EM12_OFFSET(vlp, m_column)];

	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    int cell = m_row * m_columns + m_column;
	    double complex x;

	    x = (m[cell] - el[m_row]) / er[m_row];
	    a[cell] = (m_row == m_column ? 1.0 : 0.0) + em[m_row] * x;
	    b[cell] = x;
	}
    }
}

/*
 * _vnacal_apply_common: apply a calibration to a set of measurements
 *   @vaa: common argument structure
 */
int _vnacal_apply_common(vnacal_apply_args_t vaa)
{
    vnacal_t *vcp = vaa.vaa_vcp;
    const vnacal_calibration_t *calp;
    vnacal_type_t c_type;
    vnacal_layout_t vl;
    int c_rows, c_columns, c_ports;
    double fmin, fmax;
    int segment = 0;

    /*
     * Get the calibration and validate parameters.
     */
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if ((calp = _vnacal_get_calibration(vcp, vaa.vaa_ci)) == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid calibration index: %d",
		vaa.vaa_function, vaa.vaa_ci);
        return -1;
    }
    c_type    = calp->cal_type;
    c_rows    = calp->cal_rows;
    c_columns = calp->cal_columns;
    c_ports   = MAX(c_rows, c_columns);
    _vnacal_layout(&vl, c_type, c_rows, c_columns);
    if (c_rows != c_columns && c_ports != 2) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: this function can be used with "
		"1x2, 2x1 or square calibrations only",
		vaa.vaa_function);
	return -1;
    }
    if (vaa.vaa_frequency_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid NULL frequency vector",
		vaa.vaa_function);
        return -1;
    }
    if (vaa.vaa_frequencies < 0) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid frequency count: %d",
		vaa.vaa_function, vaa.vaa_frequencies);
	return -1;
    }
    for (int i = 0; i < vaa.vaa_frequencies - 1; ++i) {
	if (vaa.vaa_frequency_vector[i] >= vaa.vaa_frequency_vector[i + 1]) {
	    _vnacal_error(vcp, VNAERR_USAGE, "%s: non-increasing frequencies",
		    vaa.vaa_function);
	    return -1;
	}
    }
    fmin = _vnacal_calibration_get_fmin_bound(calp);
    if (vaa.vaa_frequency_vector[0] < fmin) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: frequency out of bounds %.3e < %.3e", vaa.vaa_function,
		vaa.vaa_frequency_vector[0], calp->cal_frequency_vector[0]);
	return -1;
    }
    fmax = _vnacal_calibration_get_fmax_bound(calp);
    if (vaa.vaa_frequency_vector[vaa.vaa_frequencies - 1] > fmax) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: frequency out of bounds %.3e > %.3e", vaa.vaa_function,
		vaa.vaa_frequency_vector[vaa.vaa_frequencies - 1],
		calp->cal_frequency_vector[calp->cal_frequencies - 1]);
	return -1;
    }
    if (vaa.vaa_b_matrix == NULL) {
	if (vaa.vaa_m_type == 'm') {
	    _vnacal_error(vcp, VNAERR_USAGE, "%s: invalid NULL m_matrix",
		    vaa.vaa_function);
	} else {
	    _vnacal_error(vcp, VNAERR_USAGE, "%s: invalid NULL b_matrix",
		    vaa.vaa_function);
	}
	return -1;
    }
    if (vaa.vaa_b_rows != c_ports || vaa.vaa_b_columns != c_ports) {
	if (vaa.vaa_m_type == 'm') {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: m_matrix must be %d x %d with this calibration",
		    vaa.vaa_function, c_ports, c_ports);
	} else {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: b_matrix must be %d x %d with this calibration",
		    vaa.vaa_function, c_ports, c_ports);
	}
	return -1;
    }
    for (int cell = 0; cell < c_ports * c_ports; ++cell) {
	if (vaa.vaa_b_matrix[cell] == NULL) {
	    if (vaa.vaa_m_type == 'm') {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: invalid NULL entry in m_matrix", vaa.vaa_function);
	    } else {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: invalid NULL entry in b_matrix", vaa.vaa_function);
	    }
	    return -1;
	}
    }
    if (vaa.vaa_a_matrix != NULL) {
	int a_rows = c_type == VNACAL_E12 || VNACAL_IS_UE14(c_type) ?
	    1 : c_ports;

	if (vaa.vaa_a_rows != a_rows || vaa.vaa_a_columns != c_ports) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: a_matrix must be %d x %d with this calibration",
		    vaa.vaa_function, a_rows, c_ports);
	    return -1;
	}
	for (int cell = 0; cell < a_rows * c_ports; ++cell) {
	    if (vaa.vaa_a_matrix[cell] == NULL) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"%s: invalid NULL entry in a_matrix",
			vaa.vaa_function);
		return -1;
	    }
	}
    }
    if (vaa.vaa_s_parameters == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid NULL s_parameters",
		vaa.vaa_function);
	return -1;
    }

    /*
     * Set up the output matrix.
     */
    if (vnadata_init(vaa.vaa_s_parameters, vaa.vaa_frequencies,
		c_ports, c_ports, VPT_S) == -1) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "vnadata_init: %s",
		strerror(errno));
	return -1;
    }

    /*
     * For each frequency index...
     */
    for (int findex = 0; findex < vaa.vaa_frequencies; ++findex) {
	double f = vaa.vaa_frequency_vector[findex];
	double complex determinant;
	double complex t[calp->cal_error_terms];
	double complex m[c_ports * c_ports];
	double complex a[c_ports * c_ports];
	double complex b[c_ports * c_ports];
	double complex s[c_ports * c_ports];
#define M(i, j)	(m[(i) * c_ports + (j)])
#define S(i, j)	(s[(i) * c_ports + (j)])

	/*
	 * Interpolate to find the error terms for this frequency.
	 */
	for (int term = 0; term < calp->cal_error_terms; ++term) {
	    t[term] = _vnacal_rfi(calp->cal_frequency_vector,
		    calp->cal_error_term_vector[term], calp->cal_frequencies,
		    MIN(calp->cal_frequencies, VNACAL_MAX_M), &segment, f);
	}

	/*
	 * Get the measured value matrix, M.  If no A matrix was
	 * given, simply copy the values from B.
	 */
	if (vaa.vaa_a_matrix == NULL) {
	    for (int cell = 0; cell < c_ports * c_ports; ++cell) {
		m[cell] = vaa.vaa_b_matrix[cell][findex];
	    }
	/*
	 * Else if we're in UE14, the A matrix is a row vector of
	 * 1x1 matrices and the B matrix is a row vector of c_rows x 1
	 * matrices. Divide each column in B by its respective A entry.
	 */
	} else if (c_type == VNACAL_E12 || VNACAL_IS_UE14(c_type)) {
	    for (int row = 0; row < c_ports; ++row) {
		for (int column = 0; column < c_ports; ++column) {
		    int cell = row * c_ports + column;

		    M(row, column) = vaa.vaa_b_matrix[cell][findex] /
				     vaa.vaa_a_matrix[column][findex];
		}
	    }
	/*
	 * Otherwise, find M = B A^-1.
	 */
	} else {
	    for (int cell = 0; cell < c_ports * c_ports; ++cell) {
		b[cell] = vaa.vaa_b_matrix[cell][findex];
		a[cell] = vaa.vaa_a_matrix[cell][findex];
	    }
	    determinant = _vnacommon_mrdivide(m, b, a, c_ports, c_ports);
	    if (determinant == 0.0 || !isfinite(cabs(determinant))) {
		_vnacal_error(vcp, VNAERR_MATH,
			"%s: 'a' matrix is singular at frequency index %d",
			vaa.vaa_function, findex);
		return -1;
	    }
	}

	/*
	 * Build a linear system of equations with coefficient matrix A
	 * and right-hand side matrix B to solve for the S-parameters.
	 * Though we're using the same storage, the A & B matrices here
	 * are unrelated to the A & B matrices above.
	 */
	switch (c_type) {
	case VNACAL_T8:
	case VNACAL_TE10:
	    fill_t8(&vl, t, m, a, b);
	    break;

	case VNACAL_U8:
	case VNACAL_UE10:
	    fill_u8(&vl, t, m, a, b);
	    break;

	case VNACAL_T16:
	    fill_t16(&vl, t, m, a, b);
	    break;

	case VNACAL_U16:
	    fill_u16(&vl, t, m, a, b);
	    break;

	case VNACAL_UE14:
	case _VNACAL_E12_UE14:
	    fill_ue14(&vl, t, m, a, b);
	    break;

	case VNACAL_E12:
	    fill_e12(&vl, t, m, a, b);
	    break;
	}

	/*
	 * Calculate S parameters.
	 *   In T:         S = A^-1 B
	 *   In U and E12: S = B A^-1
	 */
	determinant = 0.0;
	switch (c_type) {
	case VNACAL_T8:
	case VNACAL_TE10:
	case VNACAL_T16:
	    determinant = _vnacommon_mldivide(s, a, b, c_ports, c_ports);
	    break;

	case VNACAL_U8:
	case VNACAL_UE10:
	case VNACAL_U16:
	case VNACAL_UE14:
	case _VNACAL_E12_UE14:
	case VNACAL_E12:
	    determinant = _vnacommon_mrdivide(s, b, a, c_ports, c_ports);
	    break;
	}
	if (determinant == 0.0 || !isfinite(cabs(determinant))) {
	    _vnacal_error(vcp, VNAERR_MATH,
		    "%s: solution is singular at frequency index %d",
		    vaa.vaa_function, findex);
	    return -1;
	}

	/*
	 * Store the solution.
	 */
	for (int s_row = 0; s_row < c_ports; ++s_row) {
	    for (int s_column = 0; s_column < c_ports; ++s_column) {
		(void)vnadata_set_cell(vaa.vaa_s_parameters, findex,
			s_row, s_column, S(s_row, s_column));
	    }
	}
#undef S
#undef M
    }
    return 0;
}

/*
 * vnacal_apply: apply the calibration to measured values
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @frequency_vector: vector of increasing frequency points
 *   @frequencies: number of frequencies in matrix and s_parameters
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s_parameters: caller-allocated vnadata_t structure
 */
int vnacal_apply(vnacal_t *vcp, int ci,
	const double *frequency_vector, int frequencies,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	vnadata_t *s_parameters)
{
    vnacal_apply_args_t vaa;

    (void)memset((void *)&vaa, 0, sizeof(vaa));
    vaa.vaa_function		= __func__;
    vaa.vaa_vcp			= vcp;
    vaa.vaa_ci			= ci;
    vaa.vaa_frequency_vector	= frequency_vector;
    vaa.vaa_frequencies		= frequencies;
    vaa.vaa_a_matrix		= (const double complex *const *)a;
    vaa.vaa_a_rows		= a_rows;
    vaa.vaa_a_columns		= a_columns;
    vaa.vaa_b_matrix		= (const double complex *const *)b;
    vaa.vaa_b_rows		= b_rows;
    vaa.vaa_b_columns		= b_columns;
    vaa.vaa_m_type		= 'a';
    vaa.vaa_s_parameters	= s_parameters;

    return _vnacal_apply_common(vaa);
}

/*
 * vnacal_apply_m: apply the calibration to measured values
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @frequency_vector: vector of increasing frequency points
 *   @frequencies: number of frequencies in matrix and s_parameters
 *   @m: matrix of measured voltages
 *   @m_rows: number of rows in b
 *   @m_columns: number of columns in b
 *   @s_parameters: caller-allocated vnadata_t structure
 */
int vnacal_apply_m(vnacal_t *vcp, int ci,
	const double *frequency_vector, int frequencies,
	double complex *const *m, int m_rows, int m_columns,
	vnadata_t *s_parameters)
{
    vnacal_apply_args_t vaa;

    (void)memset((void *)&vaa, 0, sizeof(vaa));
    vaa.vaa_function		= __func__;
    vaa.vaa_vcp			= vcp;
    vaa.vaa_ci			= ci;
    vaa.vaa_frequency_vector	= frequency_vector;
    vaa.vaa_frequencies		= frequencies;
    vaa.vaa_a_matrix		= NULL;
    vaa.vaa_a_rows		= 0;
    vaa.vaa_a_columns		= 0;
    vaa.vaa_b_matrix		= (const double complex *const *)m;
    vaa.vaa_b_rows		= m_rows;
    vaa.vaa_b_columns		= m_columns;
    vaa.vaa_m_type		= 'm';
    vaa.vaa_s_parameters	= s_parameters;

    return _vnacal_apply_common(vaa);
}
