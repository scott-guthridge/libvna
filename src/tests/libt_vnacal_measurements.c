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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "libt.h"
#include "libt_crand.h"
#include "libt_vnacal.h"

/*
 * libt_vnacal_sigma_n: per-freq vector of noise to add to measurements
 */
const double *libt_vnacal_sigma_n = NULL;

/*
 * libt_vnacal_sigma_t: per-freq vector of error to add to measurements
 */
const double *libt_vnacal_sigma_t = NULL;

/*
 * libt_vnacal_alloc_measurements: allocate test measurements
 *   @type: error term type
 *   @m_rows: number of rows in the measurement matrix
 *   @m_columns: number of columns in the measurement matrix
 *   @frequencies: number of frequency points
 *   @ab: true: use a, b matrices; false: use m matrix
 */
libt_vnacal_measurements_t *libt_vnacal_alloc_measurements(vnacal_type_t type,
	int m_rows, int m_columns, int frequencies, bool ab)
{
    libt_vnacal_measurements_t *tmp = NULL;

    if ((tmp = malloc(sizeof(libt_vnacal_measurements_t))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	exit(99);
    }
    (void)memset((void *)tmp, 0, sizeof(*tmp));
    if (ab) {
	const int a_rows = type != VNACAL_E12 && !VNACAL_IS_UE14(type) ?
					m_columns : 1;
	const int a_columns = m_columns;

	if ((tmp->tm_a_matrix = calloc(a_rows * a_columns,
			sizeof(double complex *))) == NULL) {
	    (void)fprintf(stderr, "%s: calloc: %s\n", progname,
		    strerror(errno));
	    libt_vnacal_free_measurements(tmp);
	    exit(99);
	}
	for (int a_cell = 0; a_cell < a_rows * a_columns; ++a_cell) {
	    if ((tmp->tm_a_matrix[a_cell] = calloc(frequencies,
			    sizeof(double complex))) == NULL) {
		(void)fprintf(stderr, "%s: calloc: %s\n", progname,
			strerror(errno));
		libt_vnacal_free_measurements(tmp);
		exit(99);
	    }
	}
	tmp->tm_a_rows    = a_rows;
	tmp->tm_a_columns = a_columns;
    }
    if ((tmp->tm_b_matrix = calloc(m_rows * m_columns,
		    sizeof(double complex *))) == NULL) {
	(void)fprintf(stderr, "%s: calloc: %s\n", progname, strerror(errno));
	libt_vnacal_free_measurements(tmp);
	exit(99);
    }
    for (int b_cell = 0; b_cell < m_rows * m_columns; ++b_cell) {
	if ((tmp->tm_b_matrix[b_cell] = calloc(frequencies,
			sizeof(double complex))) == NULL) {
	    (void)fprintf(stderr, "%s: calloc: %s\n", progname,
		    strerror(errno));
	    libt_vnacal_free_measurements(tmp);
	    exit(99);
	}
    }
    tmp->tm_b_rows    = m_rows;
    tmp->tm_b_columns = m_columns;

    return tmp;
}

/*
 * calc_m: calculate measurements given a full S matrix and error terms
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error term vector
 *   @s: s-parameter matrix s_rows x s_columns
 *   @m: measurement matrix m_rows x m_columns
 */
static int calc_m(const vnacal_layout_t *vlp, const double complex *e,
	const double complex *s, double complex *m)
{
    const int m_rows    = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int s_rows    = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);

    /*
     * Solve for M
     */
    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	{
	    const double complex *ts = &e[VL_TS_OFFSET(vlp)];/* mr x sr diag */
	    const double complex *ti = &e[VL_TI_OFFSET(vlp)];/* mr x sc diag */
	    const double complex *tx = &e[VL_TX_OFFSET(vlp)];/* mc x sr diag */
	    const double complex *tm = &e[VL_TM_OFFSET(vlp)];/* mc x sc diag */
	    const int ts_rows    = VL_TS_ROWS(vlp);
	    const int ts_columns = VL_TS_COLUMNS(vlp);
	    const int ti_rows    = VL_TI_ROWS(vlp);
	    const int ti_columns = VL_TI_COLUMNS(vlp);
	    const int tx_rows    = VL_TX_ROWS(vlp);
	    const int tx_columns = VL_TX_COLUMNS(vlp);
	    const int tm_rows    = VL_TM_ROWS(vlp);
	    const int tm_columns = VL_TM_COLUMNS(vlp);
	    double complex a[tm_rows * tm_columns]; /* mc x sc */
	    double complex b[ti_rows * ti_columns]; /* mr x sc */
	    double complex determinant;

	    assert(ts_rows    == m_rows);	/* by definition */
	    assert(ts_columns == s_rows);	/* by definition */
	    assert(ti_rows    == m_rows);	/* by definition */
	    assert(ti_columns == s_columns);	/* by definition */
	    assert(tx_rows    == m_columns);	/* by definition */
	    assert(tx_columns == s_rows);	/* by definition */
	    assert(tm_rows    == m_columns);	/* by definition */
	    assert(tm_columns == s_columns);	/* by definition */
	    assert(tm_rows    == tm_columns);	/* Tm must be square */
	    assert(m_columns  == s_columns);	/* M's must span all columns */
	    for (int a_row = 0; a_row < tm_rows; ++a_row) {
		for (int a_column = 0; a_column < tm_columns; ++a_column) {
		    const int a_cell = a_row * tm_columns + a_column;

		    a[a_cell] = 0.0;
		    if (a_row < s_rows) {
			const int s_cell = a_row * s_columns + a_column;

			a[a_cell] = tx[a_row] * s[s_cell];
		    }
		    if (a_row == a_column) {
			a[a_cell] += tm[a_row];
		    }
		}
	    }
	    for (int b_row = 0; b_row < ti_rows; ++b_row) {
		for (int b_column = 0; b_column < ti_columns; ++b_column) {
		    const int b_cell = b_row * ti_columns + b_column;

		    b[b_cell] = 0.0;
		    if (b_row < s_rows) {
			const int s_cell = b_row * s_columns + b_column;

			b[b_cell] = ts[b_row] * s[s_cell];
		    }
		    if (b_row == b_column) {
			b[b_cell] += ti[b_row];
		    }
		}
	    }
	    determinant = _vnacommon_mrdivide(m, b, a, m_rows, m_columns);
	    if (determinant == 0.0) {
		errno = EDOM;
		return -1;
	    }
	}
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
	{
	    const double complex *um = &e[VL_UM_OFFSET(vlp)]; /* sr x mr */
	    const double complex *ui = &e[VL_UI_OFFSET(vlp)]; /* sr x mc */
	    const double complex *ux = &e[VL_UX_OFFSET(vlp)]; /* sc x mr */
	    const double complex *us = &e[VL_US_OFFSET(vlp)]; /* sc x mc */
	    const int um_rows    = VL_UM_ROWS(vlp);
	    const int um_columns = VL_UM_COLUMNS(vlp);
	    const int ui_rows    = VL_UI_ROWS(vlp);
	    const int ui_columns = VL_UI_COLUMNS(vlp);
	    const int ux_rows    = VL_UX_ROWS(vlp);
	    const int ux_columns = VL_UX_COLUMNS(vlp);
	    const int us_rows    = VL_US_ROWS(vlp);
	    const int us_columns = VL_US_COLUMNS(vlp);
	    double complex a[um_rows * um_columns]; /* sr x mr */
	    double complex b[ui_rows * ui_columns]; /* sr x mc */
	    double complex determinant;

	    assert(um_rows    == s_rows);	/* by definition */
	    assert(um_columns == m_rows);	/* by definition */
	    assert(ui_rows    == s_rows);	/* by definition */
	    assert(ui_columns == m_columns);	/* by definition */
	    assert(ux_rows    == s_columns);	/* by definition */
	    assert(ux_columns == m_rows);	/* by definition */
	    assert(us_rows    == s_columns);	/* by definition */
	    assert(us_columns == m_columns);	/* by definition */
	    assert(um_rows    == um_columns);	/* Um must be square */
	    assert(m_rows     == s_rows);	/* M's must span all rows */
	    for (int a_row = 0; a_row < um_rows; ++a_row) {
		for (int a_column = 0; a_column < um_columns; ++a_column) {
		    const int a_cell = a_row * um_columns + a_column;

		    a[a_cell] = 0.0;
		    if (a_row == a_column) {
			a[a_cell] = um[a_row];
		    }
		    if (a_column < s_columns) {
			const int s_cell = a_row * s_columns + a_column;

			a[a_cell] -= s[s_cell] * ux[a_column];
		    }
		}
	    }
	    for (int b_row = 0; b_row < ui_rows; ++b_row) {
		for (int b_column = 0; b_column < ui_columns; ++b_column) {
		    const int b_cell = b_row * ui_columns + b_column;

		    b[b_cell] = 0.0;
		    if (b_column < s_columns) {
			const int s_cell = b_row * s_columns + b_column;

			b[b_cell] = us[b_column] * s[s_cell];
		    }
		    if (b_row == b_column) {
			b[b_cell] -= ui[b_row];
		    }
		}
	    }
	    determinant = _vnacommon_mldivide(m, a, b, m_rows, m_columns);
	    if (determinant == 0.0) {
		errno = EDOM;
		return -1;
	    }
	}
	break;

    case VNACAL_T16:
	{
	    const double complex *ts = &e[VL_TS_OFFSET(vlp)]; /* mr x sr */
	    const double complex *ti = &e[VL_TI_OFFSET(vlp)]; /* mr x sc */
	    const double complex *tx = &e[VL_TX_OFFSET(vlp)]; /* mc x sr */
	    const double complex *tm = &e[VL_TM_OFFSET(vlp)]; /* mc x sc */
	    const int ts_rows    = VL_TS_ROWS(vlp);
	    const int ts_columns = VL_TS_COLUMNS(vlp);
	    const int ti_rows    = VL_TI_ROWS(vlp);
	    const int ti_columns = VL_TI_COLUMNS(vlp);
	    const int tx_rows    = VL_TX_ROWS(vlp);
	    const int tx_columns = VL_TX_COLUMNS(vlp);
	    const int tm_rows    = VL_TM_ROWS(vlp);
	    const int tm_columns = VL_TM_COLUMNS(vlp);
	    double complex a[tm_rows * tm_columns]; /* mc x sc */
	    double complex b[ti_rows * ti_columns]; /* mr x sc */
	    double complex determinant;

	    assert(ts_rows    == m_rows);	/* by definition */
	    assert(ts_columns == s_rows);	/* by definition */
	    assert(ti_rows    == m_rows);	/* by definition */
	    assert(ti_columns == s_columns);	/* by definition */
	    assert(tx_rows    == m_columns);	/* by definition */
	    assert(tx_columns == s_rows);	/* by definition */
	    assert(tm_rows    == m_columns);	/* by definition */
	    assert(tm_columns == s_columns);	/* by definition */
	    assert(tm_rows    == tm_columns);	/* Tm must be square */
	    assert(m_columns  == s_columns);	/* M's must span all columns */
	    for (int a_row = 0; a_row < tm_rows; ++a_row) {
		for (int a_column = 0; a_column < tm_columns; ++a_column) {
		    const int a_cell = a_row * tm_columns + a_column;

		    a[a_cell] = 0.0;
		    for (int s_row = 0; s_row < s_rows; ++s_row) {
			const int tx_cell = a_row * s_rows + s_row;
			const int s_cell = s_row * s_columns + a_column;

			a[a_cell] += tx[tx_cell] * s[s_cell];
		    }
		    a[a_cell] += tm[a_cell];
		}
	    }
	    for (int b_row = 0; b_row < ti_rows; ++b_row) {
		for (int b_column = 0; b_column < ti_columns; ++b_column) {
		    const int b_cell = b_row * ti_columns + b_column;

		    b[b_cell] = 0.0;
		    for (int s_row = 0; s_row < s_rows; ++s_row) {
			const int ts_cell = b_row * ts_columns + s_row;
			const int s_cell = s_row * s_columns + b_column;

			b[b_cell] += ts[ts_cell] * s[s_cell];
		    }
		    b[b_cell] += ti[b_cell];
		}
	    }
	    determinant = _vnacommon_mrdivide(m, b, a, m_rows, m_columns);
	    if (determinant == 0.0 || !isnormal(cabs(determinant))) {
		errno = EDOM;
		return -1;
	    }
	}
	break;

    case VNACAL_U16:
	{
	    const double complex *um = &e[VL_UM_OFFSET(vlp)]; /* sr x mr */
	    const double complex *ui = &e[VL_UI_OFFSET(vlp)]; /* sr x mc */
	    const double complex *ux = &e[VL_UX_OFFSET(vlp)]; /* sc x mr */
	    const double complex *us = &e[VL_US_OFFSET(vlp)]; /* sc x mc */
	    const int um_rows    = VL_UM_ROWS(vlp);
	    const int um_columns = VL_UM_COLUMNS(vlp);
	    const int ui_rows    = VL_UI_ROWS(vlp);
	    const int ui_columns = VL_UI_COLUMNS(vlp);
	    const int ux_rows    = VL_UX_ROWS(vlp);
	    const int ux_columns = VL_UX_COLUMNS(vlp);
	    const int us_rows    = VL_US_ROWS(vlp);
	    const int us_columns = VL_US_COLUMNS(vlp);
	    double complex a[um_rows * um_columns]; /* sr x mr */
	    double complex b[ui_rows * ui_columns]; /* sr x mc */
	    double complex determinant;

	    assert(um_rows    == s_rows);	/* by definition */
	    assert(um_columns == m_rows);	/* by definition */
	    assert(ui_rows    == s_rows);	/* by definition */
	    assert(ui_columns == m_columns);	/* by definition */
	    assert(ux_rows    == s_columns);	/* by definition */
	    assert(ux_columns == m_rows);	/* by definition */
	    assert(us_rows    == s_columns);	/* by definition */
	    assert(us_columns == m_columns);	/* by definition */
	    assert(um_rows    == um_columns);	/* Um must be square */
	    assert(m_rows     == s_rows);	/* M's must span all rows */
	    for (int a_row = 0; a_row < um_rows; ++a_row) {
		for (int a_column = 0; a_column < um_columns; ++a_column) {
		    const int a_cell = a_row * um_columns + a_column;

		    a[a_cell] = um[a_cell];
		    for (int s_column = 0; s_column < s_columns; ++s_column) {
			const int ux_cell = s_column * ux_columns + a_column;
			const int s_cell = a_row * s_columns + s_column;

			a[a_cell] -= s[s_cell] * ux[ux_cell];
		    }
		}
	    }
	    for (int b_row = 0; b_row < ui_rows; ++b_row) {
		for (int b_column = 0; b_column < ui_columns; ++b_column) {
		    const int b_cell = b_row * ui_columns + b_column;

		    b[b_cell] = 0.0;
		    for (int s_column = 0; s_column < s_columns; ++s_column) {
			const int us_cell = s_column * us_columns + b_column;
			const int s_cell = b_row * s_columns + s_column;

			b[b_cell] += us[us_cell] * s[s_cell];
		    }
		    b[b_cell] -= ui[b_cell];
		}
	    }
	    determinant = _vnacommon_mldivide(m, a, b, m_rows, m_columns);
	    if (determinant == 0.0) {
		errno = EDOM;
		return -1;
	    }
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
	assert(m_rows == s_rows);	/* M's must span all rows */
	for (int m_column = 0; m_column < m_columns; ++m_column) {
	    const double complex *um = &e[VL_UM14_OFFSET(vlp, m_column)];
	    const double complex *ui = &e[VL_UI14_OFFSET(vlp, m_column)];
	    const double complex *ux = &e[VL_UX14_OFFSET(vlp, m_column)];
	    const double complex *us = &e[VL_US14_OFFSET(vlp, m_column)];
	    const int um_rows    = VL_UM14_ROWS(vlp);
	    const int um_columns = VL_UM14_COLUMNS(vlp);
	    const int ui_rows    = VL_UI14_ROWS(vlp);
	    const int ui_columns = VL_UI14_COLUMNS(vlp);
	    const int ux_rows    = VL_UX14_ROWS(vlp);
	    const int ux_columns = VL_UX14_COLUMNS(vlp);
	    const int us_rows    = VL_US14_ROWS(vlp);
	    const int us_columns = VL_US14_COLUMNS(vlp);
	    double complex a[um_rows * um_columns]; /* sr x mr */
	    double complex b[ui_rows * 1];	    /* sr x 1 */
	    double complex x[s_rows  * 1];
	    double complex determinant;

	    assert(um_rows    == s_rows);	/* definition */
	    assert(um_columns == m_rows);	/* definition */
	    assert(ui_rows    == s_rows);	/* definition */
	    assert(ui_columns == 1);		/* definition */
	    assert(ux_rows    == s_columns);	/* definition */
	    assert(ux_columns == m_rows);	/* definition */
	    assert(us_rows    == s_columns);	/* definition */
	    assert(us_columns == 1);		/* definition */
	    assert(um_rows    == um_columns);	/* Um must be square */
	    for (int a_row = 0; a_row < um_rows; ++a_row) {
		for (int a_column = 0; a_column < um_columns; ++a_column) {
		    const int a_cell = a_row * um_columns + a_column;

		    a[a_cell] = 0.0;
		    if (a_row == a_column) {
			a[a_cell] = um[a_row];
		    }
		    if (a_column < s_columns) {
			const int s_cell = a_row * s_columns + a_column;

			a[a_cell] -= s[s_cell] * ux[a_column];
		    }
		}
	    }
	    for (int b_row = 0; b_row < ui_rows; ++b_row) {
		b[b_row] = 0.0;
		if (m_column < s_columns) {
		    const int s_cell = b_row * s_columns + m_column;

		    b[b_row] = us[0] * s[s_cell];
		}
		if (b_row == m_column) {
		    b[b_row] -= ui[0];
		}
	    }
	    determinant = _vnacommon_mldivide(x, a, b, m_rows, 1);
	    if (determinant == 0.0) {
		errno = EDOM;
		return -1;
	    }
	    for (int m_row = 0; m_row < m_rows; ++m_row) {
		int m_cell = m_row * m_columns + m_column;

		m[m_cell] = x[m_row];
	    }
	}
	break;

    case VNACAL_E12:
	{
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		const double complex *el = &e[VL_EL12_OFFSET(vlp, m_column)];
		const double complex *er = &e[VL_ER12_OFFSET(vlp, m_column)];
		const double complex *em = &e[VL_EM12_OFFSET(vlp, m_column)];
		double complex a[s_columns * s_columns];
		double complex b[m_rows * s_columns];
		double complex x[m_rows * s_columns];
		double complex determinant;

		/*
		 * A = I - Em S
		 */
		for (int a_row = 0; a_row < s_columns; ++a_row) {
		    for (int a_column = 0; a_column < s_columns; ++a_column) {
			const int a_cell = a_row * s_columns + a_column;

			a[a_cell] = a_row == a_column ? 1.0 : 0.0;
			if (a_row < s_rows) {
			    a[a_cell] -= em[a_row] * s[a_cell];
			}
		    }
		}

		/*
		 * B = Er S
		 */
		for (int b_row = 0; b_row < m_rows; ++b_row) {
		    for (int b_column = 0; b_column < s_columns; ++b_column) {
			const int b_cell = b_row * s_columns + b_column;

			b[b_cell] = 0.0;
			if (b_row < s_rows) {
			    b[b_cell] = er[b_row] * s[b_cell];
			}
		    }
		}

		/*
		 * X = B A^-1 = Er S (I - Em S)^-1
		 */
		determinant = _vnacommon_mrdivide(x, b, a, m_rows, s_columns);
		if (determinant == 0.0) {
		    errno = EDOM;
		    return -1;
		}

		/*
		 * M(:, m_column) = El + Er S (I - Em S)^-1 Et
		 *   where Et is the m_column'th column of the identify matrix
		 */
		for (int m_row = 0; m_row < m_rows; ++m_row) {
		    const int m_cell = m_row * m_columns + m_column;
		    const int x_cell = m_row * s_columns + m_column;

		    m[m_cell] = el[m_row] + x[x_cell];
		}
	    }
	}
	break;

    default:
	abort();
    }

    /*
     * If we have leakage terms handled outside of the linear system,
     * add them here.
     */
    if (VL_TYPE(vlp) == VNACAL_TE10 || VL_TYPE(vlp) == VNACAL_UE10 ||
	    VL_IS_UE14(vlp)) {
	const double complex *el_cur = &e[VL_EL_OFFSET(vlp)];

	for (int m_row = 0; m_row < m_rows; ++m_row) {
	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		if (m_row != m_column) {
		    int m_cell = m_row * m_columns + m_column;

		    m[m_cell] += *el_cur++;
		}
	    }
	}
	assert(el_cur == &e[VL_EL_OFFSET(vlp) + VL_EL_TERMS(vlp)]);
    }
    return 0;
}

/*
 * calc_measurements_helper: form s matrix and calc m matrix
 *   @ttp: pointer to test error terms structure
 *   @tmp: caller-allocated libt_vnacal_measurements structure to hold result
 *   @s_matrix: s-parameter indices matrix describing the standard
 *   @s_matrix_rows: rows in s_matrix
 *   @s_matrix_columns: columns in s_matrix
 *   @port_map: map from standard port to VNA port
 *   @findex: frequency index
 *   @m: caller-allocated output matrix
 */
static int calc_measurements_helper(const libt_vnacal_terms_t *ttp,
	libt_vnacal_measurements_t *tmp, const int *s_matrix,
	int s_matrix_rows, int s_matrix_columns, const int *port_map,
	int findex, double complex *m)
{
    vnacal_new_t *vnp = ttp->tt_vnp;
    vnacal_t *vcp = vnp->vn_vcp;
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int s_rows = VL_S_ROWS(vlp);
    const int s_columns = VL_S_COLUMNS(vlp);
    double f = ttp->tt_frequency_vector[findex];
    double complex s[s_rows * s_columns];
    bool port_used[MAX(s_rows, s_columns)];
    bool cell_defined[s_rows * s_columns];

    /*
     * Fill in the full S matrix following the example below.  Here,
     * the VNA has five ports and the standard has three ports, which are
     * connected to VNA ports 2, 3, and 4, respectively, numbering from 1.
     * Let small s11, s12, s13, etc. refer to the ports of the standard,
     * and capital S11, S12, S13, etc. refer to the ports of the VNA.
     * The given standard matrix is rectangular indicating that s13, s23
     * and s33 parameters are unknown.	We don't know anything about the
     * s-parameters between VNA ports 1 and 5, but we assume they remain
     * consistent over the measurement and that they have no external
     * connection to ports 2, 3, 4.  Because of this, we know that S12,
     * S13, S14, S21, S25, etc. are zero.      All remaining cells of s
     * (marked with "*") are unknown, and to reflect that, we fill them
     * with random values.
     *
     *	*   0	0   0	*
     *	0   s11 s12 *	0
     *	0   s21 s22 *	0
     *	0   s31 s32 *	0
     *	*   0	0   0	*
     */
    (void)memset((void *)port_used, 0, sizeof(port_used));
    (void)memset((void *)cell_defined, 0, sizeof(cell_defined));
    for (int r = 0; r < s_matrix_rows; ++r) {
	for (int c = 0; c < s_matrix_columns; ++c) {
	    int s_row =    port_map != NULL ? port_map[r] - 1 : r;
	    int s_column = port_map != NULL ? port_map[c] - 1 : c;
	    int s_matrix_cell = r * s_matrix_columns + c;
	    int s_cell = s_row * s_columns + s_column;

	    assert(s_row >= 0 && s_row < s_rows);
	    assert(s_column >= 0 && s_column < s_columns);
	    vnacal_parameter_t *vpmrp = _vnacal_get_parameter(vcp,
		    s_matrix[s_matrix_cell]);

	    if (vpmrp == NULL) {
		(void)fprintf(stderr, "%s: _vnacal_get_parameter: %s\n",
			progname, strerror(errno));
		return -1;
	    }
	    /*
	     * TODO: handle 'unknown' parameters here
	     *
	     * First time a particular unknown is seen, create an
	     * entry in a new tt_unknown_parameters vector and fill it
	     * with an random value.  Use that value in the s matrix.
	     * If the same unknown is seen again on subsequent calls,
	     * use the same value as before.
	     */
	    s[s_cell] = _vnacal_get_parameter_value_i(vpmrp, f);
	    port_used[s_row] = true;
	    port_used[s_column] = true;
	    cell_defined[s_cell] = true;
	}
    }
    for (int s_row = 0; s_row < s_rows; ++s_row) {
	for (int s_column = 0; s_column < s_columns; ++s_column) {
	    int s_cell = s_row * s_columns + s_column;

	    if ((port_used[s_row] && !port_used[s_column]) ||
		(!port_used[s_row] && port_used[s_column])) {
		s[s_cell] = 0.0;
		cell_defined[s_cell] = true;
	    }
	}
    }
    for (int s_cell = 0; s_cell < s_rows * s_columns; ++s_cell) {
	if (!cell_defined[s_cell]) {
	    s[s_cell] = libt_crandn();
	}
    }

    /*
     * Calculate M.
     */
    if (calc_m(vlp, ttp->tt_error_term_vector[findex], s, m) == -1) {
	return -1;
    }
    return 0;
}

/*
 * libt_vnacal_calculate_measurements: calculate measurements of standard
 *   @ttp: pointer to test error terms structure
 *   @tmp: caller-allocated libt_vnacal_measurements structure to hold result
 *   @s_matrix: s-parameter indices matrix describing the standard
 *   @s_matrix_rows: rows in s_matrix
 *   @s_matrix_columns: columns in s_matrix
 *   @port_map: map from standard port to VNA port
 */
int libt_vnacal_calculate_measurements(const libt_vnacal_terms_t *ttp,
	libt_vnacal_measurements_t *tmp,
	const int *s_matrix, int s_matrix_rows, int s_matrix_columns,
	const int *port_map)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    const int b_rows = tmp->tm_b_rows;
    const int b_columns = tmp->tm_b_columns;
    const int frequencies = ttp->tt_frequencies;
    double complex **a_matrix = tmp->tm_a_matrix;
    double complex **b_matrix = tmp->tm_b_matrix;

    /*
     * If verbose, show the standard.
     */
    if (opt_v >= 2) {
	vnacal_new_t *vnp = ttp->tt_vnp;
	vnacal_t *vcp = vnp->vn_vcp;

	libt_vnacal_print_standard(vcp, s_matrix, s_matrix_rows,
		s_matrix_columns, ttp->tt_frequencies,
		ttp->tt_frequency_vector, port_map);
    }

    /*
     * For each frequency...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double complex m[b_rows * b_columns];

	/*
	 * Handle the normal case where the output M matrix has the
	 * same dimensions of the calibration M matrix.
	 */
	if (b_rows == m_rows && b_columns == m_columns) {
	    if (calc_measurements_helper(ttp, tmp, s_matrix,
			s_matrix_rows, s_matrix_columns, port_map, findex,
			m) == -1) {
		return -1;
	    }

	/*
	 * Handle the special case where the output M matrix is 2x2
	 * but the calibration matrix is either 1x2 or 2x1.
	 */
	} else {
	    int temp_map[2];
	    double complex temp_m[2];

	    /*
	     * Calculate and place the first vector.
	     */
	    assert(b_rows == 2 && b_columns == 2);
	    assert(m_rows * m_columns == 2);
	    if (calc_measurements_helper(ttp, tmp, s_matrix,
			s_matrix_rows, s_matrix_columns, port_map, findex,
			temp_m) == -1) {
		return -1;
	    }
	    m[0] = temp_m[0];
	    if (m_rows == 1) {
		m[1] = temp_m[1];
	    } else {
		m[2] = temp_m[1];
	    }
	    /*
	     * Swap the ports and calculate the second vector.  We also
	     * have to swap the resulting M values.
	     */
	    if (port_map != NULL) {
		temp_map[0] = port_map[1];
		temp_map[1] = port_map[0];
	    } else {
		temp_map[0] = 2;
		temp_map[1] = 1;
	    }
	    if (calc_measurements_helper(ttp, tmp, s_matrix,
			s_matrix_rows, s_matrix_columns, temp_map, findex,
			temp_m) == -1) {
		return -1;
	    }
	    if (m_rows == 1) {
		m[2] = temp_m[1];
	    } else {
		m[1] = temp_m[1];
	    }
	    m[3] = temp_m[0];
	}

	/*
	 * Add per-frequency random measurement error if configured.
	 */
	if (libt_vnacal_sigma_t != NULL) {
	    for (int cell = 0; cell < b_rows * b_columns; ++cell) {
		m[cell] += libt_vnacal_sigma_t[findex] * m[cell] *
		    libt_crandn();
	    }
	}
	if (libt_vnacal_sigma_n != NULL) {
	    for (int cell = 0; cell < b_rows * b_columns; ++cell) {
		m[cell] += libt_vnacal_sigma_n[findex] * libt_crandn();
	    }
	}

	/*
	 * If an A matrix was given, fill it with random values and
	 * calculate B = m * A;
	 */
	if (a_matrix == NULL) {
	    for (int cell = 0; cell < b_rows * b_columns; ++cell) {
		b_matrix[cell][findex] = m[cell];
	    }
	} else if (VL_HAS_COLUMN_SYSTEMS(vlp)) {
	    for (int b_column = 0; b_column < b_columns; ++b_column) {
		double complex a = libt_crandn();

		a_matrix[b_column][findex] = a;
		for (int m_row = 0; m_row < b_rows; ++m_row) {
		    int cell = m_row * b_columns + b_column;

		    b_matrix[cell][findex] = m[cell] * a;
		}
	    }
	} else {
	    double complex a[b_columns * b_columns];
	    double complex b[b_rows * b_columns];

	    for (int a_cell = 0; a_cell < b_columns * b_columns; ++a_cell) {
		a[a_cell] = libt_crandn();
		a_matrix[a_cell][findex] = a[a_cell];
	    }
	    _vnacommon_mmultiply(b, m, a, b_rows, b_columns, b_columns);
	    for (int b_cell = 0; b_cell < b_rows * b_columns; ++b_cell) {
		b_matrix[b_cell][findex] = b[b_cell];
	    }
	}
    }

    /*
     * If verbose, show values.
     */
    if (opt_v >= 2) {
	libt_vnacal_print_measurements(tmp, frequencies);
    }

    return 0;
}

/*
 * libt_vnacal_print_measurements: print the "measured" values
 *   @tmp: test measurements structure
 *   @frequencies: number of frequencies
 */
void libt_vnacal_print_measurements(libt_vnacal_measurements_t *tmp,
	int frequencies)
{
    (void)printf("measurements %d x %d:\n",
	    tmp->tm_b_rows, tmp->tm_b_columns);
    for (int findex = 0; findex < frequencies; ++findex) {
	(void)printf("findex %d\n", findex);
	if (tmp->tm_a_matrix != NULL) {
	    for (int row = 0; row < tmp->tm_a_rows; ++row) {
		for (int column = 0; column < tmp->tm_a_columns; ++column) {
		    int cell = row * tmp->tm_a_columns + column;

		    (void)printf("  a%d%d: %8.5f%+8.5fj\n",
			row + 1, column + 1,
			creal(tmp->tm_a_matrix[cell][findex]),
			cimag(tmp->tm_a_matrix[cell][findex]));
		}
	    }
	}
	for (int row = 0; row < tmp->tm_b_rows; ++row) {
	    for (int column = 0; column < tmp->tm_b_columns; ++column) {
		int cell = row * tmp->tm_b_columns + column;

		(void)printf("  %c%d%d: %8.5f%+8.5fj\n",
			tmp->tm_a_matrix == NULL ? 'm' : 'b',
			row + 1, column + 1,
			creal(tmp->tm_b_matrix[cell][findex]),
			cimag(tmp->tm_b_matrix[cell][findex]));
	    }
	}
    }
    (void)printf("\n");
}

/*
 * libt_vnacal_free_measurements: free a libt_vnacal_measurements_t structure
 *   @tmp: test measurements structure
 */
void libt_vnacal_free_measurements(libt_vnacal_measurements_t *tmp)
{
    if (tmp != NULL) {
	if (tmp->tm_a_matrix != NULL) {
	    for (int i = 0; i < tmp->tm_a_rows * tmp->tm_a_columns; ++i) {
		free((void *)tmp->tm_a_matrix[i]);
	    }
	    free((void *)tmp->tm_a_matrix);
	}
	if (tmp->tm_b_matrix != NULL) {
	    for (int i = 0; i < tmp->tm_b_rows * tmp->tm_b_columns; ++i) {
		free((void *)tmp->tm_b_matrix[i]);
	    }
	    free((void *)tmp->tm_b_matrix);
	}
	free((void *)tmp);
    }
}
