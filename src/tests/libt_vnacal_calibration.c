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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "libt.h"
#include "libt_vnacal.h"

/*
 * libt_vnacal_print_calibration: print solved calibration error terms
 *   @calp: pointer to calibration structure
 */
void libt_vnacal_print_calibration(vnacal_calibration_t *calp)
{
    vnacal_layout_t vl;

    (void)printf("calibration %s %d x %d",
	    _vnacal_type_to_name(calp->cal_type),
	    calp->cal_rows, calp->cal_columns);
    if (calp->cal_name != NULL) {
	(void)printf(" \"%s\":\n", calp->cal_name);
    } else {
	(void)printf(" (unnamed):\n");
    }
    _vnacal_layout(&vl, calp->cal_type, calp->cal_rows, calp->cal_columns);
    for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
	double complex **e = calp->cal_error_term_vector;

	(void)printf("f %e\n", calp->cal_frequency_vector[findex]);
	switch (VL_TYPE(&vl)) {
	case VNACAL_T8:
	case VNACAL_TE10:
	    {
		double complex **ts = &e[VL_TS_OFFSET(&vl)];
		double complex **ti = &e[VL_TI_OFFSET(&vl)];
		double complex **tx = &e[VL_TX_OFFSET(&vl)];
		double complex **tm = &e[VL_TM_OFFSET(&vl)];
		double complex **el = &e[VL_EL_OFFSET(&vl)];
		const int ts_terms = VL_TS_TERMS(&vl);
		const int ti_terms = VL_TI_TERMS(&vl);
		const int tx_terms = VL_TX_TERMS(&vl);
		const int tm_terms = VL_TM_TERMS(&vl);

		for (int i = 0; i < ts_terms; ++i) {
		    (void)printf("  ts%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(ts[i][findex]), cimag(ts[i][findex]));
		}
		for (int i = 0; i < ti_terms; ++i) {
		    (void)printf("  ti%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(ti[i][findex]), cimag(ti[i][findex]));
		}
		for (int i = 0; i < tx_terms; ++i) {
		    (void)printf("  tx%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(tx[i][findex]), cimag(tx[i][findex]));
		}
		for (int i = 0; i < tm_terms; ++i) {
		    (void)printf("  tm%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(tm[i][findex]), cimag(tm[i][findex]));
		}
		if (VL_TYPE(&vl) == VNACAL_TE10) {
		    const int el_rows    = VL_EL_ROWS(&vl);
		    const int el_columns = VL_EL_COLUMNS(&vl);
		    int term = 0;

		    for (int row = 0; row < el_rows; ++row) {
			for (int column = 0; column < el_columns; ++column) {
			    if (row != column) {
				(void)printf("  el%d%d: %8.5f%+8.5fj\n",
					row + 1, column + 1,
					creal(el[term][findex]),
					cimag(el[term][findex]));
				++term;
			    }
			}
		    }
		}
	    }
	    break;

	case VNACAL_U8:
	case VNACAL_UE10:
	    {
		double complex **um = &e[VL_UM_OFFSET(&vl)];
		double complex **ui = &e[VL_UI_OFFSET(&vl)];
		double complex **ux = &e[VL_UX_OFFSET(&vl)];
		double complex **us = &e[VL_US_OFFSET(&vl)];
		double complex **el = &e[VL_EL_OFFSET(&vl)];
		const int um_terms = VL_UM_TERMS(&vl);
		const int ui_terms = VL_UI_TERMS(&vl);
		const int ux_terms = VL_UX_TERMS(&vl);
		const int us_terms = VL_US_TERMS(&vl);

		for (int i = 0; i < um_terms; ++i) {
		    (void)printf("  um%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(um[i][findex]), cimag(um[i][findex]));
		}
		for (int i = 0; i < ui_terms; ++i) {
		    (void)printf("  ui%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(ui[i][findex]), cimag(ui[i][findex]));
		}
		for (int i = 0; i < ux_terms; ++i) {
		    (void)printf("  ux%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(ux[i][findex]), cimag(ux[i][findex]));
		}
		for (int i = 0; i < us_terms; ++i) {
		    (void)printf("  us%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1,
			    creal(us[i][findex]), cimag(us[i][findex]));
		}
		if (VL_TYPE(&vl) == VNACAL_TE10) {
		    const int el_rows    = VL_EL_ROWS(&vl);
		    const int el_columns = VL_EL_COLUMNS(&vl);
		    int term = 0;

		    for (int row = 0; row < el_rows; ++row) {
			for (int column = 0; column < el_columns; ++column) {
			    if (row != column) {
				(void)printf("  el%d%d: %8.5f%+8.5fj\n",
					row + 1, column + 1,
					creal(el[term][findex]),
					cimag(el[term][findex]));
				++term;
			    }
			}
		    }
		}
	    }
	    break;

	case VNACAL_T16:
	    {
		double complex **ts = &e[VL_TS_OFFSET(&vl)];
		double complex **ti = &e[VL_TI_OFFSET(&vl)];
		double complex **tx = &e[VL_TX_OFFSET(&vl)];
		double complex **tm = &e[VL_TM_OFFSET(&vl)];
		const int ts_rows    = VL_TS_ROWS(&vl);
		const int ts_columns = VL_TS_COLUMNS(&vl);
		const int ti_rows    = VL_TI_ROWS(&vl);
		const int ti_columns = VL_TI_COLUMNS(&vl);
		const int tx_rows    = VL_TX_ROWS(&vl);
		const int tx_columns = VL_TX_COLUMNS(&vl);
		const int tm_rows    = VL_TM_ROWS(&vl);
		const int tm_columns = VL_TM_COLUMNS(&vl);

		for (int row = 0; row < ts_rows; ++row) {
		    for (int column = 0; column < ts_columns; ++column) {
			int term = row * ts_columns + column;

			(void)printf("  ts%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ts[term][findex]),
				cimag(ts[term][findex]));
		    }
		}
		for (int row = 0; row < ti_rows; ++row) {
		    for (int column = 0; column < ti_columns; ++column) {
			int term = row * ti_columns + column;

			(void)printf("  ti%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ti[term][findex]),
				cimag(ti[term][findex]));
		    }
		}
		for (int row = 0; row < tx_rows; ++row) {
		    for (int column = 0; column < tx_columns; ++column) {
			int term = row * tx_columns + column;

			(void)printf("  tx%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(tx[term][findex]),
				cimag(tx[term][findex]));
		    }
		}
		for (int row = 0; row < tm_rows; ++row) {
		    for (int column = 0; column < tm_columns; ++column) {
			int term = row * tm_columns + column;

			(void)printf("  tm%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(tm[term][findex]),
				cimag(tm[term][findex]));
		    }
		}
	    }
	    break;

	case VNACAL_U16:
	    {
		double complex **um = &e[VL_US_OFFSET(&vl)];
		double complex **ui = &e[VL_UI_OFFSET(&vl)];
		double complex **ux = &e[VL_UX_OFFSET(&vl)];
		double complex **us = &e[VL_UM_OFFSET(&vl)];
		const int um_rows    = VL_US_ROWS(&vl);
		const int um_columns = VL_US_COLUMNS(&vl);
		const int ui_rows    = VL_UI_ROWS(&vl);
		const int ui_columns = VL_UI_COLUMNS(&vl);
		const int ux_rows    = VL_UX_ROWS(&vl);
		const int ux_columns = VL_UX_COLUMNS(&vl);
		const int us_rows    = VL_UM_ROWS(&vl);
		const int us_columns = VL_UM_COLUMNS(&vl);

		for (int row = 0; row < um_rows; ++row) {
		    for (int column = 0; column < um_columns; ++column) {
			int term = row * um_columns + column;

			(void)printf("  um%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(um[term][findex]),
				cimag(um[term][findex]));
		    }
		}
		for (int row = 0; row < ui_rows; ++row) {
		    for (int column = 0; column < ui_columns; ++column) {
			int term = row * ui_columns + column;

			(void)printf("  ui%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ui[term][findex]),
				cimag(ui[term][findex]));
		    }
		}
		for (int row = 0; row < ux_rows; ++row) {
		    for (int column = 0; column < ux_columns; ++column) {
			int term = row * ux_columns + column;

			(void)printf("  ux%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ux[term][findex]),
				cimag(ux[term][findex]));
		    }
		}
		for (int row = 0; row < us_rows; ++row) {
		    for (int column = 0; column < us_columns; ++column) {
			int term = row * us_columns + column;

			(void)printf("  us%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(us[term][findex]),
				cimag(us[term][findex]));
		    }
		}
	    }
	    break;

	case VNACAL_UE14:
	case _VNACAL_E12_UE14:
	    {
		const int m_columns = VL_M_COLUMNS(&vl);
		const int um_terms = VL_UM14_TERMS(&vl);
		const int ui_terms = VL_UI14_TERMS(&vl);
		const int ux_terms = VL_UX14_TERMS(&vl);
		const int us_terms = VL_US14_TERMS(&vl);
		const int el_rows    = VL_EL_ROWS(&vl);
		const int el_columns = VL_EL_COLUMNS(&vl);
		double complex **el = &e[VL_EL_OFFSET(&vl)];
		int term = 0;

		for (int m_column = 0; m_column < m_columns; ++m_column) {
		    double complex **um = &e[VL_UM14_OFFSET(&vl, m_column)];
		    double complex **ui = &e[VL_UI14_OFFSET(&vl, m_column)];
		    double complex **ux = &e[VL_UX14_OFFSET(&vl, m_column)];
		    double complex **us = &e[VL_US14_OFFSET(&vl, m_column)];

		    (void)printf("  m_column %d\n", m_column);
		    for (int i = 0; i < um_terms; ++i) {
			(void)printf("    um%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1,
				creal(um[i][findex]), cimag(um[i][findex]));
		    }
		    for (int i = 0; i < ui_terms; ++i) {
			(void)printf("    ui%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1,
				creal(ui[i][findex]), cimag(ui[i][findex]));
		    }
		    for (int i = 0; i < ux_terms; ++i) {
			(void)printf("    ux%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1,
				creal(ux[i][findex]), cimag(ux[i][findex]));
		    }
		    for (int i = 0; i < us_terms; ++i) {
			(void)printf("    us%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1,
				creal(us[i][findex]), cimag(us[i][findex]));
		    }
		}
		for (int row = 0; row < el_rows; ++row) {
		    for (int column = 0; column < el_columns;
			    ++column) {
			if (row != column) {
			    (void)printf("  el%d%d: %8.5f%+8.5fj\n",
				    row + 1, column + 1,
				    creal(el[term][findex]),
				    cimag(el[term][findex]));
			    ++term;
			}
		    }
		}
	    }
	    break;

	case VNACAL_E12:
	    {
		const int m_columns  = VL_M_COLUMNS(&vl);
		const int el_terms   = VL_EL12_TERMS(&vl);
		const int er_terms   = VL_ER12_TERMS(&vl);
		const int em_terms   = VL_EM12_TERMS(&vl);

		for (int m_column = 0; m_column < m_columns; ++m_column) {
		    double complex **el = &e[VL_EL12_OFFSET(&vl, m_column)];
		    double complex **er = &e[VL_ER12_OFFSET(&vl, m_column)];
		    double complex **em = &e[VL_EM12_OFFSET(&vl, m_column)];

		    (void)printf("  m_column %d\n", m_column);
		    for (int term = 0; term < el_terms; ++term) {
			(void)printf("    el%d1: %8.5f%+8.5fj\n",
				term + 1,
				creal(el[term][findex]),
				cimag(el[term][findex]));
		    }
		    for (int term = 0; term < er_terms; ++term) {
			(void)printf("    er%d%d: %8.5f%+8.5fj\n",
				term + 1, term + 1,
				creal(er[term][findex]),
				cimag(er[term][findex]));
		    }
		    for (int term = 0; term < em_terms; ++term) {
			(void)printf("    em%d%d: %8.5f%+8.5fj\n",
				term + 1, term + 1,
				creal(em[term][findex]),
				cimag(em[term][findex]));
		    }
		}
	    }
	}
    }
    if (calp->cal_properties != NULL) {
	(void)printf("properties:\n");
	libt_vnacal_print_properties(calp->cal_properties, 1);
    }
    (void)printf("\n");
}

/*
 * libt_vnacal_validate_calibration: compare calculated error terms to actual
 *   @ttp: pointer to test error terms structure
 *   @calp: pointer to calibration structure
 */
int libt_vnacal_validate_calibration(const libt_vnacal_terms_t *ttp,
	vnacal_calibration_t *calp)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;

    if (calp == NULL) {
	vnacal_new_t *vnp = ttp->tt_vnp;
	assert(vnp != NULL);
	calp = vnp->vn_calibration;
    }
    assert(calp != NULL);
    if (opt_v >= 2) {
	libt_vnacal_print_calibration(calp);
    }
    if (calp->cal_error_terms != VL_ERROR_TERMS(vlp)) {
	(void)printf("cal_error_terms (%d) != vl_error_terms (%d)\n",
		calp->cal_error_terms, VL_ERROR_TERMS(vlp));
	return -1;
    }
    for (int findex = 0; findex < ttp->tt_frequencies; ++findex) {
	for (int term = 0; term < VL_ERROR_TERMS(vlp); ++term) {
	    if (!libt_isequal(calp->cal_error_term_vector[term][findex],
			ttp->tt_error_term_vector[findex][term])) {
		if (opt_a) {
		    assert(!"data miscompare");
		}
		return -1;
	    }
	}
    }
    return 0;
}
