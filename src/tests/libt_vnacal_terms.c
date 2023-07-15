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
 * term_info_t: parameters and other state used internally
 */
typedef struct term_info {
    /* complex random number generator of off-diagonal error terms */
    libt_crand_generator_t *ti_cg0;

    /* complex random number generator of diagonal error terms */
    libt_crand_generator_t *ti_cg1;

} term_info_t;

/*
 * crand0: generate complex random numbers for off-diagonal entries
 *         which can be zero
 */
static inline double complex crand0(const term_info_t *tip)
{
    return tip->ti_cg0->cg_crand(tip->ti_cg0);
}

/*
 * crand1: generate complex random numbers for diagaonal entries which
 *         must not be zero
 */
static inline double complex crand1(const term_info_t *tip)
{
    return tip->ti_cg1->cg_crand(tip->ti_cg1);
}

/*
 * gen_e_terms: generate random error terms
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error term vector
 */
static void gen_e_terms(const vnacal_layout_t *vlp, const term_info_t *tip,
	double complex *e)
{
    const int m_columns = VL_M_COLUMNS(vlp);

    switch (VL_TYPE(vlp)) {
    case VNACAL_T8:
    case VNACAL_TE10:
	{
	    double complex *ts = &e[VL_TS_OFFSET(vlp)];
	    double complex *ti = &e[VL_TI_OFFSET(vlp)];
	    double complex *tx = &e[VL_TX_OFFSET(vlp)];
	    double complex *tm = &e[VL_TM_OFFSET(vlp)];
	    double complex *el = &e[VL_EL_OFFSET(vlp)];
	    const int ts_terms = VL_TS_TERMS(vlp);
	    const int ti_terms = VL_TI_TERMS(vlp);
	    const int tx_terms = VL_TX_TERMS(vlp);
	    const int tm_terms = VL_TM_TERMS(vlp);
	    const int el_terms = VL_EL_TERMS(vlp);
	    const int unity_offset = _vl_unity_offset(vlp, 0);

	    assert(unity_offset == VL_TM_OFFSET(vlp));
	    for (int ts_term = 0; ts_term < ts_terms; ++ts_term) {
		ts[ts_term] = crand1(tip);
	    }
	    for (int ti_term = 0; ti_term < ti_terms; ++ti_term) {
		ti[ti_term] = crand0(tip);
	    }
	    for (int tx_term = 0; tx_term < tx_terms; ++tx_term) {
		tx[tx_term] = crand0(tip);
	    }
	    for (int tm_term = 0; tm_term < tm_terms; ++tm_term) {
		tm[tm_term] = tm_term == 0 ? 1.0 : crand1(tip);
	    }
	    for (int term = 0; term < el_terms; ++term) {
		el[term] = crand0(tip);
	    }
	}
	break;

    case VNACAL_U8:
    case VNACAL_UE10:
	{
	    const int um_terms = VL_UM_TERMS(vlp);
	    const int ui_terms = VL_UI_TERMS(vlp);
	    const int ux_terms = VL_UX_TERMS(vlp);
	    const int us_terms = VL_US_TERMS(vlp);
	    const int el_terms = VL_EL_TERMS(vlp);
	    double complex *um = &e[VL_UM_OFFSET(vlp)];
	    double complex *ui = &e[VL_UI_OFFSET(vlp)];
	    double complex *ux = &e[VL_UX_OFFSET(vlp)];
	    double complex *us = &e[VL_US_OFFSET(vlp)];
	    double complex *el = &e[VL_EL_OFFSET(vlp)];
	    const int unity_offset = _vl_unity_offset(vlp, 0);

	    assert(unity_offset == VL_UM_OFFSET(vlp));
	    for (int um_term = 0; um_term < um_terms; ++um_term) {
		um[um_term] = um_term == 0 ? 1.0 : crand1(tip);
	    }
	    for (int ui_term = 0; ui_term < ui_terms; ++ui_term) {
		ui[ui_term] = crand0(tip);
	    }
	    for (int ux_term = 0; ux_term < ux_terms; ++ux_term) {
		ux[ux_term] = crand0(tip);
	    }
	    for (int us_term = 0; us_term < us_terms; ++us_term) {
		us[us_term] = crand1(tip);
	    }
	    for (int term = 0; term < el_terms; ++term) {
		el[term] = crand0(tip);
	    }
	}
	break;

    case VNACAL_T16:
	{
	    double complex *ts = &e[VL_TS_OFFSET(vlp)];
	    double complex *ti = &e[VL_TI_OFFSET(vlp)];
	    double complex *tx = &e[VL_TX_OFFSET(vlp)];
	    double complex *tm = &e[VL_TM_OFFSET(vlp)];
	    const int ts_rows    = VL_TS_ROWS(vlp);
	    const int ts_columns = VL_TS_COLUMNS(vlp);
	    const int ti_rows    = VL_TI_ROWS(vlp);
	    const int ti_columns = VL_TI_COLUMNS(vlp);
	    const int tx_rows    = VL_TX_ROWS(vlp);
	    const int tx_columns = VL_TX_COLUMNS(vlp);
	    const int tm_rows    = VL_TM_ROWS(vlp);
	    const int tm_columns = VL_TM_COLUMNS(vlp);
	    const int unity_offset = _vl_unity_offset(vlp, 0);

	    assert(unity_offset == VL_TM_OFFSET(vlp));
	    for (int ts_row = 0; ts_row < ts_rows; ++ts_row) {
		for (int ts_column = 0; ts_column < ts_columns; ++ts_column) {
		    const int ts_cell = ts_row * ts_columns + ts_column;

		    if (ts_row == ts_column) {
			ts[ts_cell] = crand1(tip);
		    } else {
			ts[ts_cell] = crand0(tip);
		    }
		}
	    }
	    for (int ti_row = 0; ti_row < ti_rows; ++ti_row) {
		for (int ti_column = 0; ti_column < ti_columns; ++ti_column) {
		    const int ti_cell = ti_row * ti_columns + ti_column;

		    ti[ti_cell] = crand0(tip);
		}
	    }
	    for (int tx_row = 0; tx_row < tx_rows; ++tx_row) {
		for (int tx_column = 0; tx_column < tx_columns; ++tx_column) {
		    const int tx_cell = tx_row * tx_columns + tx_column;

		    tx[tx_cell] = crand0(tip);
		}
	    }
	    for (int tm_row = 0; tm_row < tm_rows; ++tm_row) {
		for (int tm_column = 0; tm_column < tm_columns;
			++tm_column) {
		    const int tm_cell = tm_row * tm_columns + tm_column;

		    if (tm_cell == 0) {
			tm[tm_cell] = 1.0;
		    } else if (tm_row == tm_column) {
			tm[tm_cell] = crand1(tip);
		    } else {
			tm[tm_cell] = crand0(tip);
		    }
		}
	    }
	}
	break;

    case VNACAL_U16:
	{
	    double complex *um = &e[VL_UM_OFFSET(vlp)];
	    double complex *ui = &e[VL_UI_OFFSET(vlp)];
	    double complex *ux = &e[VL_UX_OFFSET(vlp)];
	    double complex *us = &e[VL_US_OFFSET(vlp)];
	    const int um_rows    = VL_UM_ROWS(vlp);
	    const int um_columns = VL_UM_COLUMNS(vlp);
	    const int ui_rows    = VL_UI_ROWS(vlp);
	    const int ui_columns = VL_UI_COLUMNS(vlp);
	    const int ux_rows    = VL_UX_ROWS(vlp);
	    const int ux_columns = VL_UX_COLUMNS(vlp);
	    const int us_rows    = VL_US_ROWS(vlp);
	    const int us_columns = VL_US_COLUMNS(vlp);
	    const int unity_offset = _vl_unity_offset(vlp, 0);

	    assert(unity_offset == VL_UM_OFFSET(vlp));
	    for (int um_row = 0; um_row < um_rows; ++um_row) {
		for (int um_column = 0; um_column < um_columns; ++um_column) {
		    const int um_cell = um_row * um_columns + um_column;

		    if (um_cell == 0) {
			um[um_cell] = 1.0;
		    } else if (um_row == um_column) {
			um[um_cell] = crand1(tip);
		    } else {
			um[um_cell] = crand0(tip);
		    }
		}
	    }
	    for (int ui_row = 0; ui_row < ui_rows; ++ui_row) {
		for (int ui_column = 0; ui_column < ui_columns; ++ui_column) {
		    const int ui_cell = ui_row * ui_columns + ui_column;

		    ui[ui_cell] = crand0(tip);
		}
	    }
	    for (int ux_row = 0; ux_row < ux_rows; ++ux_row) {
		for (int ux_column = 0; ux_column < ux_columns; ++ux_column) {
		    const int ux_cell = ux_row * ux_columns + ux_column;

		    ux[ux_cell] = crand0(tip);
		}
	    }
	    for (int us_row = 0; us_row < us_rows; ++us_row) {
		for (int us_column = 0; us_column < us_columns; ++us_column) {
		    const int us_cell = us_row * us_columns + us_column;

		    if (us_row == us_column) {
			us[us_cell] = crand1(tip);
		    } else {
			us[us_cell] = crand0(tip);
		    }
		}
	    }
	}
	break;

    case VNACAL_UE14:
	{
	    const int um_terms = VL_UM14_TERMS(vlp);
	    const int ui_terms = VL_UI14_TERMS(vlp);
	    const int ux_terms = VL_UX14_TERMS(vlp);
	    const int us_terms = VL_US14_TERMS(vlp);
	    const int el_terms = VL_EL_TERMS(vlp);
	    double complex *el = &e[VL_EL_OFFSET(vlp)];

	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		double complex *um = &e[VL_UM14_OFFSET(vlp, m_column)];
		double complex *ui = &e[VL_UI14_OFFSET(vlp, m_column)];
		double complex *ux = &e[VL_UX14_OFFSET(vlp, m_column)];
		double complex *us = &e[VL_US14_OFFSET(vlp, m_column)];
		const int unity_offset = _vl_unity_offset(vlp, m_column);

		for (int um_term = 0; um_term < um_terms; ++um_term) {
		    if (um_term == unity_offset) {
			um[um_term] = 1.0;
		    } else {
			um[um_term] = crand1(tip);
		    }
		}
		for (int ui_term = 0; ui_term < ui_terms; ++ui_term) {
		    ui[ui_term] = crand0(tip);
		}
		for (int ux_term = 0; ux_term < ux_terms; ++ux_term) {
		    ux[ux_term] = crand0(tip);
		}
		for (int us_term = 0; us_term < us_terms; ++us_term) {
		    us[us_term] = crand1(tip);
		}
	    }
	    for (int term = 0; term < el_terms; ++term) {
		el[term] = libt_crandn();
	    }
	}
	break;

    case VNACAL_E12:
	{
	    const int el_terms = VL_EL12_TERMS(vlp);
	    const int er_terms = VL_ER12_TERMS(vlp);
	    const int em_terms = VL_EM12_TERMS(vlp);

	    for (int m_column = 0; m_column < m_columns; ++m_column) {
		double complex *el = &e[VL_EL12_OFFSET(vlp, m_column)];
		double complex *er = &e[VL_ER12_OFFSET(vlp, m_column)];
		double complex *em = &e[VL_EM12_OFFSET(vlp, m_column)];

		for (int el_term = 0; el_term < el_terms; ++el_term) {
		    el[el_term] = crand0(tip);
		}
		for (int er_term = 0; er_term < er_terms; ++er_term) {
		    er[er_term] = crand1(tip);
		}
		for (int em_term = 0; em_term < em_terms; ++em_term) {
		    em[em_term] = crand0(tip);
		}
	    }
	}
	break;

    default:
	abort();
    }
}

/*
 * libt_vnacal_generate_error_terms: generate random error terms
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of calibration frequencies
 *   @frequency_vector: optional vector of specific frequencies
 *   @flags: misc flags (see libt_vnacal.h)
 */
libt_vnacal_terms_t *libt_vnacal_generate_error_terms(vnacal_t *vcp,
	vnacal_type_t type, int m_rows, int m_columns, int frequencies,
	const double *frequency_vector, uint32_t flags)
{
    term_info_t ti;
    libt_vnacal_terms_t *ttp = NULL;
    vnacal_layout_t *vlp;

    /*
     * Init the term_info structure.
     */
    (void)memset((void *)&ti, 0, sizeof(ti));

    if (!(flags & LIBT_CLOSE_ETERMS)) {
	/*
	 * Create the complex random number generator for off-diagonal terms
	 * with magnitude distributed to the Rayleigh distribution such that
	 * median value is 1.0.
	 */
	ti.ti_cg0 = libt_crand_generator(0.0, LIBT_IRLOG4, 0.0, 1000.0,
	                                 0.0, 360.0);

	/*
	 * Create the complex random number generator for diagonal terms with
	 * magnitude distributed according to the truncated Rice distribution
	 * with minimum value 0.1 and median value of 1.0.
	 */
	ti.ti_cg1 = libt_crand_generator(0.857148, 0.5, 0.1, 1000.0,
	                                 0.0, 360.0);

    } else { /* temp hack for poor convergence in vnacal_new_solve_auto */
	/*
	 * Create the complex random number generator for off-diagonal terms
	 * with magnitude distributed to the Rayleigh distribution with sigma
	 * 0.1, max 0.2.
	 */
	ti.ti_cg0 = libt_crand_generator(0.0, 0.1, 0.0, 0.2,
	                                 0.0, 360.0);

	/*
	 * Create the complex random number generator for diagonal terms with
	 * magnitude distributed according to the truncated Rice distribution
	 * with magnitude in [0.8 .. 2.0] and angle with 30 degrees of real.
	 */
	ti.ti_cg1 = libt_crand_generator(0.857148, 0.5, 0.8, 2.0,
	                                 0.0, 30.0);
    }
     libt_crand_generator_t *libt_crand_generator(double nu, double sigma,
        double min, double max, double rotation, double angle);



    /*
     * Create the error terms structure.
     */
    if ((ttp = malloc(sizeof(libt_vnacal_terms_t))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	exit(99);
    }
    (void)memset((void *)ttp, 0, sizeof(*ttp));
    _vnacal_layout(&ttp->tt_layout, type, m_rows, m_columns);
    vlp = &ttp->tt_layout;
    ttp->tt_frequencies = frequencies;
    if ((ttp->tt_frequency_vector = calloc(frequencies,
		    sizeof(double))) == NULL) {
	libt_vnacal_free_error_terms(ttp);
	ttp = NULL;
	goto out;
    }
    if (frequency_vector != NULL) {
	(void)memcpy((void *)ttp->tt_frequency_vector,
		(void *)frequency_vector, frequencies * sizeof(double));
    } else if (frequencies == 1) {
	ttp->tt_frequency_vector[0] = 1.0e+9;
    } else if (frequencies == 2) {
	ttp->tt_frequency_vector[0] = 0.0;
	ttp->tt_frequency_vector[1] = 1.0e+9;
    } else if (flags & LIBT_GET_2_10_GHZ) {
	for (int i = 0; i < frequencies; ++i) {
	    ttp->tt_frequency_vector[i] = 10.0e+9 *
		(double)(i + 1) / (double)frequencies;
	}
    } else {
	ttp->tt_frequency_vector[0] = 0.0;
	for (int i = 1; i < frequencies; ++i) {
	    ttp->tt_frequency_vector[i] = pow(1.0e+9,
		(double)(i - 1) / (double)(frequencies - 2));
	}
    }
    if ((ttp->tt_error_term_vector = calloc(frequencies,
		    sizeof(double complex *))) == NULL) {
	libt_vnacal_free_error_terms(ttp);
	ttp = NULL;
	goto out;
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	double complex *clfp;

	if ((clfp = calloc(VL_ERROR_TERMS(vlp),
			sizeof(double complex))) == NULL) {
	    libt_vnacal_free_error_terms(ttp);
	    ttp = NULL;
	    goto out;
	}
	ttp->tt_error_term_vector[findex] = clfp;
	gen_e_terms(vlp, &ti, clfp);
    }

    /*
     * Allocate the new calibration structure and set frequencies.
     */
    if ((ttp->tt_vnp = vnacal_new_alloc(vcp, type, m_rows, m_columns,
		    frequencies)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_new_alloc: %s\n",
		progname, strerror(errno));
	libt_vnacal_free_error_terms(ttp);
	ttp = NULL;
	goto out;
    }
    if (vnacal_new_set_frequency_vector(ttp->tt_vnp,
		ttp->tt_frequency_vector) == -1) {
	(void)fprintf(stderr, "%s: vnacal_new_set_frequency_vector: %s\n",
		progname, strerror(errno));
	libt_vnacal_free_error_terms(ttp);
	ttp = NULL;
	goto out;
    }

    /*
     * If verbose, show the error terms.
     */
    if (opt_v >= 2) {
	libt_vnacal_print_error_terms(ttp);
    }
    /*FALLTHROUGH*/

out:
    free((void *)ti.ti_cg1);
    free((void *)ti.ti_cg0);
    return ttp;
}

/*
 * libt_vnacal_print_error_terms: show the generated error terms
 *   @ttp: pointer to test error terms structure
 */
void libt_vnacal_print_error_terms(const libt_vnacal_terms_t *ttp)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;

    (void)printf("error terms %s %d x %d frequencies %d:\n",
	    vnacal_type_to_name(VL_TYPE(vlp)),
	    VL_M_ROWS(vlp), VL_M_COLUMNS(vlp), ttp->tt_frequencies);
    for (int frequency = 0; frequency < ttp->tt_frequencies; ++frequency) {
	(void)printf("f %e\n", ttp->tt_frequency_vector[frequency]);
	double complex *e  = ttp->tt_error_term_vector[frequency];

	switch (VL_TYPE(vlp)) {
	case VNACAL_T8:
	case VNACAL_TE10:
	    {
		double complex *ts = &e[VL_TS_OFFSET(vlp)];
		double complex *ti = &e[VL_TI_OFFSET(vlp)];
		double complex *tx = &e[VL_TX_OFFSET(vlp)];
		double complex *tm = &e[VL_TM_OFFSET(vlp)];
		double complex *el = &e[VL_EL_OFFSET(vlp)];
		const int ts_terms = VL_TS_TERMS(vlp);
		const int ti_terms = VL_TI_TERMS(vlp);
		const int tx_terms = VL_TX_TERMS(vlp);
		const int tm_terms = VL_TM_TERMS(vlp);

		for (int i = 0; i < ts_terms; ++i) {
		    (void)printf("  ts%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(ts[i]), cimag(ts[i]));
		}
		for (int i = 0; i < ti_terms; ++i) {
		    (void)printf("  ti%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(ti[i]), cimag(ti[i]));
		}
		for (int i = 0; i < tx_terms; ++i) {
		    (void)printf("  tx%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(tx[i]), cimag(tx[i]));
		}
		for (int i = 0; i < tm_terms; ++i) {
		    (void)printf("  tm%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(tm[i]), cimag(tm[i]));
		}
		if (VL_TYPE(vlp) == VNACAL_TE10) {
		    const int el_rows    = VL_EL_ROWS(vlp);
		    const int el_columns = VL_EL_COLUMNS(vlp);
		    int term = 0;

		    for (int row = 0; row < el_rows; ++row) {
			for (int column = 0; column < el_columns; ++column) {
			    if (row != column) {
				(void)printf("  el%d%d: %8.5f%+8.5fj\n",
					row + 1, column + 1,
					creal(el[term]), cimag(el[term]));
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
		double complex *um = &e[VL_UM_OFFSET(vlp)];
		double complex *ui = &e[VL_UI_OFFSET(vlp)];
		double complex *ux = &e[VL_UX_OFFSET(vlp)];
		double complex *us = &e[VL_US_OFFSET(vlp)];
		double complex *el = &e[VL_EL_OFFSET(vlp)];
		const int um_terms = VL_UM_TERMS(vlp);
		const int ui_terms = VL_UI_TERMS(vlp);
		const int ux_terms = VL_UX_TERMS(vlp);
		const int us_terms = VL_US_TERMS(vlp);

		for (int i = 0; i < um_terms; ++i) {
		    (void)printf("  um%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(um[i]), cimag(um[i]));
		}
		for (int i = 0; i < ui_terms; ++i) {
		    (void)printf("  ui%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(ui[i]), cimag(ui[i]));
		}
		for (int i = 0; i < ux_terms; ++i) {
		    (void)printf("  ux%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(ux[i]), cimag(ux[i]));
		}
		for (int i = 0; i < us_terms; ++i) {
		    (void)printf("  us%d%d: %8.5f%+8.5fj\n",
			    i + 1, i + 1, creal(us[i]), cimag(us[i]));
		}
		if (VL_TYPE(vlp) == VNACAL_UE10) {
		    const int el_rows    = VL_EL_ROWS(vlp);
		    const int el_columns = VL_EL_COLUMNS(vlp);
		    int term = 0;

		    for (int row = 0; row < el_rows; ++row) {
			for (int column = 0; column < el_columns; ++column) {
			    if (row != column) {
				(void)printf("  el%d%d: %8.5f%+8.5fj\n",
					row + 1, column + 1,
					creal(el[term]), cimag(el[term]));
				++term;
			    }
			}
		    }
		}
	    }
	    break;

	case VNACAL_T16:
	    {
		double complex *ts = &e[VL_TS_OFFSET(vlp)];
		double complex *ti = &e[VL_TI_OFFSET(vlp)];
		double complex *tx = &e[VL_TX_OFFSET(vlp)];
		double complex *tm = &e[VL_TM_OFFSET(vlp)];
		const int ts_rows    = VL_TS_ROWS(vlp);
		const int ts_columns = VL_TS_COLUMNS(vlp);
		const int ti_rows    = VL_TI_ROWS(vlp);
		const int ti_columns = VL_TI_COLUMNS(vlp);
		const int tx_rows    = VL_TX_ROWS(vlp);
		const int tx_columns = VL_TX_COLUMNS(vlp);
		const int tm_rows    = VL_TM_ROWS(vlp);
		const int tm_columns = VL_TM_COLUMNS(vlp);

		for (int row = 0; row < ts_rows; ++row) {
		    for (int column = 0; column < ts_columns; ++column) {
			int term = row * ts_columns + column;

			(void)printf("  ts%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ts[term]), cimag(ts[term]));
		    }
		}
		for (int row = 0; row < ti_rows; ++row) {
		    for (int column = 0; column < ti_columns; ++column) {
			int term = row * ti_columns + column;

			(void)printf("  ti%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ti[term]), cimag(ti[term]));
		    }
		}
		for (int row = 0; row < tx_rows; ++row) {
		    for (int column = 0; column < tx_columns; ++column) {
			int term = row * tx_columns + column;

			(void)printf("  tx%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(tx[term]), cimag(tx[term]));
		    }
		}
		for (int row = 0; row < tm_rows; ++row) {
		    for (int column = 0; column < tm_columns; ++column) {
			int term = row * tm_columns + column;

			(void)printf("  tm%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(tm[term]), cimag(tm[term]));
		    }
		}
	    }
	    break;

	case VNACAL_U16:
	    {
		double complex *um = &e[VL_UM_OFFSET(vlp)];
		double complex *ui = &e[VL_UI_OFFSET(vlp)];
		double complex *ux = &e[VL_UX_OFFSET(vlp)];
		double complex *us = &e[VL_US_OFFSET(vlp)];
		const int um_rows    = VL_UM_ROWS(vlp);
		const int um_columns = VL_UM_COLUMNS(vlp);
		const int ui_rows    = VL_UI_ROWS(vlp);
		const int ui_columns = VL_UI_COLUMNS(vlp);
		const int ux_rows    = VL_UX_ROWS(vlp);
		const int ux_columns = VL_UX_COLUMNS(vlp);
		const int us_rows    = VL_US_ROWS(vlp);
		const int us_columns = VL_US_COLUMNS(vlp);

		for (int row = 0; row < um_rows; ++row) {
		    for (int column = 0; column < um_columns; ++column) {
			int term = row * um_columns + column;

			(void)printf("  um%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(um[term]), cimag(um[term]));
		    }
		}
		for (int row = 0; row < ui_rows; ++row) {
		    for (int column = 0; column < ui_columns; ++column) {
			int term = row * ui_columns + column;

			(void)printf("  ui%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ui[term]), cimag(ui[term]));
		    }
		}
		for (int row = 0; row < ux_rows; ++row) {
		    for (int column = 0; column < ux_columns; ++column) {
			int term = row * ux_columns + column;

			(void)printf("  ux%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(ux[term]), cimag(ux[term]));
		    }
		}
		for (int row = 0; row < us_rows; ++row) {
		    for (int column = 0; column < us_columns; ++column) {
			int term = row * us_columns + column;

			(void)printf("  us%d%d: %8.5f%+8.5fj\n",
				row + 1, column + 1,
				creal(us[term]), cimag(us[term]));
		    }
		}
	    }
	    break;

	case VNACAL_UE14:
	case _VNACAL_E12_UE14:
	    {
		const int m_columns = VL_M_COLUMNS(vlp);
		const int um_terms = VL_UM14_TERMS(vlp);
		const int ui_terms = VL_UI14_TERMS(vlp);
		const int ux_terms = VL_UX14_TERMS(vlp);
		const int us_terms = VL_US14_TERMS(vlp);
		const int el_rows    = VL_EL_ROWS(vlp);
		const int el_columns = VL_EL_COLUMNS(vlp);
		double complex *el = &e[VL_EL_OFFSET(vlp)];
		int term = 0;

		for (int m_column = 0; m_column < m_columns; ++m_column) {
		    double complex *um = &e[VL_UM14_OFFSET(vlp, m_column)];
		    double complex *ui = &e[VL_UI14_OFFSET(vlp, m_column)];
		    double complex *ux = &e[VL_UX14_OFFSET(vlp, m_column)];
		    double complex *us = &e[VL_US14_OFFSET(vlp, m_column)];

		    (void)printf("  m_column %d\n", m_column);
		    for (int i = 0; i < um_terms; ++i) {
			(void)printf("    um%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1, creal(um[i]), cimag(um[i]));
		    }
		    for (int i = 0; i < ui_terms; ++i) {
			(void)printf("    ui%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1, creal(ui[i]), cimag(ui[i]));
		    }
		    for (int i = 0; i < ux_terms; ++i) {
			(void)printf("    ux%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1, creal(ux[i]), cimag(ux[i]));
		    }
		    for (int i = 0; i < us_terms; ++i) {
			(void)printf("    us%d%d: %8.5f%+8.5fj\n",
				i + 1, i + 1, creal(us[i]), cimag(us[i]));
		    }
		}
		for (int row = 0; row < el_rows; ++row) {
		    for (int column = 0; column < el_columns;
			    ++column) {
			if (row != column) {
			    (void)printf("  el%d%d: %8.5f%+8.5fj\n",
				    row + 1, column + 1,
				    creal(el[term]), cimag(el[term]));
			    ++term;
			}
		    }
		}
	    }
	    break;

	case VNACAL_E12:
	    {
		const int m_columns  = VL_M_COLUMNS(vlp);
		const int el_terms = VL_EL12_TERMS(vlp);
		const int er_terms = VL_ER12_TERMS(vlp);
		const int em_terms = VL_EM12_TERMS(vlp);

		for (int m_column = 0; m_column < m_columns; ++m_column) {
		    double complex *el = &e[VL_EL12_OFFSET(vlp, m_column)];
		    double complex *er = &e[VL_ER12_OFFSET(vlp, m_column)];
		    double complex *em = &e[VL_EM12_OFFSET(vlp, m_column)];

		    (void)printf("  m_column %d\n", m_column);
		    for (int term = 0; term < el_terms; ++term) {
			(void)printf("    el%d1: %8.5f%+8.5fj\n",
				term + 1,
				creal(el[term]), cimag(el[term]));
		    }
		    for (int term = 0; term < er_terms; ++term) {
			(void)printf("    er%d%d: %8.5f%+8.5fj\n",
				term + 1, term + 1,
				creal(er[term]), cimag(er[term]));
		    }
		    for (int term = 0; term < em_terms; ++term) {
			(void)printf("    em%d%d: %8.5f%+8.5fj\n",
				term + 1, term + 1,
				creal(em[term]), cimag(em[term]));
		    }
		}
	    }
	    break;

	default:
	    abort();
	}
    }
    (void)printf("\n");
}

/*
 * libt_vnacal_free_error_terms: free test error terms
 *   @ttp: pointer to test error terms structure
 */
void libt_vnacal_free_error_terms(libt_vnacal_terms_t *ttp)
{
    if (ttp != NULL) {
	vnacal_new_free(ttp->tt_vnp);
	if (ttp->tt_error_term_vector != NULL) {
	    for (int findex = ttp->tt_frequencies - 1; findex >= 0; --findex) {
		free((void *)ttp->tt_error_term_vector[findex]);
	    }
	    free((void *)ttp->tt_error_term_vector);
	}
	free((void *)ttp->tt_frequency_vector);
	free((void *)ttp);
    }
}
