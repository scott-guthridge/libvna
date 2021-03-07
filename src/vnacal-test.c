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
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <unistd.h>
#include "vnacal_internal.h"


char *progname;

#define PI	3.1415926535897932384626433832795
#define EPS	1.0e-4

#define NTRIALS		67

/*
 * test_result_type
 */
typedef enum test_result {
    T_PASS,
    T_FAIL,
    T_SKIPPED
} test_result_type;

/*
 * Test Counters
 */
static int test_count = 0;
static int fail_count = 0;

/*
 * Command Line Options
 */
static const char options[] = "av";
static const char *const usage[] = {
    "[-av]",
    NULL
};
static const char *const help[] = {
    "-a	 abort on data miscompare",
    "-v	 show verbose output",
    NULL
};
static bool opt_a = false;
static int opt_v = 0;

/*
 * crandn: generate a random complex number where real and imaginary parts
 *	are normally distributed with zero mean and unit standard deviation
 */
static double complex crandn()
{
    double u1 = (random() + 1.0) / RAND_MAX;
    double u2 = (double)random() / RAND_MAX;
    double r = sqrt(-2.0 * log(u1));
    double a = 2 * PI * u2;

    return r * (cos(a) + I * sin(a));
}

/*
 * isequal: test if x and y are approximately equal
 */
static int isequal(double complex x, double complex y)
{
    int rv;
    double d = cabs(csqrt(x * y));

    if (d < 1.0) {
	d = 1.0;
    }
    rv = cabs(x - y) / d < EPS;
    if (!rv) {
	printf("|x-y| = %f\n", cabs(x - y));
	printf("%f%+fi != %f%+fi\n", creal(x), cimag(x), creal(y), cimag(y));
    }
    return rv;
}

/*
 * cmatrix_multiply: find C = A x B
 *   @c: serialized result matrix, m x o
 *   @a: serialized A matrix, m x n
 *   @b: serialized B matrix, n x o
 *   @m: first dimension of C and A
 *   @n: second dimension of A, first dimension of B
 *   @o: second dimension of C and B
 */
static void cmatrix_multiply(double complex *c, const double complex *a,
	const double complex *b, int m, int n, int o)
{
#define A(i, j) (a[(i) * n + (j)])
#define B(i, j) (b[(i) * o + (j)])
#define C(i, j) (c[(i) * o + (j)])

    for (int i = 0; i < m; ++i) {
	for (int k = 0; k < o; ++k) {
	    double complex s = 0.0;

	    for (int j = 0; j < n; ++j) {
		s += A(i, j) * B(j, k);
	    }
	    C(i, k) = s;
	}
    }
}
#undef A
#undef B
#undef C

/*
 * error_fn: error reporting function
 *   @category: error category
 *   @message: error message
 *   @arg: (unused)
 */
static void error_fn(vnaerr_category_t category, const char *message, void *arg)
{
    (void)fprintf(stderr, "%s: %s\n", progname, message);
}

/*
 * report_test_result: report a test result
 */
static void report_test_result(const char *test_name, test_result_type result)
{
    const char *result_name;

    switch (result) {
    case T_PASS:
	result_name = "PASS";
	break;
    case T_FAIL:
	result_name = "FAIL";
	break;
    case T_SKIPPED:
	result_name = "SKIPPED";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)printf("Test %2d: %-58s %s\n", ++test_count, test_name, result_name);
    (void)fflush(stdout);
    if (result == T_FAIL) {
	++fail_count;
    }
}

/*
 * test_terms_t: expected error terms
 */
typedef struct test_terms {
    /* error term type and layout */
    vnacal_layout_t tt_layout;

    /* number and vector of test frequencies */
    int tt_frequencies;
    double *tt_frequency_vector;

    /* vector (one per frequency) of vectors of error terms */
    double complex **tt_error_term_vector;

    /* associated vnacal_new_t structure, if not NULL */
    vnacal_new_t *tt_vnp;

} test_terms_t;

/*
 * test_measurements_t: measurement matrices
 */
typedef struct test_measurements {
    double complex **tm_a_matrix;
    double complex **tm_b_matrix;
    int tm_a_rows;
    int tm_a_columns;
    int tm_b_rows;
    int tm_b_columns;
} test_measurements_t;

/*
 * print_test_error_terms: show the test error terms
 *   @ttp: pointer to test error terms structure
 */
static void print_test_error_terms(const test_terms_t *ttp)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;

    (void)printf("error terms %s %d x %d frequencies %d:\n",
	    _vnacal_type_to_name(VL_TYPE(vlp)),
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
		double complex *um = &e[VL_US_OFFSET(vlp)];
		double complex *ui = &e[VL_UI_OFFSET(vlp)];
		double complex *ux = &e[VL_UX_OFFSET(vlp)];
		double complex *us = &e[VL_UM_OFFSET(vlp)];
		const int um_rows    = VL_US_ROWS(vlp);
		const int um_columns = VL_US_COLUMNS(vlp);
		const int ui_rows    = VL_UI_ROWS(vlp);
		const int ui_columns = VL_UI_COLUMNS(vlp);
		const int ux_rows    = VL_UX_ROWS(vlp);
		const int ux_columns = VL_UX_COLUMNS(vlp);
		const int us_rows    = VL_UM_ROWS(vlp);
		const int us_columns = VL_UM_COLUMNS(vlp);

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
	}
    }
    (void)printf("\n");
}

/*
 * print_standard: show a calibration standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @s: s-parameter indices matrix describing the standard
 *   @s_rows: rows in s_matrix
 *   @s_columns: columns in s_matrix
 *   @frequencies: number of calibration frequencies
 *   @frequency_vector: vector of frequencies
 *   @port_map: map from standard port to VNA port
 */
static void print_standard(vnacal_t *vcp, const int *s,
	int s_rows, int s_columns,
	int frequencies, const double *frequency_vector,
	const int *port_map)
{
    bool has_vector = false;

    /*
     * First scan to determine if any standards are of vector type.
     */
    for (int row = 0; row < s_rows; ++row) {
	for (int column = 0; column < s_columns; ++column) {
	    int cell = row * s_columns + column;
	    vnacal_parameter_t *vpmrp;

	    vpmrp = _vnacal_get_parameter(vcp, s[cell]);
	    for (;;) {
		switch (vpmrp->vpmr_type) {
		case VNACAL_NEW:
		    abort();
		case VNACAL_SCALAR:
		    break;
		case VNACAL_VECTOR:
		    has_vector = true;
		    break;
		case VNACAL_UNKNOWN:
		case VNACAL_CORRELATED:
		    vpmrp = vpmrp->vpmr_other;
		    continue;
		}
		break;
	    }
	}
    }

    /*
     * Print
     */
    (void)printf("standard %d x %d:\n", s_rows, s_columns);
    if (has_vector) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    double f = frequency_vector[findex];

	    (void)printf("f %e\n", frequency_vector[findex]);
	    for (int row = 0; row < s_rows; ++row) {
		for (int column = 0; column < s_columns; ++column) {
		    int cell = row * s_columns + column;
		    vnacal_parameter_t *vpmrp;
		    double complex value;

		    vpmrp = _vnacal_get_parameter(vcp, s[cell]);
		    value = _vnacal_get_parameter_value(vpmrp, f);
		    (void)printf("  s%d%d: %8.5f%+8.5fj\n",
			row + 1, column + 1, creal(value), cimag(value));
		}
	    }
	}
    } else {
	for (int row = 0; row < s_rows; ++row) {
	    for (int column = 0; column < s_columns; ++column) {
		int cell = row * s_columns + column;
		vnacal_parameter_t *vpmrp;
		double complex value;

		vpmrp = _vnacal_get_parameter(vcp, s[cell]);
		value = _vnacal_get_parameter_value(vpmrp, 0.0);
		(void)printf("  s%d%d: %8.5f%+8.5fj\n",
			row + 1, column + 1, creal(value), cimag(value));
	    }
	}
    }
    if (port_map != NULL) {
	int ports = MAX(s_rows, s_columns);

	(void)printf("map:");
	for (int port = 0; port < ports; ++port) {
	    (void)printf(" %d", port_map[port]);
	}
	(void)printf("\n");
    }
    (void)printf("\n");
}

/*
 * print_test_measurements: print the "measured" values
 *   @tmp: test measurements structure
 *   @frequencies: number of frequencies
 */
static void print_test_measurements(test_measurements_t *tmp, int frequencies)
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
 * print_properties: print a property list
 */
static void print_properties(const vnaproperty_t *vprp, int indent)
{
    if (vprp == NULL) {
	for (int i = 0; i < indent; ++i) {
	    (void)printf("    ");
	}
	printf(".\n");
	return;
    }
    switch (vnaproperty_type(vprp)) {
    case VNAPROPERTY_SCALAR:
	for (int i = 0; i < indent; ++i) {
	    (void)printf("    ");
	}
	(void)printf("\"%s\"\n", vnaproperty_scalar_get(vprp));
	return;

    case VNAPROPERTY_MAP:
	{
	    const vnaproperty_map_pair_t *vmprp;

	    for (vmprp = vnaproperty_map_begin(vprp); vmprp != NULL;
		    vmprp = vnaproperty_map_next(vmprp)) {
		for (int i = 0; i < indent; ++i) {
		    (void)printf("    ");
		}
		(void)printf(".%s\n", vmprp->vmpr_key);
		print_properties(vmprp->vmpr_value, indent + 1);
	    }
	}
	return;

    case VNAPROPERTY_LIST:
	{
	    int count = vnaproperty_list_count(vprp);

	    for (int i = 0; i < count; ++i) {
		for (int i = 0; i < indent; ++i) {
		    (void)printf("    ");
		}
		(void)printf("[%d]\n", i);
		print_properties(vnaproperty_list_get(vprp, i), indent + 1);
	    }
	}
	return;

    default:
	abort();
    }
}

/*
 * print_calibration: print solved calibration error terms
 *   @calp: pointer to calibration structure
 */
static void print_calibration(vnacal_calibration_t *calp)
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
	print_properties(calp->cal_properties, 1);
    }
    (void)printf("\n");
}

/*
 * free_test_measurements: free a test_measurements_t structure
 *   @tmp: test measurements structure
 */
static void free_test_measurements(test_measurements_t *tmp)
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

/*
 * alloc_test_measurements: allocate test measurements of the given dimensions
 *   @type: error term type
 *   @m_rows: number of rows in the measurement matrix
 *   @m_columns: number of columns in the measurement matrix
 *   @frequencies: number of frequency points
 *   @ab: true: use a, b matrices; false: use m matrix
 */
static test_measurements_t *alloc_test_measurements(vnacal_type_t type,
	int m_rows, int m_columns, int frequencies, bool ab)
{
    test_measurements_t *tmp = NULL;

    if ((tmp = malloc(sizeof(test_measurements_t))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	return NULL;
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
	    free_test_measurements(tmp);
	    return NULL;
	}
	for (int a_cell = 0; a_cell < a_rows * a_columns; ++a_cell) {
	    if ((tmp->tm_a_matrix[a_cell] = calloc(frequencies,
			    sizeof(double complex))) == NULL) {
		(void)fprintf(stderr, "%s: calloc: %s\n", progname,
			strerror(errno));
		free_test_measurements(tmp);
		return NULL;
	    }
	}
	tmp->tm_a_rows    = a_rows;
	tmp->tm_a_columns = a_columns;
    }
    if ((tmp->tm_b_matrix = calloc(m_rows * m_columns,
		    sizeof(double complex *))) == NULL) {
	(void)fprintf(stderr, "%s: calloc: %s\n", progname, strerror(errno));
	free_test_measurements(tmp);
	return NULL;
    }
    for (int b_cell = 0; b_cell < m_rows * m_columns; ++b_cell) {
	if ((tmp->tm_b_matrix[b_cell] = calloc(frequencies,
			sizeof(double complex))) == NULL) {
	    (void)fprintf(stderr, "%s: calloc: %s\n", progname,
		    strerror(errno));
	    free_test_measurements(tmp);
	    return NULL;
	}
    }
    tmp->tm_b_rows    = m_rows;
    tmp->tm_b_columns = m_columns;

    return tmp;
}

/*
 * free_test_terms: free the memory allocated in gen_test_terms
 *   @ttp: pointer to test error terms structure
 */
static void free_test_terms(test_terms_t *ttp)
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

/*
 * gen_e_terms: generate random error terms
 *   @vlp: pointer to vnacal_layout_t structure
 *   @e: error term vector
 *   @sigma: standard deviation of error (if 0.0, VNA is perfect)
 */
static void gen_e_terms(const vnacal_layout_t *vlp, double complex *e,
	double sigma)
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
		ts[ts_term] = 1.0;
		if (sigma != 0.0) {
		    ts[ts_term] += sigma * crandn();
		}
	    }
	    for (int ti_term = 0; ti_term < ti_terms; ++ti_term) {
		ti[ti_term] = 0.0;
		if (sigma != 0.0) {
		    ti[ti_term] += sigma * crandn();
		}
	    }
	    for (int tx_term = 0; tx_term < tx_terms; ++tx_term) {
		tx[tx_term] = 0.0;
		if (sigma != 0.0) {
		    tx[tx_term] += sigma * crandn();
		}
	    }
	    for (int tm_term = 0; tm_term < tm_terms; ++tm_term) {
		tm[tm_term] = 1.0;
		if (sigma != 0.0 && tm_term != 0) {
		    tm[tm_term] += sigma * crandn();
		}
	    }
	    for (int term = 0; term < el_terms; ++term) {
		el[term] = crandn();
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

	    for (int um_term = 0; um_term < um_terms; ++um_term) {
		um[um_term] = 1.0;
		if (sigma != 0.0 && um_term != unity_offset) {
		    um[um_term] += sigma * crandn();
		}
	    }
	    for (int ui_term = 0; ui_term < ui_terms; ++ui_term) {
		ui[ui_term] = 0.0;
		if (sigma != 0.0) {
		    ui[ui_term] += sigma * crandn();
		}
	    }
	    for (int ux_term = 0; ux_term < ux_terms; ++ux_term) {
		ux[ux_term] = 0.0;
		if (sigma != 0.0) {
		    ux[ux_term] += sigma * crandn();
		}
	    }
	    for (int us_term = 0; us_term < us_terms; ++us_term) {
		us[us_term] = 1.0;
		if (sigma != 0.0) {
		    us[us_term] += sigma * crandn();
		}
	    }
	    for (int term = 0; term < el_terms; ++term) {
		el[term] = crandn();
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

		    ts[ts_cell] = (ts_row == ts_column) ? 1.0 : 0.0;
		    if (sigma != 0.0) {
			ts[ts_cell] += sigma * crandn();
		    }
		}
	    }
	    for (int ti_row = 0; ti_row < ti_rows; ++ti_row) {
		for (int ti_column = 0; ti_column < ti_columns; ++ti_column) {
		    const int ti_cell = ti_row * ti_columns + ti_column;

		    ti[ti_cell] = 0.0;
		    if (sigma != 0.0) {
			ti[ti_cell] += sigma * crandn();
		    }
		}
	    }
	    for (int tx_row = 0; tx_row < tx_rows; ++tx_row) {
		for (int tx_column = 0; tx_column < tx_columns; ++tx_column) {
		    const int tx_cell = tx_row * tx_columns + tx_column;

		    tx[tx_cell] = 0.0;
		    if (sigma != 0.0) {
			tx[tx_cell] += sigma * crandn();
		    }
		}
	    }
	    for (int tm_row = 0; tm_row < tm_rows; ++tm_row) {
		for (int tm_column = 0; tm_column < tm_columns;
			++tm_column) {
		    const int tm_cell = tm_row * tm_columns + tm_column;

		    tm[tm_cell] = (tm_row == tm_column) ? 1.0 : 0.0;
		    if (sigma != 0.0 && tm_cell != 0) {
			tm[tm_cell] += sigma * crandn();
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

		    um[um_cell] = (um_row == um_column) ? 1.0 : 0.0;
		    if (sigma != 0.0 && um_cell != 0) {
			um[um_cell] += sigma * crandn();
		    }
		}
	    }
	    for (int ui_row = 0; ui_row < ui_rows; ++ui_row) {
		for (int ui_column = 0; ui_column < ui_columns; ++ui_column) {
		    const int ui_cell = ui_row * ui_columns + ui_column;

		    ui[ui_cell] = 0.0;
		    if (sigma != 0.0) {
			ui[ui_cell] += sigma * crandn();
		    }
		}
	    }
	    for (int ux_row = 0; ux_row < ux_rows; ++ux_row) {
		for (int ux_column = 0; ux_column < ux_columns; ++ux_column) {
		    const int ux_cell = ux_row * ux_columns + ux_column;

		    ux[ux_cell] = 0.0;
		    if (sigma != 0.0) {
			ux[ux_cell] += sigma * crandn();
		    }
		}
	    }
	    for (int us_row = 0; us_row < us_rows; ++us_row) {
		for (int us_column = 0; us_column < us_columns; ++us_column) {
		    const int us_cell = us_row * us_columns + us_column;

		    us[us_cell] = (us_row == us_column) ? 1.0 : 0.0;
		    if (sigma != 0.0) {
			us[us_cell] += sigma * crandn();
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
		    um[um_term] = 1.0;
		    if (sigma != 0.0 && um_term != unity_offset) {
			um[um_term] += sigma * crandn();
		    }
		}
		for (int ui_term = 0; ui_term < ui_terms; ++ui_term) {
		    ui[ui_term] = 0.0;
		    if (sigma != 0.0) {
			ui[ui_term] += sigma * crandn();
		    }
		}
		for (int ux_term = 0; ux_term < ux_terms; ++ux_term) {
		    ux[ux_term] = 0.0;
		    if (sigma != 0.0) {
			ux[ux_term] += sigma * crandn();
		    }
		}
		for (int us_term = 0; us_term < us_terms; ++us_term) {
		    us[us_term] = 1.0;
		    if (sigma != 0.0) {
			us[us_term] += sigma * crandn();
		    }
		}
	    }
	    for (int term = 0; term < el_terms; ++term) {
		el[term] = crandn();
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
		    el[el_term] = 0.0;
		    if (sigma != 0.0) {
			el[el_term] += sigma * crandn();
		    }
		}
		for (int er_term = 0; er_term < er_terms; ++er_term) {
		    er[er_term] = 1.0;
		    if (sigma != 0.0) {
			er[er_term] += sigma * crandn();
		    }
		}
		for (int em_term = 0; em_term < em_terms; ++em_term) {
		    em[em_term] = 0.0;
		    if (sigma != 0.0) {
			em[em_term] += sigma * crandn();
		    }
		}
	    }
	}
	break;

    default:
	abort();
    }
}

/*
 * gen_test_terms: generate random error terms
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of calibration frequencies
 *   @stdev: standard deviation from perfect
 */
static test_terms_t *gen_test_terms(vnacal_t *vcp, vnacal_type_t type,
	int m_rows, int m_columns, int frequencies, double sigma, bool ab)
{
    test_terms_t *ttp = NULL;
    vnacal_layout_t *vlp;

    /*
     * Create the error terms structure.
     */
    if ((ttp = malloc(sizeof(test_terms_t))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	return NULL;
    }
    (void)memset((void *)ttp, 0, sizeof(*ttp));
    _vnacal_layout(&ttp->tt_layout, type, m_rows, m_columns);
    vlp = &ttp->tt_layout;
    ttp->tt_frequencies = frequencies;
    if ((ttp->tt_frequency_vector = calloc(frequencies,
		    sizeof(double))) == NULL) {
	free_test_terms(ttp);
	return NULL;
    }
    if (frequencies == 1) {
	ttp->tt_frequency_vector[0] = 1.0e+9;
    } else if (frequencies == 2) {
	ttp->tt_frequency_vector[0] = 0.0;
	ttp->tt_frequency_vector[1] = 1.0e+9;
    } else {
	ttp->tt_frequency_vector[0] = 0.0;
	for (int i = 1; i < frequencies; ++i) {
	    ttp->tt_frequency_vector[i] = pow(1.0e+9,
		(double)(i - 1) / (double)(frequencies - 2));
	}
    }
    if ((ttp->tt_error_term_vector = calloc(frequencies,
		    sizeof(double complex *))) == NULL) {
	free_test_terms(ttp);
	return NULL;
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	double complex *clfp;

	if ((clfp = calloc(VL_ERROR_TERMS(vlp),
			sizeof(double complex))) == NULL) {
	    free_test_terms(ttp);
	    return NULL;
	}
	ttp->tt_error_term_vector[findex] = clfp;
	gen_e_terms(vlp, clfp, sigma);

    }

    /*
     * Allocate the new calibration structure and set frequencies.
     */
    if ((ttp->tt_vnp = vnacal_new_alloc(vcp, type, m_rows, m_columns,
		    frequencies)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_new_alloc: %s\n",
		progname, strerror(errno));
	free_test_terms(ttp);
	return NULL;
    }
    if (vnacal_new_set_frequency_vector(ttp->tt_vnp,
		ttp->tt_frequency_vector) == -1) {
	(void)fprintf(stderr, "%s: vnacal_new_set_frequency_vector: %s\n",
		progname, strerror(errno));
	free_test_terms(ttp);
	return NULL;
    }

    /*
     * If verbose, show the error terms.
     */
    if (opt_v >= 2) {
	print_test_error_terms(ttp);
    }
    return ttp;
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
	    if (determinant == 0.0 || !isfinite(cabs(determinant))) {
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
 *   @tmp: caller-allocated test_measurements structure to hold result
 *   @s_matrix: s-parameter indices matrix describing the standard
 *   @s_matrix_rows: rows in s_matrix
 *   @s_matrix_columns: columns in s_matrix
 *   @port_map: map from standard port to VNA port
 *   @findex: frequency index
 *   @m: caller-allocated output matrix
 */
static int calc_measurements_helper(const test_terms_t *ttp,
	test_measurements_t *tmp, const int *s_matrix,
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
	    s[s_cell] = _vnacal_get_parameter_value(vpmrp, f);
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
	    s[s_cell] = crandn();
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
 * calc_measurements: calculate measurements given error terms and a standard
 *   @ttp: pointer to test error terms structure
 *   @tmp: caller-allocated test_measurements structure to hold result
 *   @s_matrix: s-parameter indices matrix describing the standard
 *   @s_matrix_rows: rows in s_matrix
 *   @s_matrix_columns: columns in s_matrix
 *   @port_map: map from standard port to VNA port
 */
static int calc_measurements(const test_terms_t *ttp, test_measurements_t *tmp,
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

	print_standard(vcp, s_matrix, s_matrix_rows, s_matrix_columns,
		ttp->tt_frequencies, ttp->tt_frequency_vector, port_map);
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
	 * If an A matrix was given, fill it with random values and
	 * replace B with B * A.
	 */
	if (a_matrix == NULL) {
	    for (int cell = 0; cell < b_rows * b_columns; ++cell) {
		b_matrix[cell][findex] = m[cell];
	    }
	} else if (VL_HAS_COLUMN_SYSTEMS(vlp)) {
	    for (int b_column = 0; b_column < b_columns; ++b_column) {
		double complex a = crandn();

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
		a[a_cell] = crandn();
		a_matrix[a_cell][findex] = a[a_cell];
	    }
	    cmatrix_multiply(b, m, a, b_rows, b_columns, b_columns);
	    for (int b_cell = 0; b_cell < b_rows * b_columns; ++b_cell) {
		b_matrix[b_cell][findex] = b[b_cell];
	    }
	}
    }

    /*
     * If verbose, show values.
     */
    if (opt_v >= 2) {
	print_test_measurements(tmp, frequencies);
    }

    return 0;
}

/*
 * gen_random_parameters: generate n random scalar parameters
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @n: number of parameters to generate
 *   @vector: caller-supplied vector to fill with parameter indices
 */
static int gen_random_parameters(vnacal_t *vcp, int n, int *vector)
{
    for (int i = 0; i < n; ++i) {
	if ((vector[i] = vnacal_make_scalar_parameter(vcp, crandn())) == -1) {
	    while (--i >= 0) {
		(void)vnacal_delete_parameter(vcp, vector[i]);
	    }
	    return -1;
	}
    }
    return 0;
}

/*
 * add_single_reflect: measure a single reflect standard on the given port
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @s11: reflection parameter
 *   @port: port to measure
 */
static int add_single_reflect(const test_terms_t *ttp,
	test_measurements_t *tmp, int s11, int port)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    const double complex **a = (const double complex **)tmp->tm_a_matrix;
    const double complex **b = (const double complex **)tmp->tm_b_matrix;

    if (calc_measurements(ttp, tmp, &s11, 1, 1, &port) == -1) {
	return -1;
    }
    if (a != NULL) {
	if (vnacal_new_add_single_reflect(vnp,
		    a, a_rows, m_columns,
		    b, m_rows, m_columns,
		    s11, port) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_single_reflect_m(vnp,
		    b, m_rows, m_columns,
		    s11, port) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * add_double_reflect: measure a double reflect standard on the given port
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @s11: first reflection parameter
 *   @s22: second reflection parameter
 *   @port1: first port
 *   @port2: second port
 */
static int add_double_reflect(const test_terms_t *ttp,
	test_measurements_t *tmp, int s11, int s22, int port1, int port2)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    const double complex **a = (const double complex **)tmp->tm_a_matrix;
    const double complex **b = (const double complex **)tmp->tm_b_matrix;
    int s_matrix[2][2] = {
	{ s11, VNACAL_ZERO },
	{ VNACAL_ZERO, s22 }
    };
    int port_map[2] = { port1, port2 };

    if (calc_measurements(ttp, tmp, &s_matrix[0][0], 2, 2, port_map) == -1) {
	return -1;
    }
    if (a != NULL) {
	if (vnacal_new_add_double_reflect(vnp,
		    a, a_rows, m_columns,
		    b, m_rows, m_columns,
		    s11, s22, port1, port2) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_double_reflect_m(vnp,
		    b, m_rows, m_columns,
		    s11, s22, port1, port2) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * add_through: measure a through standard between the given ports
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @port1: first port
 *   @port2: second port
 */
static int add_through(const test_terms_t *ttp, test_measurements_t *tmp,
	int port1, int port2)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    static int s_matrix[2][2] = { { VNACAL_MATCH, VNACAL_ONE },
				  { VNACAL_ONE, VNACAL_MATCH } };
    int port_map[2] = { port1, port2 };
    const double complex **a = (const double complex **)tmp->tm_a_matrix;
    const double complex **b = (const double complex **)tmp->tm_b_matrix;

    if (calc_measurements(ttp, tmp, &s_matrix[0][0], 2, 2, port_map) == -1) {
	return -1;
    }
    if (a != NULL) {
	if (vnacal_new_add_through(vnp,
		    a, a_rows, m_columns,
		    b, m_rows, m_columns,
		    port1, port2) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_through_m(vnp,
		    b, m_rows, m_columns,
		    port1, port2) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * add_line: measure a line standard between the given ports
 *   @ttp: pointer to test error terms structure
 *   @s_2x2: 2x2 parameter matrix
 *   @port1: first port
 *   @port2: second port
 */
static int add_line(const test_terms_t *ttp, test_measurements_t *tmp,
	const int *s_2x2, int port1, int port2)
{
    const vnacal_layout_t *vlp = &ttp->tt_layout;
    const int m_rows = VL_M_ROWS(vlp);
    const int m_columns = VL_M_COLUMNS(vlp);
    vnacal_new_t *vnp = ttp->tt_vnp;
    int a_rows = VL_HAS_COLUMN_SYSTEMS(vlp) ? 1 : m_columns;
    int port_map[2] = { port1, port2 };
    const double complex **a = (const double complex **)tmp->tm_a_matrix;
    const double complex **b = (const double complex **)tmp->tm_b_matrix;

    if (calc_measurements(ttp, tmp, s_2x2, 2, 2, port_map) == -1) {
	return -1;
    }
    if (a != NULL) {
	if (vnacal_new_add_line(vnp,
		    a, a_rows, m_columns,
		    b, m_rows, m_columns,
		    s_2x2, port1, port2) == -1) {
	    return -1;
	}
    } else {
	if (vnacal_new_add_line_m(vnp,
		    b, m_rows, m_columns,
		    s_2x2, port1, port2) == -1) {
	    return -1;
	}
    }
    return 0;
}

/*
 * validate_error_parameters: compare calculated error terms to actual
 *   @ttp: pointer to test error terms structure
 *   @calp: pointer to calibration structure
 */
static int validate_error_parameters(const test_terms_t *ttp,
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
	print_calibration(calp);
    }
    if (calp->cal_error_terms != VL_ERROR_TERMS(vlp)) {
	(void)printf("cal_error_terms (%d) != vl_error_terms (%d)\n",
		calp->cal_error_terms, VL_ERROR_TERMS(vlp));
	return -1;
    }
    for (int findex = 0; findex < ttp->tt_frequencies; ++findex) {
	for (int term = 0; term < VL_ERROR_TERMS(vlp); ++term) {
	    if (!isequal(calp->cal_error_term_vector[term][findex],
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

/*
 * run_solt_trial_helper: add short, open and load calibrations on port
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @port: port to calibrate
 */
static test_result_type run_solt_trial_helper(test_terms_t *ttp,
	test_measurements_t *tmp, int port)
{
    if (add_single_reflect(ttp, tmp, VNACAL_SHORT, port) == -1) {
	return T_FAIL;
    }
    if (add_single_reflect(ttp, tmp, VNACAL_OPEN, port) == -1) {
	return T_FAIL;
    }
    if (add_single_reflect(ttp, tmp, VNACAL_MATCH, port) == -1) {
	return T_FAIL;
    }
    return T_PASS;
}

/*
 * run_vnacal_new_solt_trial: test 8-12 parameter SOLT calibration
 *   @trial: test trial
 *   @type: error term type
 *   @rows: number of VNA ports that detect signal
 *   @columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 *   @ab: true: use a, b matrices; false: use m matrix
 */
static test_result_type run_vnacal_new_solt_trial(int trial,
	vnacal_type_t type, int rows, int columns, int frequencies, bool ab)
{
    vnacal_t *vcp = NULL;
    test_terms_t *ttp = NULL;
    test_measurements_t *tmp = NULL;
    int diagonals = MIN(rows, columns);
    int ports = MAX(rows, columns);
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_new: trial %3d size %d x %d "
		"type %-4s %s SOLT\n",
	    trial, rows, columns, _vnacal_type_to_name(type), ab ? "AB" : "M ");
    }

    /*
     * Create the calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate random error parameters.
     */
    if ((ttp = gen_test_terms(vcp, type, rows, columns,
		    frequencies, 1.0, ab)) == NULL) {
	(void)fprintf(stderr, "%s: gen_test_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate the test measurement matrices.
     */
    if ((tmp = alloc_test_measurements(type, rows, columns,
		    frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * For E12 and UE14, we have to do short, open and load calibration
     * on every diagonal port.  For the others, we can choose any one
     * diagonal port.
     */
    if (type == VNACAL_E12 || VNACAL_IS_UE14(type)) {
	for (int port = 1; port <= diagonals; ++port) {
	    if ((result = run_solt_trial_helper(ttp, tmp, port)) != T_PASS) {
		goto out;
	    }
	}
    } else {
	int port = random() % diagonals + 1;

	if ((result = run_solt_trial_helper(ttp, tmp, port)) != T_PASS) {
	    goto out;
	}
    }

    /*
     * Do through tests between every diagonal port and every other port.
     */
    for (int port1 = 1; port1 <= diagonals; ++port1) {
	for (int port2 = port1 + 1; port2 <= ports; ++port2) {
	    if (add_through(ttp, tmp, port1, port2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	}
    }

    /*
     * Solve for the error parameters and check.
     */
    if (vnacal_new_solve(ttp->tt_vnp) == -1) {
	(void)fprintf(stderr, "%s: vnacal_solve: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (validate_error_parameters(ttp, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    free_test_measurements(tmp);
    free_test_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_solt: run SOLT tests for 8-12 term parameters
 */
static void test_vnacal_new_solt()
{
    static const int sizes[] = { 1, 2, 3, 4 };
    static const vnacal_type_t types[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10, VNACAL_UE14, VNACAL_E12
    };
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
        for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
            for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
                int rows = sizes[si];
                int columns = sizes[sj];

		for (int ti = 0; ti < sizeof(types) / sizeof(types[0]); ++ti) {
		    vnacal_type_t type = types[ti];

		    if (type == VNACAL_T8 || type == VNACAL_TE10) {
			if (rows > columns) {
			    continue;
			}
		    } else {
			if (rows < columns) {
			    continue;
			}
		    }
		    result = run_vnacal_new_solt_trial(trial, type,
			    rows, columns, 2, false);
		    if (result != T_PASS)
			goto out;
		    result = run_vnacal_new_solt_trial(trial, type,
			    rows, columns, 2, true);
		    if (result != T_PASS)
			goto out;
		}
            }
        }
    }
    result = T_PASS;

out:
    report_test_result("vnacal_new SOLT", result);
}

#define MM	 0
#define MO	 1
#define MS	 2
#define OM	 3
#define OO	 4
#define OS	 5
#define SM	 6
#define SO	 7
#define SS	 8
#define T	 9
#define D	10

static const char *const standard_names[] = {
    /*  0 */ "MM",
    /*  1 */ "MO",
    /*  2 */ "MS",
    /*  3 */ "OM",
    /*  4 */ "OO",
    /*  5 */ "OS",
    /*  6 */ "SM",
    /*  7 */ "SO",
    /*  8 */ "SS",
    /*  9 */ "T",
    /* 10 */ "D"
};

/*
 * silvonen_table: 16-Term Silvonen Calibration Sequences
 *
 * From: Silvonen, Kimmo. (1994). New five-standard calibration procedures
 * for network analyzers and wafer probes. NASA STI/Recon Technical Report N.
 */
static const int silvonen_table[][6] = {
    { T, MM, SS, OO, SO, -1 },
    { T, MM, SS, OO, OS, -1 },
    { T, MM, SS, OO, SM, -1 },
    { T, MM, SS, OO, MS, -1 },
    { T, MM, SS, OO, OM, -1 },
    { T, MM, SS, OO, MO, -1 },
    { T, MM, SS, SO, OS, -1 },
    { T, MM, SS, SO, MS, -1 },
    { T, MM, SS, SO, MO, -1 },
    { T, MM, SS, OS, SM, -1 },
    { T, MM, SS, OS, OM, -1 },
    { T, MM, SS, OM, MO, -1 },
    { T, MM, SS, SM, MS, -1 },
    { T, MM, SS, SM, MO, -1 },
    { T, MM, SS, MS, OM, -1 },
    { T, MM, OS, SM, MS, -1 },
    { T, MM, OS, SM, MO, -1 },
    { T, MM, OS, MS, OM, -1 },
    { T, MM, OO, SO, OS, -1 },
    { T, MM, OO, SO, SM, -1 },
    { T, MM, OO, SO, OM, -1 },
    { T, MM, OO, OS, MS, -1 },
    { T, MM, OO, OS, MO, -1 },
    { T, MM, OO, SM, MS, -1 },
    { T, MM, OO, SM, MO, -1 },
    { T, MM, OO, MS, OM, -1 },
    { T, MM, OO, OM, MO, -1 },
    { T, MM, SO, OS, SM, -1 },
    { T, MM, SO, OS, MS, -1 },
    { T, MM, SO, OS, OM, -1 },
    { T, MM, SO, OS, MO, -1 },
    { T, MM, SO, SM, MS, -1 },
    { T, MM, SO, SM, MO, -1 },
    { T, MM, SO, MS, OM, -1 },
    { T, MM, SO, OM, MO, -1 },
    { T, MM, OS, OM, MO, -1 },
    { T, SS, OO, SO, MS, -1 },
    { T, SS, OO, SO, OM, -1 },
    { T, SS, OO, OS, SM, -1 },
    { T, SS, OO, OS, MO, -1 },
    { T, SS, OO, SM, MS, -1 },
    { T, SS, OO, SM, OM, -1 },
    { T, SS, OO, MS, MO, -1 },
    { T, SS, OO, OM, MO, -1 },
    { T, SS, SO, OS, OM, -1 },
    { T, SS, SO, OS, MO, -1 },
    { T, SS, SO, MS, OM, -1 },
    { T, SS, SO, MS, MO, -1 },
    { T, SS, SO, OM, MO, -1 },
    { T, SS, OS, SM, OM, -1 },
    { T, SS, OS, SM, MO, -1 },
    { T, SS, OS, OM, MO, -1 },
    { T, SS, SM, MS, OM, -1 },
    { T, SS, SM, MS, MO, -1 },
    { T, SS, SM, OM, MO, -1 },
    { T, SS, MS, OM, MO, -1 },
    { T, SO, OS, MS, OM, -1 },
    { T, SO, OS, MS, MO, -1 },
    { T, OO, SO, OS, SM, -1 },
    { T, OO, SO, OS, MS, -1 },
    { T, OO, SO, SM, MS, -1 },
    { T, OO, SO, SM, OM, -1 },
    { T, OO, SO, MS, OM, -1 },
    { T, OO, OS, SM, MS, -1 },
    { T, OO, OS, SM, MO, -1 },
    { T, OO, OS, MS, MO, -1 },
    { T, OO, SM, MS, OM, -1 },
    { T, OO, SM, MS, MO, -1 },
    { T, OO, SM, OM, MO, -1 },
    { T, OO, MS, OM, MO, -1 },
    { T, SO, OS, SM, OM, -1 },
    { T, SO, OS, SM, MO, -1 },
    { T, SO, SM, MS, OM, -1 },
    { T, SO, SM, MS, MO, -1 },
    { T, SO, SM, OM, MO, -1 },
    { T, SO, MS, OM, MO, -1 },
    { T, OS, SM, MS, OM, -1 },
    { T, OS, SM, MS, MO, -1 },
    { T, OS, SM, OM, MO, -1 },
    { T, OS, MS, OM, MO, -1 },
    { T, D,  MM, SS, SO, -1 },
    { T, D,  MM, SS, OS, -1 },
    { T, D,  MM, SS, SM, -1 },
    { T, D,  MM, SS, MS, -1 },
    { T, D,  MM, SS, OM, -1 },
    { T, D,  MM, SS, MO, -1 },
    { T, D,  MM, OO, SO, -1 },
    { T, D,  MM, OO, OS, -1 },
    { T, D,  MM, OO, SM, -1 },
    { T, D,  MM, OO, MS, -1 },
    { T, D,  MM, OO, OM, -1 },
    { T, D,  MM, OO, MO, -1 },
    { T, D,  MM, SO, SM, -1 },
    { T, D,  MM, SO, MS, -1 },
    { T, D,  MM, SO, OM, -1 },
    { T, D,  MM, SO, MO, -1 },
    { T, D,  MM, OS, SM, -1 },
    { T, D,  MM, OS, MS, -1 },
    { T, D,  MM, OS, OM, -1 },
    { T, D,  MM, OS, MO, -1 },
    { T, D,  MM, SM, MS, -1 },
    { T, D,  MM, SM, MO, -1 },
    { T, D,  MM, MS, OM, -1 },
    { T, D,  MM, OM, MO, -1 },
    { T, D,  OO, SO, OS, -1 },
    { T, D,  OO, SO, SM, -1 },
    { T, D,  OO, SO, MS, -1 },
    { T, D,  OO, SO, OM, -1 },
    { T, D,  OO, OS, SM, -1 },
    { T, D,  OO, OS, MS, -1 },
    { T, D,  OO, OS, MO, -1 },
    { T, D,  OO, SM, MS, -1 },
    { T, D,  OO, SM, OM, -1 },
    { T, D,  OO, SM, MO, -1 },
    { T, D,  OO, MS, OM, -1 },
    { T, D,  OO, MS, MO, -1 },
    { T, D,  OO, OM, MO, -1 },
    { T, D,  SO, OS, SM, -1 },
    { T, D,  SO, OS, MS, -1 },
    { T, D,  SO, OS, OM, -1 },
    { T, D,  SO, OS, MO, -1 },
    { T, D,  SO, SM, MS, -1 },
    { T, D,  SO, SM, OM, -1 },
    { T, D,  SO, SM, MO, -1 },
    { T, D,  SO, MS, OM, -1 },
    { T, D,  SO, MS, MO, -1 },
    { T, D,  SO, OM, MO, -1 },
    { T, D,  OS, SM, MS, -1 },
    { T, D,  OS, SM, OM, -1 },
    { T, D,  OS, SM, MO, -1 },
    { T, D,  OS, MS, OM, -1 },
    { T, D,  OS, MS, MO, -1 },
    { T, D,  OS, OM, MO, -1 },
    { T, D,  SM, MS, OM, -1 },
    { T, D,  SM, MS, MO, -1 },
    { T, D,  SM, OM, MO, -1 },
    { T, D,  MS, OM, MO, -1 },
    { T, D,  SS, OO, SO, -1 },
    { T, D,  SS, OO, OS, -1 },
    { T, D,  SS, OO, SM, -1 },
    { T, D,  SS, OO, MS, -1 },
    { T, D,  SS, OO, OM, -1 },
    { T, D,  SS, OO, MO, -1 },
    { T, D,  SS, SO, OS, -1 },
    { T, D,  SS, SO, MS, -1 },
    { T, D,  SS, SO, OM, -1 },
    { T, D,  SS, SO, MO, -1 },
    { T, D,  SS, OS, SM, -1 },
    { T, D,  SS, OS, OM, -1 },
    { T, D,  SS, OS, MO, -1 },
    { T, D,  SS, SM, MS, -1 },
    { T, D,  SS, SM, OM, -1 },
    { T, D,  SS, SM, MO, -1 },
    { T, D,  SS, MS, OM, -1 },
    { T, D,  SS, MS, MO, -1 },
    { T, D,  SS, OM, MO, -1 },
};

/*
 * test_vnacal_new_table_entry: add calibration standards from a table row
 *   @trial: test trial
 *   @type: error term type
 *   @rows: number of VNA ports that detect signal
 *   @columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 *   @ab: true: use a, b matrices; false: use m matrix
 */
static test_result_type test_vnacal_new_table_entry(int trial,
	vnacal_type_t type, int frequencies, const int *table_entry, bool ab)
{
    vnacal_t *vcp = NULL;
    test_terms_t *ttp = NULL;
    test_measurements_t *tmp = NULL;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_new: trial %3d size 2 x 2 type %s %s:",
	    trial, _vnacal_type_to_name(type), ab ? "AB" : "M ");
	for (const int *ip = table_entry; *ip != -1; ++ip) {
	    (void)printf(" %s", standard_names[*ip]);
	}
	(void)printf("\n");
    }

    /*
     * Create the calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate random error parameters.
     */
    if ((ttp = gen_test_terms(vcp, type, 2, 2,
		    frequencies, 1.0, ab)) == NULL) {
	(void)fprintf(stderr, "%s: gen_test_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate the test measurement matrices.
     */
    if ((tmp = alloc_test_measurements(type, 2, 2, frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add standards based on the table.
     */
    for (const int *ip = table_entry; *ip != -1; ++ip) {
	switch (*ip) {
	case MM:
	    if (add_double_reflect(ttp, tmp, VNACAL_MATCH, VNACAL_MATCH,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case MO:
	    if (add_double_reflect(ttp, tmp, VNACAL_MATCH, VNACAL_OPEN,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case MS:
	    if (add_double_reflect(ttp, tmp, VNACAL_MATCH, VNACAL_SHORT,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case OM:
	    if (add_double_reflect(ttp, tmp, VNACAL_OPEN, VNACAL_MATCH,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case OO:
	    if (add_double_reflect(ttp, tmp, VNACAL_OPEN, VNACAL_OPEN,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case OS:
	    if (add_double_reflect(ttp, tmp, VNACAL_OPEN, VNACAL_SHORT,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case SM:
	    if (add_double_reflect(ttp, tmp, VNACAL_SHORT, VNACAL_MATCH,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case SO:
	    if (add_double_reflect(ttp, tmp, VNACAL_SHORT, VNACAL_OPEN,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case SS:
	    if (add_double_reflect(ttp, tmp, VNACAL_SHORT, VNACAL_SHORT,
			1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case T:
	    if (add_through(ttp, tmp, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case D:
	    {
		const double complex g_through = -I * sqrt(2.0) / 2.0;
		const double complex g_reflect = sqrt(2.0) / 2.0;
		int p1, p2;
		int s[2][2];

		if ((p1 = vnacal_make_scalar_parameter(vcp, g_through)) == -1) {
		    result = T_FAIL;
		    goto out;
		}
		if ((p2 = vnacal_make_scalar_parameter(vcp, g_reflect)) == -1) {
		    vnacal_delete_parameter(vcp, p1);
		    result = T_FAIL;
		    goto out;
		}
		s[0][0] = p2;
		s[0][1] = p1;
		s[1][0] = p1;
		s[1][1] = p2;
		if (add_line(ttp, tmp, &s[0][0], 1, 2) == -1) {
		    vnacal_delete_parameter(vcp, p2);
		    vnacal_delete_parameter(vcp, p1);
		    result = T_FAIL;
		    goto out;
		}
		vnacal_delete_parameter(vcp, p2);
		vnacal_delete_parameter(vcp, p1);
	    }
	    break;

	default:
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Solve for the error parameters and check.
     */
    if (vnacal_new_solve(ttp->tt_vnp) == -1) {
	(void)fprintf(stderr, "%s: vnacal_solve: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (validate_error_parameters(ttp, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    free_test_measurements(tmp);
    free_test_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_solt: run 16 parameter 2 port tests from Silvonen table
 */
static void test_vnacal_new_silvonen16()
{
    static const vnacal_type_t types[] = {
	VNACAL_T16, VNACAL_U16
    };
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= 10; ++trial) {
	for (int entry = 0; entry < sizeof(silvonen_table) /
		sizeof(silvonen_table[0]); ++entry) {
	    for (int ti = 0; ti < sizeof(types) / sizeof(types[0]); ++ti) {
		vnacal_type_t type = types[ti];

		result = test_vnacal_new_table_entry(trial, type, 2,
			silvonen_table[entry], false);
		if (result != T_PASS)
		    goto out;
		result = test_vnacal_new_table_entry(trial, type, 2,
			silvonen_table[entry], true);
		if (result != T_PASS)
		    goto out;
	    }
	}
    }
    result = T_PASS;

out:
    report_test_result("vnacal_new Silvonen 16", result);
}

#define DIVROUND(k, n)	(((k) + (n) - 1) / (n))

/*
 * calc_n_needed_random_standards: calculate the number standards needed
 *   @type: error term type
 *   @m_rows: rows in the measurement matrix
 *   @m_columns: columns in the measurement matrix
 *   @add_all_match: set to true if all-match standard needed
 *
 * This function may sometimes overestimate T8, U8, T16 and U16 where
 * we add an extra standard.
 */
static int calc_n_needed_random_standards(vnacal_type_t type,
	int m_rows, int m_columns, bool *add_all_match)
{
    const int ports = MAX(m_rows, m_columns);
    int standards = 0;

    *add_all_match = false;
    if (ports == 1) {			/* special-case */
	return 3;
    }
    switch (type) {
    case VNACAL_T8:
    case VNACAL_U8:
	{
	    int terms = 2 * (m_rows + m_columns) - 1;

	    standards = DIVROUND(terms, m_rows * m_columns) + 1;
	}
	break;

    case VNACAL_TE10:
    case VNACAL_UE10:
	{
	    int terms = 2 * (m_rows + m_columns) - 1;

	    *add_all_match = true;
	    standards = DIVROUND(terms, m_rows * m_columns);
	}
	break;

    case VNACAL_T16:
    case VNACAL_U16:
	{
	    int terms = (m_rows + m_columns) * 2 * ports - 1;

	    standards = DIVROUND(terms, m_rows * m_columns) + 1;
	}
	break;

    case VNACAL_UE14:
    case _VNACAL_E12_UE14:
    case VNACAL_E12:
	{
	    int terms = m_columns * (2 * m_rows + 1);

	    *add_all_match = true;
	    standards = DIVROUND(terms, m_rows * m_columns);
	}
	break;

    default:
	abort();
    }
    return standards;
}

/*
 * make_random_calibration: make a random calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 *   @ab: true: use a, b matrices; false: use m matrix
 */
static test_terms_t *make_random_calibration(vnacal_t *vcp, vnacal_type_t type,
	int m_rows, int m_columns, int frequencies, bool ab)
{
    test_terms_t *ttp = NULL;
    const int ports = MAX(m_rows, m_columns);
    bool add_all_match;
    int standards;
    test_measurements_t *tmp = NULL;
    const double complex **a;
    const double complex **b;

    /*
     * Generate random error parameters.
     */
    if ((ttp = gen_test_terms(vcp, type, m_rows, m_columns,
		    frequencies, 1.0, false)) == NULL) {
	(void)fprintf(stderr, "%s: gen_test_terms: %s\n",
		progname, strerror(errno));
	goto error;
    }

    /*
     * Calculate the number of standards needed.
     */
    standards = calc_n_needed_random_standards(type, m_rows, m_columns,
	    &add_all_match);

    /*
     * Allocate the measurements matrices.
     */
    if ((tmp = alloc_test_measurements(type, m_rows, m_columns,
		    frequencies, ab)) == NULL) {
	goto error;
    }
    a = (const double complex **)tmp->tm_a_matrix;
    b = (const double complex **)tmp->tm_b_matrix;

    /*
     * If needed, add an all match matrix.
     */
    if (add_all_match) {
	int s[ports * ports];

	for (int i = 0; i < ports * ports; ++i) {
	    s[i] = VNACAL_MATCH;
	}
	if (calc_measurements(ttp, tmp, s, ports, ports, NULL) == -1) {
	    goto error;
	}
	if (ab) {
	    if (vnacal_new_add_mapped_matrix(ttp->tt_vnp,
			a, tmp->tm_a_rows, tmp->tm_a_columns,
			b, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	} else {
	    if (vnacal_new_add_mapped_matrix_m(ttp->tt_vnp,
			b, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	}
    }

    /*
     * Add random standards.
     */
    for (int standard = 0; standard < standards; ++standard) {
	int s[ports * ports];

	if (gen_random_parameters(vcp, ports * ports, s) == -1) {
	    goto error;
	}
	if (calc_measurements(ttp, tmp, s, ports, ports, NULL) == -1) {
	    goto error;
	}
	if (ab) {
	    if (vnacal_new_add_mapped_matrix(ttp->tt_vnp,
			a, tmp->tm_a_rows, tmp->tm_a_columns,
			b, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	} else {
	    if (vnacal_new_add_mapped_matrix_m(ttp->tt_vnp,
			b, m_rows, m_columns,
			s, ports, ports, NULL) == -1) {

		goto error;
	    }
	}
	for (int i = 0; i < ports * ports; ++i) {
	    if (vnacal_delete_parameter(vcp, s[i]) == -1) {
		goto error;
	    }
	}
    }
    free_test_measurements(tmp);
    tmp = NULL;

    /*
     * Solve for the error parameters and check.
     */
    if (vnacal_new_solve(ttp->tt_vnp) == -1) {
	(void)fprintf(stderr, "%s: vnacal_solve: %s\n",
		progname, strerror(errno));
	goto error;
    }
    if (validate_error_parameters(ttp, NULL) == -1) {
	goto error;
    }
    return ttp;

error:
    free_test_measurements(tmp);
    free_test_terms(ttp);
    return NULL;
}

/*
 * run_vnacal_new_random_trial: calibrate with random multi-port standards
 *   @trial: test trial
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 *   @ab: true: use a, b matrices; false: use m matrix
 */
static test_result_type run_vnacal_new_random_trial(int trial,
	vnacal_type_t type, int m_rows, int m_columns,
	int frequencies, bool ab)
{
    vnacal_t *vcp = NULL;
    test_terms_t *ttp = NULL;
    bool add_all_match;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	int standards;

	standards = calc_n_needed_random_standards(type, m_rows, m_columns,
		&add_all_match);
	(void)printf("Test vnacal_new: trial %3d size %d x %d "
		"type %-4s %s %2d random standards%s\n",
		trial, m_rows, m_columns, _vnacal_type_to_name(type),
		ab ? "AB" : "M ", standards, add_all_match ? "+match" : "");
    }

    /*
     * Create the calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Make the calibration, solve and check.
     */
    if ((ttp = make_random_calibration(vcp, type, m_rows, m_columns,
		    frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    free_test_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_random: test vnacal_new_* with random multi-port standards
 */
static void test_vnacal_new_random()
{
    static const int sizes[] = { 1, 2, 3, 4 };
    static const vnacal_type_t types[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10, VNACAL_T16, VNACAL_U16,
	VNACAL_UE14, VNACAL_E12
    };
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= 12; ++trial) {
        for (int si = 0; si < sizeof(sizes) / sizeof(int); ++si) {
            for (int sj = 0; sj < sizeof(sizes) / sizeof(int); ++sj) {
                int rows = sizes[si];
                int columns = sizes[sj];

		for (int ti = 0; ti < sizeof(types) / sizeof(types[0]); ++ti) {
		    vnacal_type_t type = types[ti];

		    if (type == VNACAL_T8 || type == VNACAL_TE10 ||
			    type == VNACAL_T16) {
			if (rows > columns) {
			    continue;
			}
		    } else {
			if (rows < columns) {
			    continue;
			}
		    }
		    result = run_vnacal_new_random_trial(trial, type,
			    rows, columns, 2, false);
		    if (result != T_PASS)
			goto out;
		    result = run_vnacal_new_random_trial(trial, type,
			    rows, columns, 2, true);
		    if (result != T_PASS)
			goto out;
		}
            }
        }
    }
    result = T_PASS;

out:
    report_test_result("vnacal_new random standards", result);
}

/*
 * run_vnacal_apply_trial: test vnacal_apply
 *   @trial: test trial
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 */
static test_result_type run_vnacal_apply_trial(int trial,
	vnacal_type_t type, int m_rows, int m_columns,
	int frequencies, bool ab)
{
    const int ports = MAX(m_rows, m_columns);
    vnacal_t *vcp = NULL;
    test_terms_t *ttp = NULL;
    test_measurements_t *tmp = NULL;
    int ci = -1;
    int s[ports * ports];
    const double complex **a;
    const double complex **b;
    vnadata_t *vdp = NULL;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_apply: trial %3d size %d x %d type %-4s %s\n",
	    trial, m_rows, m_columns, _vnacal_type_to_name(type),
	    ab ? "AB" : "M ");
    }

    /*
     * Create the calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Make the requested calibration.
     */
    if ((ttp = make_random_calibration(vcp, type, m_rows, m_columns,
		    frequencies, false)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate a test measurement structure to hold the DUT measurements.
     */
    if ((tmp = alloc_test_measurements(type, ports, ports,
		    frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }
    a = (const double complex **)tmp->tm_a_matrix;
    b = (const double complex **)tmp->tm_b_matrix;

    /*
     * Add it to the vnacal_t structure.
     */
    if ((ci = vnacal_add_calibration(vcp, "cal1", ttp->tt_vnp)) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Create random s-parameters for the DUT.
     */
    if (gen_random_parameters(vcp, ports * ports, s) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (calc_measurements(ttp, tmp, s, ports, ports, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Create a vnadata_t structure to hold the result.
     */
    if ((vdp = vnadata_alloc()) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Apply the calibration.
     */
    if (ab) {
	if (vnacal_apply(vcp, ci, ttp->tt_frequency_vector, ttp->tt_frequencies,
		    a, tmp->tm_a_rows, tmp->tm_b_columns,
		    b, ports, ports, vdp) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    } else {
	if (vnacal_apply_m(vcp, ci, ttp->tt_frequency_vector,
		    ttp->tt_frequencies, b, ports, ports, vdp) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Check the result.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = ttp->tt_frequency_vector[findex];

	if (opt_v >= 2) {
	    (void)printf("findex %d  f %e\n", findex, f);
	    (void)printf("  expected s parameters:\n");
	    for (int s_row = 0; s_row < ports; ++s_row) {
		(void)printf("  ");
		for (int s_column = 0; s_column < ports; ++s_column) {
		    int s_cell = s_row * ports + s_column;
		    vnacal_parameter_t *vpmrp;
		    double complex v;

		    if ((vpmrp = _vnacal_get_parameter(vcp,
				    s[s_cell])) == NULL) {
			result = T_FAIL;
			goto out;
		    }
		    v = _vnacal_get_parameter_value(vpmrp, f);

		    (void)printf(" %8.5f%+8.5fj", creal(v), cimag(v));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	    (void)printf("  computed s parameters:\n");
	    for (int s_row = 0; s_row < ports; ++s_row) {
		(void)printf("  ");
		for (int s_column = 0; s_column < ports; ++s_column) {
		    double complex v;

		    v = vnadata_get_cell(vdp, findex, s_row, s_column);
		    (void)printf(" %8.5f%+8.5fj", creal(v), cimag(v));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	}
	for (int s_row = 0; s_row < ports; ++s_row) {
	    for (int s_column = 0; s_column < ports; ++s_column) {
		int s_cell = s_row * ports + s_column;
		vnacal_parameter_t *vpmrp;
		double complex expected, actual;

		if ((vpmrp = _vnacal_get_parameter(vcp, s[s_cell])) == NULL) {
		    result = T_FAIL;
		    goto out;
		}
		expected = _vnacal_get_parameter_value(vpmrp, f);
		actual = vnadata_get_cell(vdp, findex, s_row, s_column);
		if (!isequal(actual, expected)) {
		    if (opt_a) {
			assert(!"data miscompare");
		    }
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    for (int i = 0; i < ports * ports; ++i) {
	(void)vnacal_delete_parameter(vcp, s[i]);
    }
    result = T_PASS;

out:
    vnadata_free(vdp);
    free_test_measurements(tmp);
    free_test_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_solt: test vnacal_new_* with random multi-port standards
 */
static void test_vnacal_apply()
{
    static const int sizes[][2] = {
	{  1,  1 },
	{  1,  2 },
	{  2,  1 },
	{  2,  2 },
	{  3,  3 },
	{  4,  4 },
	{  5,  5 }
    };
    static const vnacal_type_t types[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10, VNACAL_T16, VNACAL_U16,
	VNACAL_UE14, VNACAL_E12
    };
    test_result_type result = T_SKIPPED;

    for (int trial = 1; trial <= 12; ++trial) {
	for (int s = 0; s < sizeof(sizes) / sizeof(sizes[0]); ++s) {
	    int rows = sizes[s][0];
	    int columns = sizes[s][1];

	    for (int ti = 0; ti < sizeof(types) / sizeof(types[0]); ++ti) {
		vnacal_type_t type = types[ti];

		if (type == VNACAL_T8 || type == VNACAL_TE10 ||
			type == VNACAL_T16) {
		    if (rows > columns) {
			continue;
		    }
		} else {
		    if (rows < columns) {
			continue;
		    }
		}
		result = run_vnacal_apply_trial(trial, type,
			rows, columns, 2, false);
		if (result != T_PASS)
		    goto out;
		result = run_vnacal_apply_trial(trial, type,
			rows, columns, 2, true);
		if (result != T_PASS)
		    goto out;
	    }
        }
    }
    result = T_PASS;

out:
    report_test_result("vnacal_apply", result);
}

#if 0
/************************************************************
 * Code needed for vnacal_map_* functions (yet to be written)
 * that replace the older (broken) vnacal_apply_t functions.
 ************************************************************/

/*
 * map_type: calibration and DUT matrix dimensions and port maps
 */
typedef struct apply_test_case {
    int atc_vrows;		/* calibration matrix rows */
    int atc_vcolumns;		/* calibration matrix columns */
    int atc_drows;		/* DUT matrix rows */
    int atc_dcolumns;		/* DUT matrix columns */
    const int *const *atc_maps;	/* optional vector of maps */
} apply_test_case_type;

/*
 * apply_test_cases: VNA dimensions, DUT dimensions and port maps to test
 *	The idea here is to determine the S parameters of a DUT with
 *	more ports than the VNA.  Unused DUT ports must be terminated.
 */
static const apply_test_case_type apply_test_cases[] = {
    { 1, 1, 1, 1, NULL },
    { 1, 2, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1 },
	    NULL
	}
    },
    { 1, 2, 1, 2, NULL },
    { 1, 2, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    NULL
	}
    },
    { 1, 2, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    NULL
	}
    },
    { 1, 2, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    NULL
	}
    },
    { 1, 2, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    NULL
	}
    },
    { 1, 2, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    NULL
	}
    },
    { 1, 2, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    NULL
	}
    },
    { 1, 2, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    NULL
	}
    },
    { 1, 2, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    NULL
	}
    },
    { 1, 2, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    NULL
	}
    },
    { 1, 2, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 2, 3 },
	    NULL
	}
    },
    { 1, 2, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    (const int []){ 3, 0 },
	    NULL
	}
    },
    { 1, 2, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 3, 0 },
	    (const int []){ 3, 1 },
	    NULL
	}
    },
    { 1, 2, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 3, 0 },
	    (const int []){ 3, 1 },
	    (const int []){ 3, 2 },
	    NULL
	}
    },
    { 1, 2, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 2, 3 },
	    (const int []){ 3, 0 },
	    (const int []){ 3, 1 },
	    (const int []){ 3, 2 },
	    NULL
	}
    },
    { 1, 3, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1 },
	    NULL
	}
    },
    { 1, 3, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 1, 3, 1, 3, NULL },
    { 1, 3, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 1, 3, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    (const int []){ 1, 0, -1 },
	    NULL
	}
    },
    { 1, 3, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    (const int []){ 1, 0, -1 },
	    NULL
	}
    },
    { 1, 3, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    NULL
	}
    },
    { 1, 3, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 1, 0, 3 },
	    NULL
	}
    },
    { 1, 3, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    NULL
	}
    },
    { 1, 3, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    NULL
	}
    },
    { 1, 3, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    NULL
	}
    },
    { 1, 3, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 1, 0, 3 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 2, 0, 3 },
	    NULL
	}
    },
    { 1, 3, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 3, 0, 1 },
	    NULL
	}
    },
    { 1, 3, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 3, 0, 1 },
	    NULL
	}
    },
    { 1, 3, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 3, 0, 1 },
	    (const int []){ 3, 0, 2 },
	    NULL
	}
    },
    { 1, 3, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 1, 0, 3 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 2, 0, 3 },
	    (const int []){ 3, 0, 1 },
	    (const int []){ 3, 0, 2 },
	    NULL
	}
    },
    { 1, 4, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 1, 4, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 1, 4, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 1, 4, 1, 4, NULL },
    { 1, 4, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    (const int []){ 1, 0, -1, -1 },
	    NULL
	}
    },
    { 1, 4, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    (const int []){ 1, 0, -1, -1 },
	    NULL
	}
    },
    { 1, 4, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    NULL
	}
    },
    { 1, 4, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    NULL
	}
    },
    { 1, 4, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    (const int []){ 2, 0, 1, -1 },
	    NULL
	}
    },
    { 1, 4, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    (const int []){ 2, 0, 1, -1 },
	    NULL
	}
    },
    { 1, 4, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    (const int []){ 2, 0, 1, -1 },
	    NULL
	}
    },
    { 1, 4, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    NULL
	}
    },
    { 1, 4, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 1, 4, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 1, 4, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 1, 4, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 2, 1, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1 },
	    NULL
	}
    },
    { 2, 1, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    NULL
	}
    },
    { 2, 1, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    NULL
	}
    },
    { 2, 1, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    (const int []){ 3, 0 },
	    NULL
	}
    },
    { 2, 1, 2, 1, NULL },
    { 2, 1, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    NULL
	}
    },
    { 2, 1, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    NULL
	}
    },
    { 2, 1, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 1, 0 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 3, 0 },
	    (const int []){ 3, 1 },
	    NULL
	}
    },
    { 2, 1, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    NULL
	}
    },
    { 2, 1, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    NULL
	}
    },
    { 2, 1, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    NULL
	}
    },
    { 2, 1, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 3, 0 },
	    (const int []){ 3, 1 },
	    (const int []){ 3, 2 },
	    NULL
	}
    },
    { 2, 1, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    NULL
	}
    },
    { 2, 1, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    NULL
	}
    },
    { 2, 1, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 2, 3 },
	    NULL
	}
    },
    { 2, 1, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 0 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 0 },
	    (const int []){ 2, 1 },
	    (const int []){ 2, 3 },
	    (const int []){ 3, 0 },
	    (const int []){ 3, 1 },
	    (const int []){ 3, 2 },
	    NULL
	}
    },
    { 2, 2, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1 },
	    NULL
	}
    },
    { 2, 2, 1, 2, NULL },
    { 2, 2, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    NULL
	}
    },
    { 2, 2, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    NULL
	}
    },
    { 2, 2, 2, 1, NULL },
    { 2, 2, 2, 2, NULL },
    { 2, 2, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 2 },
	    NULL
	}
    },
    { 2, 2, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    NULL
	}
    },
    { 2, 2, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    NULL
	}
    },
    { 2, 2, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 2 },
	    NULL
	}
    },
    { 2, 2, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 1, 2 },
	    NULL
	}
    },
    { 2, 2, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 3 },
	    NULL
	}
    },
    { 2, 2, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    NULL
	}
    },
    { 2, 2, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    NULL
	}
    },
    { 2, 2, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 3 },
	    NULL
	}
    },
    { 2, 2, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1 },
	    (const int []){ 0, 2 },
	    (const int []){ 0, 3 },
	    (const int []){ 1, 2 },
	    (const int []){ 1, 3 },
	    (const int []){ 2, 3 },
	    NULL
	}
    },
    { 2, 3, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1 },
	    NULL
	}
    },
    { 2, 3, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 2, 3, 1, 3, NULL },
    { 2, 3, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 2, 3, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 2, 3, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 2, 3, 2, 3, NULL },
    { 2, 3, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 2, 3, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 1 },
	    NULL
	}
    },
    { 2, 3, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 1 },
	    NULL
	}
    },
    { 2, 3, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 1 },
	    NULL
	}
    },
    { 2, 3, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 3 },
	    (const int []){ 1, 2, 3 },
	    NULL
	}
    },
    { 2, 3, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 2, 3, 0 },
	    NULL
	}
    },
    { 2, 3, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 2, 1 },
	    (const int []){ 1, 3, 0 },
	    NULL
	}
    },
    { 2, 3, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 2, 3, 0 },
	    (const int []){ 2, 3, 1 },
	    NULL
	}
    },
    { 2, 3, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 2, 3, 0 },
	    (const int []){ 2, 3, 1 },
	    NULL
	}
    },
    { 2, 4, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 2, 4, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 2, 4, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 2, 4, 1, 4, NULL },
    { 2, 4, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 2, 4, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 2, 4, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 2, 4, 2, 4, NULL },
    { 2, 4, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 0, 2, 1, -1 },
	    NULL
	}
    },
    { 2, 4, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 0, 2, 1, -1 },
	    NULL
	}
    },
    { 2, 4, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 0, 2, 1, -1 },
	    NULL
	}
    },
    { 2, 4, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 2, 1, 3 },
	    NULL
	}
    },
    { 2, 4, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 2, 4, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 2, 4, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 2, 4, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 3, 1, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1 },
	    NULL
	}
    },
    { 3, 1, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    (const int []){ 1, 0, -1 },
	    NULL
	}
    },
    { 3, 1, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    NULL
	}
    },
    { 3, 1, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 3, 0, 1 },
	    NULL
	}
    },
    { 3, 1, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 1, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    (const int []){ 1, 0, -1 },
	    NULL
	}
    },
    { 3, 1, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    NULL
	}
    },
    { 3, 1, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 3, 0, 1 },
	    NULL
	}
    },
    { 3, 1, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    NULL
	}
    },
    { 3, 1, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    NULL
	}
    },
    { 3, 1, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    NULL
	}
    },
    { 3, 1, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 3, 0, 1 },
	    (const int []){ 3, 0, 2 },
	    NULL
	}
    },
    { 3, 1, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 1, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 1, 0, 3 },
	    NULL
	}
    },
    { 3, 1, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 1, 0, 3 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 2, 0, 3 },
	    NULL
	}
    },
   { 3, 1, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 1, 0, 2 },
	    (const int []){ 1, 0, 3 },
	    (const int []){ 2, 0, 1 },
	    (const int []){ 2, 0, 3 },
	    (const int []){ 3, 0, 1 },
	    (const int []){ 3, 0, 2 },
	    NULL
	}
    },
    { 3, 2, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1 },
	    NULL
	}
    },
    { 3, 2, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 2, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 1 },
	    NULL
	}
    },
    { 3, 2, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 2, 3, 0 },
	    NULL
	}
    },
    { 3, 2, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 2, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 2, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 1 },
	    NULL
	}
    },
    { 3, 2, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 2, 1 },
	    (const int []){ 1, 3, 0 },
	    NULL
	}
    },
    { 3, 2, 3, 1, NULL },
    { 3, 2, 3, 2, NULL },
    { 3, 2, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 1 },
	    NULL
	}
    },
    { 3, 2, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 2, 3, 0 },
	    (const int []){ 2, 3, 1 },
	    NULL
	}
    },
    { 3, 2, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 2, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 2, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 2, 3 },
	    (const int []){ 1, 2, 3 },
	    NULL
	}
    },
    { 3, 2, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 2, 3, 0 },
	    (const int []){ 2, 3, 1 },
	    NULL
	}
    },
    { 3, 3, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1 },
	    NULL
	}
    },
    { 3, 3, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 3, 1, 3, NULL },
    { 3, 3, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 3, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 3, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1 },
	    NULL
	}
    },
    { 3, 3, 2, 3, NULL },
    { 3, 3, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 3, 3, 1, NULL },
    { 3, 3, 3, 2, NULL },
    { 3, 3, 3, 3, NULL },
    { 3, 3, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 0, 2, 3 },
	    NULL
	}
    },
    { 3, 3, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 3, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    NULL
	}
    },
    { 3, 3, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 0, 2, 3 },
	    NULL
	}
    },
    { 3, 3, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2 },
	    (const int []){ 0, 1, 3 },
	    (const int []){ 0, 2, 3 },
	    NULL
	}
    },
    { 3, 4, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 3, 4, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 3, 4, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 3, 4, 1, 4, NULL },
    { 3, 4, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 3, 4, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 3, 4, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 3, 4, 2, 4, NULL },
    { 3, 4, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 3, 4, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 3, 4, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 3, 4, 3, 4, NULL },
    { 3, 4, 4, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 3, 4, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 3, 4, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 3, 4, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 4, 1, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 4, 1, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    (const int []){ 1, 0, -1, -1 },
	    NULL
	}
    },
    { 4, 1, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    (const int []){ 2, 0, 1, -1 },
	    NULL
	}
    },
    { 4, 1, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 4, 1, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 1, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    (const int []){ 1, 0, -1, -1 },
	    NULL
	}
    },
    { 4, 1, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    (const int []){ 2, 0, 1, -1 },
	    NULL
	}
    },
    { 4, 1, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 4, 1, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 1, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    NULL
	}
    },
    { 4, 1, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 1, 0, 2, -1 },
	    (const int []){ 2, 0, 1, -1 },
	    NULL
	}
    },
    { 4, 1, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 4, 1, 4, 1, NULL },
    { 4, 1, 4, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    NULL
	}
    },
    { 4, 1, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    NULL
	}
    },
    { 4, 1, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 1, 0, 2, 3 },
	    (const int []){ 2, 0, 1, 3 },
	    (const int []){ 3, 0, 1, 2 },
	    NULL
	}
    },
    { 4, 2, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 4, 2, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 2, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 0, 2, 1, -1 },
	    NULL
	}
    },
    { 4, 2, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 4, 2, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 2, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 2, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 0, 2, 1, -1 },
	    NULL
	}
    },
    { 4, 2, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 4, 2, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 2, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 2, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    (const int []){ 0, 2, 1, -1 },
	    NULL
	}
    },
    { 4, 2, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 4, 2, 4, 1, NULL },
    { 4, 2, 4, 2, NULL },
    { 4, 2, 4, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 2, 1, 3 },
	    NULL
	}
    },
    { 4, 2, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 2, 3, 0, 1 },
	    NULL
	}
    },
    { 4, 3, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 4, 3, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 3, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 3, 1, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 4, 3, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 3, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 3, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 3, 2, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 4, 3, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 3, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 3, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 3, 3, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 4, 3, 4, 1, NULL },
    { 4, 3, 4, 2, NULL },
    { 4, 3, 4, 3, NULL },
    { 4, 3, 4, 4,
	(const int *const[]){
	    (const int []){ 0, 1, 2, 3 },
	    (const int []){ 0, 1, 3, 2 },
	    NULL
	}
    },
    { 4, 4, 1, 1,
	(const int *const[]){
	    (const int []){ 0, -1, -1, -1 },
	    NULL
	}
    },
    { 4, 4, 1, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 4, 1, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 4, 1, 4, NULL },
    { 4, 4, 2, 1,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 4, 2, 2,
	(const int *const[]){
	    (const int []){ 0, 1, -1, -1 },
	    NULL
	}
    },
    { 4, 4, 2, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 4, 2, 4, NULL },
    { 4, 4, 3, 1,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 4, 3, 2,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 4, 3, 3,
	(const int *const[]){
	    (const int []){ 0, 1, 2, -1 },
	    NULL
	}
    },
    { 4, 4, 3, 4, NULL },
    { 4, 4, 4, 1, NULL },
    { 4, 4, 4, 2, NULL },
    { 4, 4, 4, 3, NULL },
    { 4, 4, 4, 4, NULL },
};
#define N_APPLY_CASES	\
	(sizeof(apply_test_cases) / sizeof(apply_test_case_type))

/*
 * test_vnacal_apply_helper
 *   @trial: trial number
 *   @vrows: rows in calibration matrix
 *   @vcolumns: columns in calibration matrix
 *   @drows: rows in DUT matrix
 *   @dcolumns: columns in DUT matrix
 *   @frequencies: number of frequency points
 *   @map_flag: map DUT ports to VNA ports (required if DUT larger than cal)
 */
static test_result_type test_vnacal_apply_helper(int trial, int frequencies,
	const apply_test_case_type *atcp)
{
    int vrows    = atcp->atc_vrows;
    int vcolumns = atcp->atc_vcolumns;
    int drows    = atcp->atc_drows;
    int dcolumns = atcp->atc_dcolumns;
    int vports = MAX(vrows, vcolumns);
    vnacal_calibration_t *cmsp = NULL;
    vnacal_apply_t *vap = NULL;
    double complex *(*error_terms)[3] = NULL;
    vnacal_t *vcp = NULL;
    double complex **actual_matrix   = NULL;
    double complex **measured_matrix = NULL;
    const int *const *map;
    vnadata_t *output_matrix = NULL;
    test_result_type result = T_SKIPPED;
#define S(i, j)	\
	(actual_matrix[(i) * dcolumns + (j)])	/* drows x dcolumns */
#define E(i, j, t) \
	(error_terms[(i) * vcolumns + (j)][t])	/* vrows x vcolumns */
#define M(i, j)	\
	(measured_matrix[(i) * vcolumns + (j)])	/* vrows x vcolumns */

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_oapply: trial %3d cal size (%d x %d) "
		"S size (%d x %d) map %d\n", trial,
		vrows, vcolumns, drows, dcolumns,
		atcp->atc_maps != NULL);
    }

    /*
     * Create calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate the error terms and calibration measurements.
     */
    if ((cmsp = vnacal_calset_alloc(vcp, "test", VNACAL_UE14,
		    vrows, vcolumns, frequencies)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_calset_alloc: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((error_terms = gen_test_terms(cmsp)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Create a new vnacal_t based on the calibration measurements.
     */
    if (vnacal_add_calibration(vcp, cmsp) == -1) {
	(void)fprintf(stderr, "%s: vnacal_add_calibration: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (opt_v != 0) {
	vnacal_oetermset_t *etsp = vcp->vc_set_vector[0];

	(void)printf("error terms:\n");
	(void)printf("R C F ET\n");
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int ci = 0; ci < vrows; ++ci) {
		for (int cj = 0; cj < vcolumns; ++cj) {
		    int c_cell = ci * vcolumns + cj;
		    double complex **epp = error_terms[c_cell];
		    vnacal_error_terms_t *etp;

		    etp = &etsp->ets_error_term_matrix[c_cell];
		    for (int k = 0; k < 3; ++k) {
			(void)printf("%d %d %d %-6s %+e%+ei %+e%+ei\n",
				ci, cj, findex, error_term_names[k],
				creal(etp->et_data_vectors[k][findex]),
				cimag(etp->et_data_vectors[k][findex]),
				creal(epp[k][findex]),
				cimag(epp[k][findex]));
		    }
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * Generate the "actual" S-parameters.
     */
    actual_matrix = alloc_matrix_of_vectors(drows * dcolumns, frequencies);
    if (actual_matrix == NULL) {
	result = T_FAIL;
	goto out;
    }
    for (int row = 0; row < drows; ++row) {
	for (int column = 0; column < dcolumns; ++column) {
	    for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
		actual_matrix[row * dcolumns + column][findex] = crandn();
	    }
	}
    }
    if (opt_v != 0) {
	(void)printf("actual_matrix:\n");
	(void)printf("R C F\n");
	for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
	    for (int row = 0; row < drows; ++row) {
		for (int column = 0; column < dcolumns; ++column) {
		    double complex v;

		    v = actual_matrix[row * dcolumns + column][findex];
		    (void)printf("%d %d %d %+e%+ei\n",
			row, column, findex, creal(v), cimag(v));
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * Create the vnacal_apply_t.
     */
    vap = vnacal_apply_alloc(vcp, /*calibration=*/0,
	    drows, dcolumns, frequencies);
    if (vap == NULL) {
	(void)fprintf(stderr, "%s: vnacal_apply_alloc: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_apply_set_frequency_vector(vap,
		calp->cal_frequency_vector) == -1) {
	(void)fprintf(stderr, "%s: vnacal_apply_set_frequency_vector: "
		"%s\n", progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate the vrows x vcolumns measured matrix.
     */
    measured_matrix = alloc_matrix_of_vectors(vrows * vcolumns, frequencies);
    if (measured_matrix == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * For each port map...
     */
    map = atcp->atc_maps;
    do {
	/*
	 * If opt_v, show the map.
	 */
	if (opt_v != 0 && map != NULL) {
	    (void)printf("map:\n");
	    for (int i = 0; i < vports; ++i) {
		(void)printf(" %d", (*map)[i]);
	    }
	    (void)printf("\n\n");
	}

	/*
	 * For each frequency...
	 */
	for (int findex = 0; findex < calp->cal_frequencies; ++findex) {

	    if (opt_v != 0) {
		(void)printf("findex %d:\n", findex);
	    }

	    /*
	     * Init measured to zero.
	     */
	    for (int i = 0; i < vrows * vcolumns; ++i) {
		measured_matrix[i][findex] = 0.0;
	    }

	    /*
	     * For each vcolumn (each driven port) find the corresponding
	     * column in the measured_matrix.
	     */
	    for (int vcolumn = 0; vcolumn < vcolumns; ++vcolumn) {
		double complex a[vrows * vrows];
		double complex x[vrows];
		double complex b[vrows];
		double complex d;
		int dcolumn = (map != NULL) ? (*map)[vcolumn] : vcolumn;
#define A(i, j) (a[(i) * vrows + (j)])

		/*
		 * We start by forming a (drows x drows) matrix A and column
		 * vector b which we'll use to solve for column vector x.
		 * From x, we can easily calculate the VNA measurements for
		 * the current port mapping.
		 *
		 * Initialize A to the identity matrix and b to the
		 * mapped row and column of the actual S matrix.
		 */
		for (int i = 0; i < vrows; ++i) {
		    for (int j = 0; j < vrows; ++j) {
			A(i, j) = (i == j) ? 1.0 : 0.0;
		    }
		}
		for (int vrow = 0; vrow < vrows; ++vrow) {
		    int drow = (map != NULL) ? (*map)[vrow] : vrow;

		    if (drow >= 0 && drow < drows &&
			    dcolumn >= 0 && dcolumn < dcolumns) {
			b[vrow] = S(drow, dcolumn)[findex];
		    } else {
			b[vrow] = 0.0;
		    }
		}
		/*
		 * Make A = I - S E, where E is a diagonal matrix made of
		 * the e11/e22 error terms for this column.
		 */
		for (int i = 0; i < vrows; ++i) {	/* row in A */
		    int ii = (map != NULL) ? (*map)[i] : i;

		    if (ii < 0 || ii >= drows) {
			continue;
		    }
		    for (int j = 0; j < vrows; ++j) {	/* each vrow */
			int jj = (map != NULL) ? (*map)[j] : j;

			if (jj < 0 || jj >= dcolumns) {
			    continue;
			}
			A(i, j) -= S(ii, jj)[findex] * E(j, vcolumn, 2)[findex];
		    }
		}
		if (opt_v != 0) {
		    (void)printf("vcolumn %d dcolumn %d:\n", vcolumn, dcolumn);
		    (void)printf("a:\n");
		    cmatrix_print(a, vrows, vrows);
		    (void)printf("b:\n");
		    cmatrix_print(b, vrows, 1);
		}
		/*
		 * Find x = A^-1 b.
		 */
		d = _vnacommon_mldivide(x, a, b, vrows, 1);
		if (cabs(d) <= EPS)	 {
		    (void)fprintf(stderr, "%s: test_vnacal_mrdivide: warning: "
			    "skipping nearly singular test matrix\n",
			    progname);
		    result = T_SKIPPED;
		    goto out;
		}
		if (opt_v != 0) {
		    (void)printf("x:\n");
		    cmatrix_print(x, vrows, 1);
		}
		/*
		 * From x, calculate the "measured" S-parameters for this
		 * column.
		 */
		for (int vrow = 0; vrow < vrows; ++vrow) {
		    double complex e00    = E(vrow, vcolumn, 0)[findex];
		    double complex e10e01 = E(vrow, vcolumn, 1)[findex];

		    M(vrow, vcolumn)[findex] = e00 + e10e01 * x[vrow];
		}
#undef A
	    }
	}
	if (opt_v != 0) {
	    (void)printf("measured_matrix:\n");
	    (void)printf("R C F\n");
	    for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
		for (int vrow = 0; vrow < vrows; ++vrow) {
		    for (int vcolumn = 0; vcolumn < vcolumns; ++vcolumn) {
			double complex v = M(vrow, vcolumn)[findex];

			(void)printf("%d %d %d %+e%+ei\n",
				vrow, vcolumn, findex, creal(v), cimag(v));
		    }
		}
	    }
	    (void)printf("\n");
	}
	if (vnacal_apply_add_matrix(vap,
		    (const double complex *const *)&M(0, 0),
		    map != NULL ? *map : NULL) == -1) {
	    (void)fprintf(stderr, "%s: vnacal_apply_add_matrix: "
		    "%s\n", progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	/*
	 * If there's no port map and the DUT matrix has the same
	 * dimension as the calibration matrix, test using the simple
	 * vnacal_oapply function.
	 */
	if (map == NULL && vrows == drows && vcolumns == dcolumns) {
	    if ((output_matrix = vnadata_alloc()) == NULL)  {
		(void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    if (vnacal_oapply(vcp, /*ci*/0, frequencies,
			calp->cal_frequency_vector,
			(const double complex *const *)&M(0, 0),
			output_matrix) == -1) {
		(void)fprintf(stderr, "%s: vnacal_oapply: "
			"%s\n", progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    if (opt_v != 0) {
		(void)printf("computed_vector (vnacal_oapply):\n");
		(void)printf("R C F\n");
		for (int findex = 0; findex < calp->cal_frequencies;
			++findex) {
		    for (int row = 0; row < drows; ++row) {
			for (int column = 0; column < dcolumns; ++column) {
			    double complex v;

			    v = vnadata_get_cell(output_matrix, findex,
				    row, column);
			    (void)printf("%d %d %d %+e%+ei\n",
				    row, column, findex,
				    creal(v), cimag(v));
			}
		    }
		}
		(void)printf("\n");
	    }
	    for (int i = 0; i < drows; ++i) {
		for (int j = 0; j < dcolumns; ++j) {
		    for (int findex = 0; findex < calp->cal_frequencies;
			    ++findex) {
			double complex v;
			double dy;

			v = vnadata_get_cell(output_matrix, findex, i, j);
			dy = cabs(v - actual_matrix[i * dcolumns + j][findex]);
			if (dy >= EPS) {
			    if (opt_a) {
				assert(!"data miscompare");
			    }
			    result = T_FAIL;
			    goto out;
			}
		    }
		}
	    }
	    (void)vnadata_free(output_matrix);
	    output_matrix = NULL;
	}
    } while (map != NULL && *++map != NULL);

    /*
     * Get the computed S-parameters.
     */
    if ((output_matrix = vnadata_alloc()) == NULL)  {
	(void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_apply_get_data(vap, output_matrix) == -1) {
	(void)fprintf(stderr, "%s: vnacal_apply_get_data: "
		"%s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (opt_v != 0) {
	(void)printf("computed_vector (vnacal_apply_get_data):\n");
	(void)printf("R C F\n");
	for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
	    for (int row = 0; row < drows; ++row) {
		for (int column = 0; column < dcolumns; ++column) {
		    double complex v;

		    v = vnadata_get_cell(output_matrix, findex, row, column);
		    (void)printf("%d %d %d %+e%+ei\n",
			row, column, findex, creal(v), cimag(v));
		}
	    }
	}
	(void)printf("\n");
    }

    /*
     * Check the result.
     */
    for (int i = 0; i < drows; ++i) {
	for (int j = 0; j < dcolumns; ++j) {
	    for (int findex = 0; findex < calp->cal_frequencies; ++findex) {
		double complex v;
		double dy;

		v = vnadata_get_cell(output_matrix, findex, i, j);
		dy = cabs(v - actual_matrix[i * dcolumns + j][findex]);
		if (dy >= EPS) {
		    if (opt_a) {
			assert(!"data miscompare");
		    }
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    result = T_PASS;

out:
    if (vap != NULL) {
	vnacal_apply_free(vap);
	vap = NULL;
    }
    free_matrix_of_vectors(measured_matrix, vrows * vcolumns);
    free_matrix_of_vectors(actual_matrix,   drows * dcolumns);
    if (cmsp != NULL) {
	if (error_terms != NULL) {
	    free_test_terms(error_terms, cmsp);
	    error_terms = NULL;
	}
	vnacal_calset_free(cmsp);
	cmsp = NULL;
    }
    if (vcp != NULL) {
	vnacal_free(vcp);
	vcp = NULL;
    }
    return result;
}
#undef M
#undef E
#undef S

/*
 * test_vnacal_map_apply: test vnacal_map
 */
static void test_vnacal_map_apply()
{
    test_result_type result = T_SKIPPED;
    bool pass = false;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
	for (const apply_test_case_type *atcp = apply_test_cases;
		atcp < &apply_test_cases[N_APPLY_CASES]; ++atcp) {
	    result = test_vnacal_apply_helper(trial,2, atcp);
	    switch (result) {
		case T_PASS:
		    pass = true;
		    break;
		case T_SKIPPED:
		    break;
		case T_FAIL:
		default:
		    goto out;
	    }
	}
    }
    result = pass ? T_PASS : T_SKIPPED;

out:
    report_test_result("vnacal_oapply", result);
}
#endif /* end of mapped vnacal_apply tests */

/*
 * Test Strings for vnacal_property_set
 */
static const char property_foo_value[] = "1234567890";
static const char property_bar_value[] = "abcdefghijkl\nmnopqrstuvwxyz";
static const char property3_value[] = "αβγδεζηθικλμνξοπρστυφχψω";

/*
 * run_vnacal_save_load_trial
 */
static test_result_type run_vnacal_save_load_trial(int trial)
{
    static const int dimension_table[][2] = {
	{ 1, 1 }, { 1, 2 }, { 1, 3 }, { 1, 4 }, { 2, 2 },
	{ 2, 3 }, { 2, 4 }, { 3, 3 }, { 3, 4 }, { 4, 4 },
    };
    static const vnacal_type_t type_table[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10,
	VNACAL_T16, VNACAL_U16, VNACAL_UE14, VNACAL_E12
    };
    test_terms_t *ttp_table[8] =
        { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
    vnacal_t *vcp = NULL;
    const int types = sizeof(type_table) / sizeof(vnacal_type_t);
    const int dimensions = sizeof(dimension_table) / sizeof(int [2]);
    test_result_type result = T_SKIPPED;
    const char *cp_temp;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_save, vnacal_load: trial %d\n", trial);
    }

    /*
     * Create calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Generate error terms and save to the vnacal_t structure.
     */
    for (int tindex = 0; tindex < types; ++tindex) {
	vnacal_type_t type = type_table[tindex];
	int frequencies = random() % 3 + 1;
	int dindex = random() % dimensions;
	int m_rows = 0, m_columns = 0;

	switch (type) {
	case VNACAL_T8:
	case VNACAL_TE10:
	case VNACAL_T16:
	    m_rows = dimension_table[dindex][0];
	    m_columns = dimension_table[dindex][1];
	    break;

	case VNACAL_U8:
	case VNACAL_UE10:
	case VNACAL_UE14:
	case _VNACAL_E12_UE14:
	case VNACAL_U16:
	case VNACAL_E12:
	    m_rows = dimension_table[dindex][1];
	    m_columns = dimension_table[dindex][0];
	    break;
	}
	if ((ttp_table[tindex] = make_random_calibration(vcp, type,
			m_rows, m_columns, frequencies, /*ab*/false)) == NULL) {
	    result = T_FAIL;
	    goto out;
	}
	if (vnacal_add_calibration(vcp, _vnacal_type_to_name(type),
		    ttp_table[tindex]->tt_vnp) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	vnacal_new_free(ttp_table[tindex]->tt_vnp);
	ttp_table[tindex]->tt_vnp = NULL;
    }

    /*
     * Set test properties.
     */
    if (vnacal_property_set(vcp, -1, "global_property=47") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 0, "foo=999999999999") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 0, "bar=%s", property_bar_value) == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 0, "foo=%s", property_foo_value) == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 1, "baz=!!!") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_set(vcp, 1, "property3=%s", property3_value) == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnacal_property_delete(vcp, 1, "baz") == -1) {
	(void)fprintf(stderr, "%s: vnacal_property_delete: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    {
	int rows = ttp_table[0]->tt_layout.vl_m_rows;
	int columns = ttp_table[0]->tt_layout.vl_m_columns;

	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		int cell = row * columns + column;
		int value = (cell + 1) % (rows * columns);

		if (vnacal_property_set(vcp, 0, "switches[%d][%d]=%d",
			    row, column, value) == -1) {
		    (void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
			    progname, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    {
	int rows = ttp_table[1]->tt_layout.vl_m_rows;
	int columns = ttp_table[1]->tt_layout.vl_m_columns;

	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		int cell = row * columns + column;
		int value = (cell + 3) % (rows * columns);

		if (vnacal_property_set(vcp, 1, "switches[%d][%d]=%d",
			    row, column, value) == -1) {
		    (void)fprintf(stderr, "%s: vnacal_property_set: %s\n",
			    progname, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }

    /*
     * Save
     */
    if (vnacal_save(vcp, "vnacal-test.vnacal") == -1) {
	(void)fprintf(stderr, "%s: vnacal_save: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    vnacal_free(vcp);
    vcp = NULL;

    /*
     * Load
     */
    if ((vcp = vnacal_load("vnacal-test.vnacal", error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_load: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Validate error parameters.
     */
    for (int tindex = 0; tindex < types; ++tindex) {
	vnacal_type_t type = type_table[tindex];
	vnacal_calibration_t *calp = NULL;
	int ci = -1;

	if ((ci = vnacal_find_calibration(vcp,
			_vnacal_type_to_name(type))) == -1) {
	    (void)fprintf(stderr, "%s: vnacal_find_calibration: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	    (void)fprintf(stderr, "%s: _vnacal_get_calibration: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (validate_error_parameters(ttp_table[tindex], calp) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Validate properties.
     */
    if ((cp_temp = vnacal_property_get(vcp, -1, "global_property")) == NULL) {
	(void)printf("property \"global_property\" not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, "47") != 0) {
	(void)printf("expected \"47\" for property \"global_property\"; "
		"found \"%s\"\n", cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 0, "foo")) == NULL) {
	(void)printf("property \"foo\" in calibration 0 not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, property_foo_value) != 0) {
	(void)printf("expected \"%s\" for property \"foo\"; found \"%s\"\n",
		property_foo_value, cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 0, "bar")) == NULL) {
	(void)printf("property \"bar\" in calibration 0 not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, property_bar_value) != 0) {
	(void)printf("expected \"%s\" for property \"bar\"; found \"%s\"\n",
		property_bar_value, cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 0, "baz")) != NULL) {
	(void)printf("property \"baz\" not expected in calibration 0; "
		"found it with value \"%s\"\n", cp_temp);
	result = T_FAIL;
	goto out;
    }
    if ((cp_temp = vnacal_property_get(vcp, 1, "property3")) == NULL) {
	(void)printf("property \"property3\" in calibration 1 not found\n");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(cp_temp, property3_value) != 0) {
	(void)printf("expected \"%s\" for property \"property3\"; "
		"found \"%s\"\n", property3_value, cp_temp);
	result = T_FAIL;
	goto out;
    }
    {
	int rows = ttp_table[1]->tt_layout.vl_m_rows;
	int columns = ttp_table[1]->tt_layout.vl_m_columns;

	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		int cell = row * columns + column;
		int value = (cell + 3) % (rows * columns);

		if ((cp_temp = vnacal_property_get(vcp, 1, "switches[%d][%d]",
			    row, column)) == NULL) {
		    (void)fprintf(stderr, "%s: vnacal_property_get: "
			    "switches[%d][%d] in calibration 0 not found: %s\n",
			    progname, row, column, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
		if (atoi(cp_temp) != value) {
		    (void)fprintf(stderr, "%s: ci 1: expected %d for property "
			    "switches[%d][%d]; found \"%s\"\n",
			    progname, value, row, column, cp_temp);
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    {
	int rows = ttp_table[0]->tt_layout.vl_m_rows;
	int columns = ttp_table[0]->tt_layout.vl_m_columns;

	for (int row = 0; row < rows; ++row) {
	    for (int column = 0; column < columns; ++column) {
		int cell = row * columns + column;
		int value = (cell + 1) % (rows * columns);

		if ((cp_temp = vnacal_property_get(vcp, 0, "switches[%d][%d]",
			    row, column)) == NULL) {
		    (void)fprintf(stderr, "%s: vnacal_property_get: "
			    "switches[%d][%d] in calibration 0 not found: %s\n",
			    progname, row, column, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
		if (atoi(cp_temp) != value) {
		    (void)fprintf(stderr, "%s: ci 0: expected %d for property "
			    "switches[%d][%d]; found \"%s\"\n",
			    progname, value, row, column, cp_temp);
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    result = T_PASS;

out:
    for (int tindex = 0; tindex < types; ++tindex) {
	free_test_terms(ttp_table[tindex]);
    }
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_save_load
 */
static void test_vnacal_save_load()
{
    test_result_type result = T_SKIPPED;

    for (int trial = 0; trial < 5; ++trial) {
	result = run_vnacal_save_load_trial(trial);
	if (result != T_PASS)
	    goto out;
    }
    result = T_PASS;

out:
    report_test_result("vnacal_save/vnacal_load", result);
}

/*
 * CV2_F: number of freuquencies for the compat_V2 test
 */
#define CV2_F	11

/*
 * compat_V2_frequency_vector: frequency vector for compat_V2 test
 */
static const double compat_V2_frequency_vector[CV2_F] = {
    1.000000e+05,
    1.584893e+05,
    2.511886e+05,
    3.981072e+05,
    6.309573e+05,
    1.000000e+06,
    1.584893e+06,
    2.511886e+06,
    3.981072e+06,
    6.309573e+06,
    1.000000e+07
};

/*
 * compat_V2_measured: "measured" s-parameters for old "VNACAL 2.0" format
 *
 * These tables were generated using the E12 example with 11 calibration
 * and 11 measurement points from 100 kHz to 10 MHz.
 */
static const double complex compat_V2_measured[4][CV2_F] = {
    {	/* s11 */
	-3.926540e-03 + 4.341532e-04*I,
	-9.773344e-03 + 1.726204e-03*I,
	-2.397616e-02 + 6.848069e-03*I,
	-5.649961e-02 + 2.696819e-02*I,
	-1.171773e-01 + 1.031452e-01*I,
	-1.379310e-01 + 3.448276e-01*I,
	+2.440045e-01 + 6.724525e-01*I,
	+8.548239e-01 + 3.846625e-01*I,
	+9.586034e-01 - 2.523617e-01*I,
	+6.399219e-01 - 7.672778e-01*I,
	+1.061320e-01 - 9.942893e-01*I
    },
    {	/* s12 */
	+9.939136e-01 - 1.099960e-01*I,
	+9.845757e-01 - 1.742980e-01*I,
	+9.604276e-01 - 2.759055e-01*I,
	+8.958771e-01 - 4.339283e-01*I,
	+7.175210e-01 - 6.559979e-01*I,
	+2.891602e-01 - 8.052561e-01*I,
	-1.570320e-01 - 5.267873e-01*I,
	-1.809236e-01 - 1.774419e-01*I,
	-9.240888e-02 - 4.711767e-02*I,
	-3.972649e-02 - 1.194356e-02*I,
	-1.625128e-02 - 3.004438e-03*I
    },
    {	/* s21 */
	+9.939350e-01 - 1.098983e-01*I,
	+9.847092e-01 - 1.739230e-01*I,
	+9.612490e-01 - 2.745518e-01*I,
	+9.006954e-01 - 4.299166e-01*I,
	+7.414183e-01 - 6.526327e-01*I,
	+3.448276e-01 - 8.620690e-01*I,
	-2.383455e-01 - 6.568568e-01*I,
	-3.176208e-01 - 1.429263e-01*I,
	-1.275371e-01 + 3.357539e-02*I,
	-2.705835e-02 + 3.244345e-02*I,
	-1.185832e-03 + 1.110938e-02*I
    },
    {	/* s22 */
	+6.013177e-03 - 6.654756e-04*I,
	+1.496251e-02 - 2.648791e-03*I,
	+3.666232e-02 - 1.053212e-02*I,
	+8.590211e-02 - 4.160766e-02*I,
	+1.728184e-01 - 1.580003e-01*I,
	+1.749419e-01 - 4.871800e-01*I,
	-2.386401e-01 - 8.005541e-01*I,
	-6.906384e-01 - 6.773475e-01*I,
	-8.860722e-01 - 4.517927e-01*I,
	-9.568318e-01 - 2.876665e-01*I,
	-9.832025e-01 - 1.817685e-01*I
    }
};

/*
 * compat_V2_m: measurement matrix for vnacal_apply_m
 */
static const double complex *compat_V2_m[4] = {
    compat_V2_measured[0], compat_V2_measured[1],
    compat_V2_measured[2], compat_V2_measured[3],
};

/*
 * compat_V2_expected: expected s-parameters for old "VNACAL 2.0" format
 */
static const double complex compat_V2_expected[4][CV2_F] = {
    {	/* s11 */
	-4.974876e-03 + 4.999875e-04*I,
	-1.239974e-02 + 1.990222e-03*I,
	-3.052222e-02 + 7.916587e-03*I,
	-7.250960e-02 + 3.135099e-02*I,
	-1.533550e-01 + 1.208076e-01*I,
	-2.000000e-01 + 4.000000e-01*I,
	+1.247191e-01 + 7.723058e-01*I,
	+6.206602e-01 + 7.235185e-01*I,
	+8.601119e-01 + 4.945027e-01*I,
	+9.473713e-01 + 3.161807e-01*I,
	+9.796082e-01 + 1.999200e-01*I
    },
    {	/* s12 */
	+9.949751e-01 - 9.999750e-02*I,
	+9.872848e-01 - 1.584643e-01*I,
	+9.674892e-01 - 2.509389e-01*I,
	+9.150093e-01 - 3.956228e-01*I,
	+7.704206e-01 - 6.069102e-01*I,
	+4.000000e-01 - 8.000000e-01*I,
	-9.930313e-02 - 6.149210e-01*I,
	-1.967360e-01 - 2.293399e-01*I,
	-1.085388e-01 - 6.240202e-02*I,
	-4.759378e-02 - 1.588420e-02*I,
	-1.959216e-02 - 3.998401e-03*I
    },
    {	/* s21 */

	+9.949751e-01 - 9.999750e-02*I,
	+9.872848e-01 - 1.584643e-01*I,
	+9.674892e-01 - 2.509389e-01*I,
	+9.150093e-01 - 3.956228e-01*I,
	+7.704206e-01 - 6.069102e-01*I,
	+4.000000e-01 - 8.000000e-01*I,
	-9.930313e-02 - 6.149210e-01*I,
	-1.967360e-01 - 2.293399e-01*I,
	-1.085388e-01 - 6.240202e-02*I,
	-4.759378e-02 - 1.588420e-02*I,
	-1.959216e-02 - 3.998401e-03*I
    },
    {	/* s22 */
	+4.974876e-03 - 4.999875e-04*I,
	+1.239974e-02 - 1.990222e-03*I,
	+3.052222e-02 - 7.916587e-03*I,
	+7.250960e-02 - 3.135099e-02*I,
	+1.533550e-01 - 1.208076e-01*I,
	+2.000000e-01 - 4.000000e-01*I,
	-1.247191e-01 - 7.723058e-01*I,
	-6.206602e-01 - 7.235185e-01*I,
	-8.601119e-01 - 4.945027e-01*I,
	-9.473713e-01 - 3.161807e-01*I,
	-9.796082e-01 - 1.999200e-01*I
    }
};

/*
 * test_vnacal_load_compat_e2: test compatibility load of old E term format
 */
static void test_vnacal_load_compat_V2()
{
    vnacal_t *vcp = NULL;
    vnadata_t *vdp = NULL;
    test_result_type result = T_SKIPPED;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_load VNACAL 2.0 format\n");
    }

    /*
     * Load the old format save file.
     */
    if ((vcp = vnacal_load("compat-V2.vnacal", error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_load: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Create a vnadata_t structure to hold the result.
     */
    if ((vdp = vnadata_alloc()) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Apply the calibration.
     */
    if (vnacal_apply_m(vcp, 0, compat_V2_frequency_vector, CV2_F,
		compat_V2_m, 2, 2, vdp) == -1) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Check the result.
     */
    for (int findex = 0; findex < CV2_F; ++findex) {
	double f = compat_V2_frequency_vector[findex];

	if (opt_v >= 2) {
	    (void)printf("findex %d  f %e\n", findex, f);
	    (void)printf("  computed s parameters:\n");
	    for (int s_row = 0; s_row < 2; ++s_row) {
		(void)printf("  ");
		for (int s_column = 0; s_column < 2; ++s_column) {
		    double complex v;

		    v = vnadata_get_cell(vdp, findex, s_row, s_column);
		    (void)printf(" %8.5f%+8.5fj", creal(v), cimag(v));
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	}
	for (int s_row = 0; s_row < 2; ++s_row) {
	    for (int s_column = 0; s_column < 2; ++s_column) {
		int s_cell = s_row * 2 + s_column;
		double complex expected, actual;

		expected = compat_V2_expected[s_cell][findex];
		actual = vnadata_get_cell(vdp, findex, s_row, s_column);
		if (!isequal(actual, expected)) {
		    if (opt_a) {
			assert(!"data miscompare");
		    }
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    }
    result = T_PASS;

out:
    vnadata_free(vdp);
    vnacal_free(vcp);
    report_test_result("vnacal_load VNACAL 2.0 format", result);
}

/*
 * print_usage: print a usage message and exit
 */
static void print_usage()
{
    const char *const *cpp;

    for (cpp = usage; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s: usage %s\n", progname, *cpp);
    }
    for (cpp = help; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s\n", *cpp);
    }
    exit(2);
    /*NOTREACHED*/
}

/*
 * main
 */
int
main(int argc, char **argv)
{
    /*
     * Parse Options
     */
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case -1:
	    break;

	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	default:
	    print_usage();
	    /*NOTREACHED*/
	}
	break;
    }
    test_vnacal_new_solt();
    test_vnacal_new_silvonen16();
    test_vnacal_new_random();
    test_vnacal_apply();
    test_vnacal_save_load();
    test_vnacal_load_compat_V2();

    exit(fail_count != 0);
    /*NOTREACHED*/
}
