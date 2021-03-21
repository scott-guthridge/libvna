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

#include "src/archdep.h"

#include <assert.h>
#include <complex.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "src/vnacal_internal.h"
#include "test.h"
#include "vnacaltest.h"


#define NTRIALS		50

/*
 * Command Line Options
 */
char *progname;
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
bool opt_a = false;
int  opt_v = 0;

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
 * run_vnacal_apply_trial: test vnacal_apply
 *   @trial: test trial
 *   @type: error term type
 *   @m_rows: number of VNA ports that detect signal
 *   @m_columns: number of VNA ports that generate signal
 *   @frequencies: number of test frequenciens
 */
static test_result_t run_vnacal_apply_trial(int trial,
	vnacal_type_t type, int m_rows, int m_columns,
	int frequencies, bool ab)
{
    const int ports = MAX(m_rows, m_columns);
    vnacal_t *vcp = NULL;
    test_vnacal_terms_t *ttp = NULL;
    test_vnacal_measurements_t *tmp = NULL;
    int ci = -1;
    int s[ports * ports];
    vnadata_t *vdp = NULL;
    test_result_t result = T_FAIL;

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
    if ((tmp = test_vnacal_alloc_measurements(type, ports, ports,
		    frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }

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
    if (test_vnacal_generate_random_parameters(vcp, s, ports * ports) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (test_vnacal_calculate_measurements(ttp, tmp, s, ports, ports,
		NULL) == -1) {
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
		    tmp->tm_a_matrix, tmp->tm_a_rows, tmp->tm_b_columns,
		    tmp->tm_b_matrix, ports, ports, vdp) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    } else {
	if (vnacal_apply_m(vcp, ci, ttp->tt_frequency_vector,
		    ttp->tt_frequencies, tmp->tm_b_matrix,
		    ports, ports, vdp) == -1) {
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
		    v = _vnacal_get_parameter_value_i(vpmrp, f);

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
		expected = _vnacal_get_parameter_value_i(vpmrp, f);
		actual = vnadata_get_cell(vdp, findex, s_row, s_column);
		if (!test_isequal(actual, expected)) {
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
    test_vnacal_free_measurements(tmp);
    test_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_apply: test vnacal_apply
 */
static test_result_t test_vnacal_apply()
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
    test_result_t result = T_FAIL;

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
    test_report(result);
    return result;
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
    exit(99);
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
	}
	break;
    }
    test_init_isequal();
    exit(test_vnacal_apply());
}

#if 0
/************************************************************
 * Code needed for vnacal_map_* functions (yet to be written)
 * that replace the older (broken) vnacal_apply_t functions.
 * This will go into a different file eventually, but park it
 * here for now.
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
static test_result_t test_vnacal_apply_helper(int trial, int frequencies,
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
    test_result_t result = T_FAIL;
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
	(void)printf("Test vnacal_opply: trial %3d cal size (%d x %d) "
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
    if ((error_terms = test_vnacal_generate_error_terms(cmsp)) == NULL) {
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
		    vnacal_terms_t *etp;

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
		actual_matrix[row * dcolumns + column][findex] = test_crandn();
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
		    test_print_cmatrix(a, vrows, vrows);
		    (void)printf("b:\n");
		    test_print_cmatrix(b, vrows, 1);
		}
		/*
		 * Find x = A^-1 b.
		 */
		d = _vnacommon_mldivide(x, a, b, vrows, 1);
		if (cabs(d) <= EPS)	 {
		    (void)fprintf(stderr, "%s: test_vnacal_mrdivide: warning: "
			    "skipping nearly singular test matrix\n",
			    progname);
		    result = T_FAIL;
		    goto out;
		}
		if (opt_v != 0) {
		    (void)printf("x:\n");
		    test_print_cmatrix(x, vrows, 1);
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
	 * vnacal_apply function.
	 */
	if (map == NULL && vrows == drows && vcolumns == dcolumns) {
	    if ((output_matrix = vnadata_alloc()) == NULL)  {
		(void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    if (vnacal_apply(vcp, /*ci*/0, frequencies,
			calp->cal_frequency_vector,
			(const double complex *const *)&M(0, 0),
			output_matrix) == -1) {
		(void)fprintf(stderr, "%s: vnacal_apply: "
			"%s\n", progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    if (opt_v != 0) {
		(void)printf("computed_vector (vnacal_apply):\n");
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
	    test_vnacal_free_error_terms(error_terms, cmsp);
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
static test_result_t test_vnacal_map_apply()
{
    test_result_t result = T_FAIL;
    bool pass = false;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
	for (const apply_test_case_type *atcp = apply_test_cases;
		atcp < &apply_test_cases[N_APPLY_CASES]; ++atcp) {
	    result = test_vnacal_apply_helper(trial, 2, atcp);
	    if (result != T_PASS)
		goto out;
	}
    }
    result = T_PASS;

out:
    test_report(result);
    return result;
}
#endif /* end of mapped vnacal_apply tests */
