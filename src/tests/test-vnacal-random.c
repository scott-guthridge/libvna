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

#include "archdep.h"

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
#include "vnacal_internal.h"
#include "libt.h"
#include "libt_vnacal.h"


#define NTRIALS		67

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
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)fprintf(stderr, "%s: %s\n", progname, message);
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
static libt_result_t run_vnacal_new_random_trial(int trial,
	vnacal_type_t type, int m_rows, int m_columns,
	int frequencies, bool ab)
{
    vnacal_t *vcp = NULL;
    libt_vnacal_terms_t *ttp = NULL;
    bool add_all_match;
    libt_result_t result = T_FAIL;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	int standards;

	standards = libt_vnacal_calc_needed_standards(type, m_rows, m_columns,
		&add_all_match);
	(void)printf("Test vnacal_new: trial %3d size %d x %d "
		"type %-4s %s %2d random standards%s\n",
		trial, m_rows, m_columns, vnacal_type_to_name(type),
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
    if ((ttp = libt_vnacal_make_random_calibration(vcp, type, m_rows, m_columns,
		    frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    libt_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_random: test vnacal_new_* with random multi-port standards
 */
static libt_result_t test_vnacal_new_random()
{
    static const int sizes[] = { 1, 2, 3, 4 };
    static const vnacal_type_t types[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10, VNACAL_T16, VNACAL_U16,
	VNACAL_UE14, VNACAL_E12
    };
    libt_result_t result = T_FAIL;

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
    libt_report(result);
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
    libt_isequal_init();
    exit(test_vnacal_new_random());
}
