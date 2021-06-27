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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnacal_internal.h"
#include "test.h"
#include "vnacaltest.h"


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
 *   @category: error category
 *   @message: error message
 *   @arg: (unused)
 */
static void error_fn(vnaerr_category_t category, const char *message, void *arg)
{
    (void)fprintf(stderr, "%s: %s\n", progname, message);
}

/*
 * run_solt_trial_helper: add short, open and load calibrations on port
 *   @ttp: pointer to test error terms structure
 *   @tmp: test measurements structure
 *   @port: port to calibrate
 */
static test_result_t run_solt_trial_helper(test_vnacal_terms_t *ttp,
	test_vnacal_measurements_t *tmp, int port)
{
    if (test_vnacal_add_single_reflect(ttp, tmp, VNACAL_SHORT, port) == -1) {
	return T_FAIL;
    }
    if (test_vnacal_add_single_reflect(ttp, tmp, VNACAL_OPEN, port) == -1) {
	return T_FAIL;
    }
    if (test_vnacal_add_single_reflect(ttp, tmp, VNACAL_MATCH, port) == -1) {
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
static test_result_t run_vnacal_new_solt_trial(int trial, vnacal_type_t type,
	int rows, int columns, int frequencies, bool ab)
{
    vnacal_t *vcp = NULL;
    test_vnacal_terms_t *ttp = NULL;
    test_vnacal_measurements_t *tmp = NULL;
    int diagonals = MIN(rows, columns);
    int ports = MAX(rows, columns);
    test_result_t result = T_FAIL;

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
    if ((ttp = test_vnacal_generate_error_terms(vcp, type, rows, columns,
		    frequencies, NULL, 1.0, ab)) == NULL) {
	(void)fprintf(stderr, "%s: test_vnacal_generate_error_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate the test measurement matrices.
     */
    if ((tmp = test_vnacal_alloc_measurements(type, rows, columns,
		    frequencies, ab)) == NULL) {
	result = T_ERROR;
	goto out;
    }

    /*
     * For E12 and UE14, we have to do short, open and load calibration
     * on every diagonal port.  For the others, we can choose any one
     * diagonal port.
     */
    if (type == VNACAL_E12 || VNACAL_IS_UE14(type)) {
	for (int port = 1; port <= diagonals; ++port) {
	    result = run_solt_trial_helper(ttp, tmp, port);
	    if (result != T_PASS)
		goto out;
	}
    } else {
	int port = random() % diagonals + 1;

	result = run_solt_trial_helper(ttp, tmp, port);
	if (result != T_PASS)
	    goto out;
    }

    /*
     * Do through tests between every diagonal port and every other port.
     */
    for (int port1 = 1; port1 <= diagonals; ++port1) {
	for (int port2 = port1 + 1; port2 <= ports; ++port2) {
	    if (test_vnacal_add_through(ttp, tmp, port1, port2) == -1) {
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
    if (test_vnacal_validate_calibration(ttp, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    test_vnacal_free_measurements(tmp);
    test_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_solt: run SOLT tests for 8-12 term parameters
 */
static test_result_t test_vnacal_new_solt()
{
    static const int sizes[] = { 1, 2, 3, 4 };
    static const vnacal_type_t types[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10, VNACAL_UE14, VNACAL_E12
    };
    test_result_t result = T_FAIL;

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
    exit(test_vnacal_new_solt());
}
