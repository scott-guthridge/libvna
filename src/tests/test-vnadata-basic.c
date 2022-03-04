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
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnadata.h"
#include "libt.h"
#include "libt_vnadata.h"


#define N_TRIALS	5

/*
 * Options
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
int opt_v = 0;

/*
 * error_fn: error reporting function
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)printf("error: %s: %s\n", progname, message);
}



/*
 * run_trial
 */
static int run_trial(int trial, vnadata_t *vdp, vnadata_parameter_type_t type,
	int rows, int columns, int frequencies, libt_vnadata_z0_type_t z0_type,
	libt_vnadata_fill_method_t fill_method)
{
    libt_vnadata_t *tdp = NULL;
    libt_result_t result = T_FAIL;

    if (opt_v >= 1) {
	const char *type_name = vnadata_get_type_name(type);

	(void)printf("Test vndata_basic: trial %2d type %-3s size %d x %d "
		"f %d %s %s\n",
		trial, (type == VPT_UNDEF) ? "-" : type_name,
		rows, columns, frequencies, libt_vnadata_z0_names[z0_type],
		libt_vnadata_fill_names[fill_method]);
	(void)fflush(stdout);
    }
    if ((tdp = libt_vnadata_create(type, rows, columns,
		    frequencies, z0_type)) == NULL) {
	result = T_FAIL;
	goto out;
    }
    if ((result = libt_vnadata_fill(tdp, vdp, fill_method)) != T_PASS) {
	goto out;
    }
    if ((result = libt_vnadata_validate(tdp, vdp)) != T_PASS) {
	goto out;
    }
    result = T_PASS;

out:
    libt_vnadata_free(tdp);
    return result;
}

/*
 * test_vnadata_basic_helper: run inner loops
 */
static libt_result_t test_vnadata_basic_helper(int trial,
	vnadata_t *vdp, vnadata_parameter_type_t type, int rows, int columns)
{
    const int frequency_list[] = { 0, 1, 2, 10 };
    libt_vnadata_z0_type_t z0_type;
    libt_vnadata_fill_method_t fill_method;
    libt_result_t result = T_SKIPPED;

    for (int i = 0; i < sizeof(frequency_list) / sizeof(int); ++i) {
	for (z0_type = 0; z0_type < Z0_NTYPES; ++z0_type) {
	    for (fill_method = 0; fill_method < FM_NMETHODS; ++fill_method) {
		result = run_trial(trial, vdp, type, rows, columns,
			frequency_list[i], z0_type, fill_method);
		if (result != T_PASS) {
		    return result;
		}
	    }
	}
    }
    return T_PASS;
}

/*
 * test_vnadata_basic: run basic tests on vnadata
 */
static libt_result_t test_vnadata_basic()
{
    libt_result_t result = T_SKIPPED;
    vnadata_t *vdp = NULL;

    /*
     * Allocate the vndata structure.  We'll use the same structure
     * through all the trials to make sure init/resize works.
     */
    if ((vdp = vnadata_alloc(error_fn, NULL)) == NULL) {
	libt_fail("vnadata_alloc: returned NULL\n");
	result = T_FAIL;
	goto out;
    }
    for (int trial = 0; trial < N_TRIALS; ++trial) {
	for (vnadata_parameter_type_t type = 0; type < VPT_NTYPES; ++type) {
	    switch (type) {
	    case VPT_UNDEF:
	    case VPT_S:
		for (int rows = 0; rows < 10; ++rows) {
		    for (int columns = 0; columns < 10; ++columns) {
			result = test_vnadata_basic_helper(trial, vdp, type,
				rows, columns);
			if (result != T_PASS) {
			    goto out;
			}
		    }
		}
		break;

	    case VPT_Z:
	    case VPT_Y:
		for (int ports = 0; ports < 10; ++ports) {
		    result = test_vnadata_basic_helper(trial, vdp, type,
			    ports, ports);
		    if (result != T_PASS) {
			goto out;
		    }
		}
		break;

	    case VPT_H:
	    case VPT_G:
	    case VPT_A:
	    case VPT_B:
	    case VPT_T:
		result = test_vnadata_basic_helper(trial, vdp, type, 2, 2);
		if (result != T_PASS) {
		    goto out;
		}
		break;

	    case VPT_ZIN:
		for (int ports = 0; ports < 10; ++ports) {
		    result = test_vnadata_basic_helper(trial, vdp,
			    type, 1, ports);
		    if (result != T_PASS) {
			goto out;
		    }
		}
		break;

	    default:
		assert(!"unhandled case in switch");
	    }
	}
    }
    result = T_PASS;

out:
    vnadata_free(vdp);
    libt_report(result);;
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
    exit(2);
}

/*
 * main: test program
 */
    int
main(int argc, char **argv)
{
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	case -1:
	    break;

	default:
	    print_usage();
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
    }
    libt_isequal_init();
    exit(test_vnadata_basic());
}
