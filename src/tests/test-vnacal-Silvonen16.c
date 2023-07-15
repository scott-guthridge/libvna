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
    (void)printf("%s: %s\n", progname, message);
}

/*
 * Calibration Standard Types
 */
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
static libt_result_t test_vnacal_new_table_entry(int trial, vnacal_type_t type,
	int frequencies, const int *table_entry, bool ab)
{
    vnacal_t *vcp = NULL;
    libt_vnacal_terms_t *ttp = NULL;
    libt_vnacal_measurements_t *tmp = NULL;
    libt_result_t result = T_FAIL;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal_new: trial %3d size 2 x 2 type %s %s:",
	    trial, vnacal_type_to_name(type), ab ? "AB" : "M ");
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
    if ((ttp = libt_vnacal_generate_error_terms(vcp, type, 2, 2,
		    frequencies, NULL, 1.0, 0)) == NULL) {
	(void)fprintf(stderr, "%s: libt_vnacal_generate_error_terms: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Allocate the test measurement matrices.
     */
    if ((tmp = libt_vnacal_alloc_measurements(type, 2, 2,
		    frequencies, ab)) == NULL) {
	result = T_FAIL;
	goto out;
    }

    /*
     * Add standards based on the table.
     */
    for (const int *ip = table_entry; *ip != -1; ++ip) {
	switch (*ip) {
	case MM:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_MATCH, VNACAL_MATCH, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case MO:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_MATCH, VNACAL_OPEN, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case MS:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_MATCH, VNACAL_SHORT, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case OM:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_OPEN, VNACAL_MATCH, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case OO:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_OPEN, VNACAL_OPEN, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case OS:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_OPEN, VNACAL_SHORT, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case SM:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_SHORT, VNACAL_MATCH, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case SO:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_SHORT, VNACAL_OPEN, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case SS:
	    if (libt_vnacal_add_double_reflect(ttp, tmp,
			VNACAL_SHORT, VNACAL_SHORT, 1, 2) == -1) {
		result = T_FAIL;
		goto out;
	    }
	    break;

	case T:
	    if (libt_vnacal_add_through(ttp, tmp, 1, 2) == -1) {
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
		if (libt_vnacal_add_line(ttp, tmp, &s[0][0], 1, 2) == -1) {
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
    if (libt_vnacal_validate_calibration(ttp, NULL) == -1) {
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    libt_vnacal_free_measurements(tmp);
    libt_vnacal_free_error_terms(ttp);
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_new_solt: run 16 parameter 2 port tests from Silvonen table
 */
static libt_result_t test_vnacal_new_silvonen16()
{
    static const vnacal_type_t types[] = {
	VNACAL_T16, VNACAL_U16
    };
    libt_result_t result = T_FAIL;

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
    exit(test_vnacal_new_silvonen16());
}
