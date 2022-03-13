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
 * Test Strings for vnacal_property_set
 */
static const char property_foo_value[] = "1234567890";
static const char property_bar_value[] = "abcdefghijkl\nmnopqrstuvwxyz";
static const char property3_value[] = "αβγδεζηθικλμνξοπρστυφχψω";

/*
 * run_vnacal_save_load_trial
 */
static libt_result_t run_vnacal_save_load_trial(int trial)
{
    static const int dimension_table[][2] = {
	{ 1, 1 }, { 1, 2 }, { 1, 3 }, { 1, 4 }, { 2, 2 },
	{ 2, 3 }, { 2, 4 }, { 3, 3 }, { 3, 4 }, { 4, 4 },
    };
    static const vnacal_type_t type_table[] = {
	VNACAL_T8, VNACAL_U8, VNACAL_TE10, VNACAL_UE10,
	VNACAL_T16, VNACAL_U16, VNACAL_UE14, VNACAL_E12
    };
    libt_vnacal_terms_t *ttp_table[8] =
        { NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL };
    vnacal_t *vcp = NULL;
    const int types = sizeof(type_table) / sizeof(vnacal_type_t);
    const int dimensions = sizeof(dimension_table) / sizeof(int [2]);
    libt_result_t result = T_FAIL;
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

	default:
	    abort();
	}
	if ((ttp_table[tindex] = libt_vnacal_make_random_calibration(vcp, type,
			m_rows, m_columns, frequencies, /*ab*/false)) == NULL) {
	    result = T_FAIL;
	    goto out;
	}
	if (vnacal_add_calibration(vcp, vnacal_type_to_name(type),
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
    if (vnacal_save(vcp, "test-vnacal.vnacal") == -1) {
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
    if ((vcp = vnacal_load("test-vnacal.vnacal", error_fn, NULL)) == NULL) {
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
			vnacal_type_to_name(type))) == -1) {
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
	if (libt_vnacal_validate_calibration(ttp_table[tindex], calp) == -1) {
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
			    "switches[%d][%d] in calibration 1 not found: %s\n",
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
	libt_vnacal_free_error_terms(ttp_table[tindex]);
    }
    vnacal_free(vcp);
    return result;
}

/*
 * test_vnacal_save_load
 */
static libt_result_t test_vnacal_save_load()
{
    libt_result_t result = T_FAIL;

    for (int trial = 0; trial < 5; ++trial) {
	result = run_vnacal_save_load_trial(trial);
	if (result != T_PASS)
	    goto out;
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
    if (libt_isequal_eps < 0.00001) {	/* save uses 6 digits by default */
	libt_isequal_eps = 0.00001;
    }
    exit(test_vnacal_save_load());
}
