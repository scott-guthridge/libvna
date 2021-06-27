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
#include <errno.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnaproperty_internal.h"
#include "test.h"


/*
 * Options
 */
char *progname;
static const char options[] = "v";
static const char *const usage[] = {
    "[-v]",
    NULL
};
static const char *const help[] = {
    "-v	 show verbose output",
    NULL
};
static bool opt_v = false;

/*
 * test_vnaproperty_expr
 */
static test_result_t test_vnaproperty_expr()
{
    vnaproperty_t *root = NULL;
    vnaproperty_type_t type;
    int count;
    const char *value;
    const char **keys;
    test_result_t result = T_SKIPPED;

    if (vnaproperty_expr_set(&root, ".=scalar-only") == -1) {
	(void)printf("%s: vnaproperty_expr_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((value = vnaproperty_expr_get(root, ".")) == NULL) {
	(void)printf("%s: vnaproperty_expr_get: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, "scalar-only") != 0) {
	(void)printf("%s: expected value \"scalar-only\", found \"%s\"\n",
		progname, value);
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_expr_set(&root, "A=B") == -1) {
	(void)printf("%s: vnaproperty_expr_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((value = vnaproperty_expr_get(root, "A")) == NULL) {
	(void)printf("%s: vnaproperty_expr_get: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, "B") != 0) {
	(void)printf("%s: expected value \"B\", found \"%s\"\n",
		progname, value);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 3; ++i) {
	for (int j = 0; j < 4; ++j) {
	    if (vnaproperty_expr_set(&root, "matrix[%d][%d]=%d,%d",
			i, j, i, j) == -1) {
		(void)printf("%s: vnaproperty_expr_set: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	}
    }
    if (vnaproperty_expr_set(&root, "foo.bar=baz") == -1) {
	(void)printf("%s: vnaproperty_expr_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((type = vnaproperty_expr_type(root, ".")) == -1) {
	(void)printf("%s: vnaproperty_expr_type: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (type != VNAPROPERTY_MAP) {
	(void)printf("%s: expected type MAP, found type %d\n", progname, type);
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_expr_count(root, ".")) == -1) {
	(void)printf("%s: vnaproperty_expr_count: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 3) {
	(void)printf("%s: expected top-level count 3, found %d\n",
		progname, count);
	result = T_FAIL;
	goto out;
    }
    if ((keys = vnaproperty_expr_keys(root, ".")) == NULL) {
	(void)printf("%s: vnaproperty_expr_keys: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (strcmp(keys[0], "A") != 0) {
	(void)printf("%s: expected first key \"A\" found \"%s\"\n",
		progname, keys[0]);
	result = T_FAIL;
	goto out;
    }
    if (strcmp(keys[1], "matrix") != 0) {
	(void)printf("%s: expected first key \"matrix\" found \"%s\"\n",
		progname, keys[0]);
	result = T_FAIL;
	goto out;
    }
    if (strcmp(keys[2], "foo") != 0) {
	(void)printf("%s: expected first key \"foo\" found \"%s\"\n",
		progname, keys[0]);
	result = T_FAIL;
	goto out;
    }
    if (keys[3] != NULL) {
	(void)printf("%s: expected NULL, found 0x%08lX\n",
		progname, (long)keys[3]);
	result = T_FAIL;
	goto out;
    }
    free((void *)keys);
    keys = NULL;
    result = T_PASS;
    /* delete matrix columns 1 and 3 (zero-based) */
    for (int i = 0; i < 3; ++i) {
	if (vnaproperty_expr_delete(&root, "matrix[%d][3]", i) == -1) {
	    (void)printf("%s: vnaproperty_expr_set: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_expr_delete(&root, "matrix[%d][1]", i) == -1) {
	    (void)printf("%s: vnaproperty_expr_set: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
    }
    /* delete matrix row 1 */
    if (vnaproperty_expr_delete(&root, "matrix[1]") == -1) {
	(void)printf("%s: vnaproperty_expr_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    /* check */
    if ((count = vnaproperty_expr_count(root, "matrix")) == -1) {
	(void)printf("%s: vnaproperty_expr_count: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 2) {
	(void)printf("%s: expected matrix row count 2, found %d\n",
		progname, count);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 2; ++i) {
	if ((count = vnaproperty_expr_count(root, "matrix[%d]", i)) == -1) {
	    (void)printf("%s: vnaproperty_expr_count: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (count != 2) {
	    (void)printf("%s: expected matrix column count 2, found %d\n",
		    progname, count);
	    result = T_FAIL;
	    goto out;
	}
	for (int j = 0; j < 2; ++j) {
	    char buf[6*sizeof(int)+2];

	    if ((value = vnaproperty_expr_get(root, "matrix[%d][%d]",
			    i, j)) == NULL) {
		(void)printf("%s: vnaproperty_expr_get: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    (void)sprintf(buf, "%d,%d", 2*i, 2*j);
	    if (strcmp(value, buf) != 0) {
		(void)printf("%s: expected value \"%s\", found \"%s\"\n",
			progname, buf, value);
		result = T_FAIL;
		goto out;
	    }
	}
    }
    if (vnaproperty_expr_set(&root, "foo[0].bar=zap") == -1) {
	(void)printf("%s: vnaproperty_expr_set: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((type = vnaproperty_expr_type(root, "foo")) == -1) {
	(void)printf("%s: vnaproperty_expr_type: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (type != VNAPROPERTY_LIST) {
	(void)printf("%s: expected type LIST, found type %d\n", progname, type);
	result = T_FAIL;
	goto out;
    }
    if ((value = vnaproperty_expr_get(root, "foo[0].bar")) == NULL) {
	(void)printf("%s: vnaproperty_expr_get: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, "zap") != 0) {
	(void)printf("%s: expected value \"zap\" found \"%s\"\n",
		progname, value);
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
    vnaproperty_free(root);
    root = NULL;
    test_report(result);;
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
	case 'v':
	    opt_v = true;
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
    exit(test_vnaproperty_expr());
}
