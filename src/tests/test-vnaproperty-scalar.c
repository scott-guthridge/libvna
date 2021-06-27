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
 * test_vnaproperty_scalar
 */
static test_result_t test_vnaproperty_scalar()
{
    vnaproperty_t *scalar;
    const char *value;
    static const char text1[] = "abcdefghijklmnopqrstuvwxyz";
    static const char text2[] = "0123456789";
    char buf[64];
    test_result_t result = T_SKIPPED;

    (void)strcpy(buf, text1);
    scalar = vnaproperty_scalar_alloc(buf);
    if (scalar == NULL) {
	(void)printf("vnaproperty_scalar_alloc: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
	(void)printf("vnaproperty_type(scalar) != VNAPROPERTY_SCALAR "
		"(1)\n");
	result = T_FAIL;
	goto out;
    }
    (void)strcpy(buf, text2);
    value = vnaproperty_scalar_get(scalar);
    if (value == NULL) {
	(void)printf("vnaproperty_get_value(scalar) == NULL (1)");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, text1) != 0) {
	(void)printf("vnaproperty_get_value(scalar) ne text1");
	result = T_FAIL;
	goto out;
    }
    if (vmaproperty_scalar_set(scalar, buf) == -1) {
	(void)printf("vmaproperty_scalar_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    value = vnaproperty_scalar_get(scalar);
    if (value == NULL) {
	(void)printf("vnaproperty_get_value(scalar) == NULL (2)");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, text2) != 0) {
	(void)printf("vnaproperty_get_value(scalar) ne text2 (1)");
	result = T_FAIL;
	goto out;
    }
    vnaproperty_hold(scalar);
    vnaproperty_free(scalar);
    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
	(void)printf("vnaproperty_type(scalar) != VNAPROPERTY_SCALAR "
		"(2)\n");
	result = T_FAIL;
	goto out;
    }
    if (value == NULL) {
	(void)printf("vnaproperty_get_value(scalar) == NULL (3)");
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, text2) != 0) {
	(void)printf("vnaproperty_get_value(scalar) ne text2 (2)");
	result = T_FAIL;
	goto out;
    }
    vnaproperty_free(scalar);
    if (vnaproperty_type(scalar) == VNAPROPERTY_SCALAR) {
	(void)printf("still a scalar after free!");
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

out:
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
    exit(test_vnaproperty_scalar());
}
