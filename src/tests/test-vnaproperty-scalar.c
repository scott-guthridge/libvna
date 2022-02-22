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
#include "vnaproperty.h"
#include "libt.h"


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
bool opt_a = false;
int opt_v = 0;

/*
 * test_vnaproperty_scalar
 */
static libt_result_t test_vnaproperty_scalar()
{
    vnaproperty_t *root = NULL;
    int type = -1;
    const char *value;
    static const char text1[] = "abcdefghijklmnopqrstuvwxyz";
    static const char text2[] = "~";
    libt_result_t result = T_SKIPPED;

    if (vnaproperty_set(&root, ".=%s", text1) == -1) {
	(void)printf("1: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((type = vnaproperty_type(root, ".")) != 's') {
	(void)printf("2: vnaproperty_type: 0x%04X != 's'\n",
		(int)type);
	result = T_FAIL;
	goto out;
    }
    errno = 0;
    if ((value = vnaproperty_get(root, ".")) == NULL) {
	(void)printf("3: vnaproperty_get: NULL: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, text1) != 0) {
	(void)printf("3: vnaproperty_get: mismatch: \"%s\" != \"%s\"",
		value, text1);
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_set(&root, ".#") == -1) {
	(void)printf("5: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (root != NULL) {
	(void)printf("6: root not NULL after set .#\n");
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_set(&root, ".=%s", text2) == -1) {
	(void)printf("7: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((type = vnaproperty_type(root, ".")) != 's') {
	(void)printf("8: vnaproperty_type: 0x%04X != 's'\n",
		(int)type);
	result = T_FAIL;
	goto out;
    }
    errno = 0;
    if ((value = vnaproperty_get(root, ".")) == NULL) {
	(void)printf("9: vnaproperty_get: NULL: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (strcmp(value, text2) != 0) {
	(void)printf("10: vnaproperty_get: mismatch: \"%s\" != \"%s\"",
		value, text2);
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_delete(&root, ".") == -1) {
	(void)printf("11: vnaproperty_delete: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (root != NULL) {
	(void)printf("12: root not NULL after delete .\n");
	result = T_FAIL;
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
    exit(test_vnaproperty_scalar());
}
