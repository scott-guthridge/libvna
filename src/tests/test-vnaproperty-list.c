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
 * test_vnaproperty_list
 */
static libt_result_t test_vnaproperty_list()
{
    vnaproperty_t *root = NULL;
    int type = -1;
    int count;
    libt_result_t result = T_SKIPPED;

    /*
     * Test alloc, type and count of empty list.
     */
    if (vnaproperty_set_subtree(&root, "[]") == NULL) {
	(void)printf("1: vnaproperty_set_subtree: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((type = vnaproperty_type(root, ".")) != 'l') {
	(void)printf("2: vnaproperty_type: 0x%04X != 'l'\n",
		(int)type);
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, ".")) != 0) {
	(void)printf("3: vnaproperty_count: %d != 0\n", count);
	result = T_FAIL;
	goto out;
    }

    /*
     * Test append.
     */
    for (int i = 0; i < 100; ++i) {
	if ((count = vnaproperty_count(root, "[]")) == -1) {
	    (void)printf("4[%d]: vnaproperty_count: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (count != i) {
	    (void)printf("5[%d]: vnaproperty_count: %d != %d\n", i, count, i);
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_set(&root, "[+]=%d", i) == -1) {
	    (void)printf("6[%d]: vnaproperty_set: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
    }
    if ((count = vnaproperty_count(root, ".")) == -1) {
	(void)printf("7: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 100) {
	(void)printf("8: vnaproperty_count: %d != %d\n", count, 100);
	result = T_FAIL;
	goto out;
    }

    /*
     * Test get.
     */
    for (int i = 0; i < 100; ++i) {
	const char *value;

	errno = 0;
	if ((value = vnaproperty_get(root, "[%d]", i)) == NULL) {
	    (void)printf("10[%d]: vnaproperty_get: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (atoi(value) != i) {
	    (void)printf("11[%d]: value %s != %d\n", i, value, i);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test set with invalid index.
     */
    if (vnaproperty_set(&root, "[-1]=invalid") != -1) {
	(void)printf("20: expected set out of bounds to fail\n");
	result = T_FAIL;
	goto out;
    }
    if (errno != EINVAL) {
	(void)printf("21: %s: expected EINVAL\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Test set in the middle.
     *   starting state: 0..99
     *   ending state:   0..49 "fifty" 51..99
     */
    if (vnaproperty_set(&root, "[50]=fifty") == -1) {
	(void)printf("30: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, "[]")) == -1) {
	(void)printf("31: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 100) {
	(void)printf("32: vnaproperty_count: %d != %d\n", count, 100);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 100; ++i) {
	const char *value;

	errno = 0;
	if ((value = vnaproperty_get(root, "[%d]", i)) == NULL) {
	    (void)printf("33[%d]: vnaproperty_get: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (i != 50) {
	    if (atoi(value) != i) {
		(void)printf("34[%d]: value %s != %d\n", i, value, i);
		result = T_FAIL;
		goto out;
	    }
	} else if (strcmp(value, "fifty") != 0) {
	    (void)printf("35[%d]: value %s != fifty\n", i, value);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test setting past the end.
     *   starting state: 0..49 "fifty" 51..99
     *   ending state:   0..49 "fifty" 51..99 ~ ~ "hundred-two"
     */
    if (vnaproperty_set(&root, "[102]=hundred-two") == -1) {
	(void)printf("40: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, "[]")) == -1) {
	(void)printf("41: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 103) {
	(void)printf("42: vnaproperty_count: %d != %d\n", count, 103);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 104; ++i) {
	const char *value;
	char expected[50];

	errno = 0;
	value = vnaproperty_get(root, "[%d]", i);
	if (i == 100 || i == 101 || i == 103) {
	    if (value != NULL) {
		(void)printf("43[%d]: expected NULL; found \"%s\"\n",
			i, value);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 103) {
		if (errno != ENOENT) {
		    (void)printf("44[%d]: %s: expected ENOENT\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    } else {
		if (errno != 0) {
		    (void)printf("45[%d]: %s: expected no error\n",
			i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	    continue;
	}
	if (value == NULL) {
	    (void)printf("46[%d]: vnaproperty_get: %s\n",
		    i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	switch (i) {
	case 50:
	    (void)strcpy(expected, "fifty");
	    break;
	case 102:
	    (void)strcpy(expected, "hundred-two");
	    break;
	default:
	    (void)sprintf(expected, "%d", i);
	    break;
	}
	if (strcmp(value, expected) != 0) {
	    (void)printf("47[%d]: \"%s\" != \"%s\"\n", i, value, expected);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test insert in the middle.
     *	 starting state: 0..49 "fifty" 51..99 ~ ~ "hundred-two"
     *   ending state: 0..50 [51]="fifty" [52]=51..[100]=99
     *			     [101]=~ [102]=~ [103]="hundred-two"
     */
    if (vnaproperty_set(&root, "[50+]=50") == -1) {
	(void)printf("50: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, ".")) == -1) {
	(void)printf("51: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 104) {
	(void)printf("52: vnaproperty_count: %d != %d\n", count, 104);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 105; ++i) {
	const char *value;
	char expected[50];

	errno = 0;
	value = vnaproperty_get(root, "[%d]", i);
	if (i == 101 || i == 102 || i == 104) {
	    if (value != NULL) {
		(void)printf("53[%d]: expected NULL; found \"%s\"\n",
			i, value);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 104) {
		if (errno != ENOENT) {
		    (void)printf("54[%d]: %s: expected ENOENT\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    } else {
		if (errno != 0) {
		    (void)printf("55[%d]: %s: expected no error\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	    continue;
	}
	if (value == NULL) {
	    (void)printf("56[%d]: vnaproperty_get: %s\n",
		    i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	switch (i) {
	case 51:
	    (void)strcpy(expected, "fifty");
	    break;
	case 103:
	    (void)strcpy(expected, "hundred-two");
	    break;
	default:
	    (void)sprintf(expected, "%d", i <= 51 ? i : i - 1);
	    break;
	}
	if (strcmp(value, expected) != 0) {
	    (void)printf("57[%d]: \"%s\" != \"%s\"\n", i, value, expected);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test insert at the end.
     *   starting state: 0..50 [51]="fifty" [52]=51..[100]=99
     *	[101]=~ [102]=~ [103]="hundred-two"
     *   ending state:   0..50 [51]="fifty" [52]=51..[100]=99
     *	[101]=~ [102]=~ [103]="hundred-two" [104]="one-o-four"
     */
    if (vnaproperty_set(&root, "[104+]=one-o-four") == -1) {
	(void)printf("60: vnaproperty_set: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, ".")) == -1) {
	(void)printf("61: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 105) {
	(void)printf("62: vnaproperty_count: %d != %d\n", count, 105);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 106; ++i) {
	const char *value;
	char expected[50];

	errno = 0;
	value = vnaproperty_get(root, "[%d]", i);
	if (i == 101 || i == 102 || i == 105) {
	    if (value != NULL) {
		(void)printf("63[%d]: expected NULL; found \"%s\"\n",
			i, value);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 105) {
		if (errno != ENOENT) {
		    (void)printf("64[%d]: %s: expected ENOENT\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    } else {
		if (errno != 0) {
		    (void)printf("65[%d]: %s: expected no error\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	    continue;
	}
	if (value == NULL) {
	    (void)printf("66[%d]: vnaproperty_get: %s\n",
		    i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	switch (i) {
	case 51:
	    (void)strcpy(expected, "fifty");
	    break;
	case 103:
	    (void)strcpy(expected, "hundred-two");
	    break;
	case 104:
	    (void)strcpy(expected, "one-o-four");
	    break;
	default:
	    (void)sprintf(expected, "%d", i <= 51 ? i : i - 1);
	    break;
	}
	if (strcmp(value, expected) != 0) {
	    (void)printf("67[%d]: \"%s\" != \"%s\"\n", i, value, expected);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test delete in the middle.
     *   starting state: 0..50 [51]="fifty" [52]=51..[100]=99
     *	   [101]=~ [102]=~ [103]="hundred-two" [104]="one-o-four"
     *   ending state:   0..99 [100]=~ [101]=~ [102]="hundred-two"
     *     [103]="one-o-four"
     */
    if (vnaproperty_delete(&root, "[51]") == -1) {
	(void)printf("70: vnaproperty_delete: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, ".")) == -1) {
	(void)printf("71: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 104) {
	(void)printf("72: vnaproperty_count: %d != %d\n", count, 104);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 105; ++i) {
	const char *value;
	char expected[50];

	errno = 0;
	value = vnaproperty_get(root, "[%d]", i);
	if (i == 100 || i == 101 || i == 104) {
	    if (value != NULL) {
		(void)printf("73[%d]: expected NULL; found \"%s\"\n",
			i, value);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 104) {
		if (errno != ENOENT) {
		    (void)printf("74[%d]: %s: expected ENOENT\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    } else {
		if (errno != 0) {
		    (void)printf("75[%d]: %s: expected no error\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	    continue;
	}
	if (value == NULL) {
	    (void)printf("76[%d]: vnaproperty_get: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	switch (i) {
	case 102:
	    (void)strcpy(expected, "hundred-two");
	    break;
	case 103:
	    (void)strcpy(expected, "one-o-four");
	    break;
	default:
	    (void)sprintf(expected, "%d", i);
	    break;
	}
	if (strcmp(value, expected) != 0) {
	    (void)printf("77[%d]: \"%s\" != \"%s\"\n", i, value, expected);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test delete at the end.
     *   starting state: 0..99 [100]=~ [101]=~ [102]="hundred-two"
     *     [103]="one-o-four"
     *   ending state:   0..99 [100]=~ [101]=~ [102]="hundred-two"
     */
    if (vnaproperty_delete(&root, "[103]") == -1) {
	(void)printf("80: vnaproperty_delete: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, ".")) == -1) {
	(void)printf("81: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != 103) {
	(void)printf("82: vnaproperty_count: %d != %d\n", count, 104);
	result = T_FAIL;
	goto out;
    }
    for (int i = 0; i < 105; ++i) {
	const char *value;
	char expected[50];

	errno = 0;
	value = vnaproperty_get(root, "[%d]", i);
	if (i == 100 || i == 101 || i == 103 || i == 104) {
	    if (value != NULL) {
		(void)printf("83[%d]: expected NULL; found \"%s\"\n",
			i, value);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 103 || i == 104) {
		if (errno != ENOENT) {
		    (void)printf("84[%d]: %s: expected ENOENT\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    } else {
		if (errno != 0) {
		    (void)printf("85[%d]: %s: expected no error\n",
			    i, strerror(errno));
		    result = T_FAIL;
		    goto out;
		}
	    }
	    continue;
	}
	if (value == NULL) {
	    (void)printf("86[%d]: vnaproperty_get: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	switch (i) {
	case 102:
	    (void)strcpy(expected, "hundred-two");
	    break;
	default:
	    (void)sprintf(expected, "%d", i);
	    break;
	}
	if (strcmp(value, expected) != 0) {
	    (void)printf("87[%d]: \"%s\" != \"%s\"\n", i, value, expected);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test delete all.
     */
    if (vnaproperty_delete(&root, "[]") == -1) {
	(void)printf("90: vnaproperty_delete: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (root != NULL) {
	(void)printf("91: expected NULL after delete .\n");
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
    exit(test_vnaproperty_list());
}
