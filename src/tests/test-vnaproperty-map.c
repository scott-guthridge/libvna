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
 * collection of words randomly chosen from /usr/share/dict/words
 * plus a few special cases
 */
static const char *words[] = {
    "done",
    "unbrilliantly",
    "Sextonville",
    "seconal",
    "rock-bestudded",
    "preorganically",
    "Praxitelean",
    "neurotoxia",
    "suisimilar",
    "outgives",
    "insidiation",
    "proadoption",
    "prepontine",
    "sororize",
    "ZZZ",
    "preestimates",
    "cognatus",
    "Bundaberg",
    "Ennosigaeus",
    "postcommunion",
    "Cardin",
    "fanaticalness",
    "zoisite",
    "prospeculation",
    "fillock",
    "oreman",
    "nimming",
    "Wattenscheid",
    "imitator",
    "Evert",
    "tropaeolaceous",
    "This is a phrase\\.",
    "\\[specials and trailing spaces\\]\\ \\ "
};
#define N_WORDS		(sizeof(words) / sizeof(char *))

/*
 * test_vnaproperty_map
 */
static libt_result_t test_vnaproperty_map()
{
    vnaproperty_t *root = NULL;
    int type = -1;
    int count, rv;
    const char **keys = NULL;
    libt_result_t result = T_SKIPPED;

    /*
     * Test alloc, type, count and keys of an empty map.
     */
    if (vnaproperty_set_subtree(&root, "{}") == NULL) {
	(void)printf("1: vnaproperty_set_subtree: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((type = vnaproperty_type(root, ".")) != 'm') {
	(void)printf("2: vnaproperty_type: 0x%04X != 'm'\n",
		(int)type);
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, ".")) != 0) {
	(void)printf("3: vnaproperty_count: %d != 0\n", count);
	result = T_FAIL;
	goto out;
    }
    if ((keys = vnaproperty_keys(root, "{}")) == NULL) {
	(void)printf("4: vnaproperty_keys: returned NULL\n");
	result = T_FAIL;
	goto out;
    }
    if (keys[0] != NULL) {
	(void)printf("5: keys[0] (%s) != NULL\n", keys[0]);
	result = T_FAIL;
	goto out;
    }
    free((void *)keys);
    keys = NULL;

    /*
     * Test set.
     */
    for (int i = 0; i < N_WORDS; ++i) {
	if (vnaproperty_set(&root, "%s=%d", words[i], i) == -1) {
	    (void)printf("10[%d]: vnaproperty_set: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if ((count = vnaproperty_count(root, ".")) == -1) {
	    (void)printf("11[%d]: vnaproperty_count: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (count != i + 1) {
	    (void)printf("12[%d]: vnaproperty_count: %d != %d\n",
		    i, count, i + 1);
	    result = T_FAIL;
	    goto out;
	}
    }
    if ((count = vnaproperty_count(root, ".")) == -1) {
	(void)printf("13: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != N_WORDS) {
	(void)printf("14: vnaproperty_count: %d != %d\n", count, (int)N_WORDS);
	result = T_FAIL;
	goto out;
    }

    /*
     * Test get.
     */
    for (int i = 0; i < N_WORDS; ++i) {
	const char *value;

	errno = 0;
	if ((value = vnaproperty_get(root, "%s", words[i])) == NULL) {
	    (void)printf("20[%d]: vnaproperty_get: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (atoi(value) != i) {
	    (void)printf("21[%d]: vnaproperty_get: %s != %d\n", i, value, i);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test get of non-existent.
     */
    errno = 0;
    if (vnaproperty_get(root, "NotInList") != NULL || errno != ENOENT) {
	(void)printf("30: vnaproperty_get: %s (NotInList)\n",
		strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Test change via set.
     */
    for (int i = N_WORDS-1; i >= 0; --i) {
	if (vnaproperty_set(&root, "%s=%d", words[i], -i) == -1) {
	    (void)printf("40[%d]: vnaproperty_set: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if ((count = vnaproperty_count(root, ".")) == -1) {
	    (void)printf("41: vnaproperty_count: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (count != N_WORDS) {
	    (void)printf("42[%d]: vnaproperty_count: %d != %d\n", i, count, i);
	    result = T_FAIL;
	    goto out;
	}
    }
    for (int i = 0; i < N_WORDS; ++i) {
	const char *value;

	errno = 0;
	if ((value = vnaproperty_get(root, "%s", words[i])) == NULL) {
	    (void)printf("43[%d]: vnaproperty_get: %s\n", i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (atoi(value) != -i) {
	    (void)printf("44[%d]: vnaproperty_get: %s != %d\n", i, value, -i);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test delete by deleting all the odd words.
     */
    for (int i = 0; i < N_WORDS / 2; ++i) {
	if (vnaproperty_delete(&root, "%s", words[2 * i + 1]) == -1) {
	    (void)printf("50[%d]: vnaproperty_delete: %s\n",
		    i, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if ((count = vnaproperty_count(root, "{}")) == -1) {
	    (void)printf("51: vnaproperty_count: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (count != N_WORDS - i - 1) {
	    (void)printf("52[%d]: vnaproperty_count: %d != %d\n", i, count, i);
	    result = T_FAIL;
	    goto out;
	}
    }
    for (int i = 0; i < N_WORDS; ++i) {
	const char *value;

	errno = 0;
	value = vnaproperty_get(root, "%s", words[i]);
	if (i & 1) {
	    if (value != NULL) {
		(void)printf("53[%d]: deleted element \"%s\" should be NULL\n",
			i, value);
		result = T_FAIL;
		goto out;
	    }
	    if (errno != ENOENT) {
		(void)printf("54[%d]: %s: errno should be ENOENT\n",
			i, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    continue;
	}
	if (value == NULL) {
	    (void)printf("55[%d]: vnaproperty_get: unexpected NULL\n", i);
	    result = T_FAIL;
	    goto out;
	}
	if (atoi(value) != -i) {
	    (void)printf("56[%d]: vnaproperty_get: %s != %d\n", i, value, -i);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test delete of non-existent.
     */
    errno = 0;
    if ((rv = vnaproperty_delete(&root, "NotInList")) != -1) {
	(void)printf("60: delete of non-existent returned %d", rv);
	result = T_FAIL;
	goto out;
    }
    if (errno != ENOENT) {
	(void)printf("61: %s: errno should be ENOENT", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if ((count = vnaproperty_count(root, "{}")) == -1) {
	(void)printf("62: vnaproperty_count: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (count != (N_WORDS + 1) / 2) {
	(void)printf("63: vnaproperty_count: %d != %d\n",
		count, (int)(N_WORDS + 1) / 2);
	result = T_FAIL;
	goto out;
    }
    result = T_PASS;

    /*
     * Test keys and quote_key.
     */
    if ((keys = vnaproperty_keys(root, ".")) == NULL) {
	(void)printf("70: vnaproperty_keys: returned NULL\n");
	result = T_FAIL;
	goto out;
    }
    count = 0;
    for (const char **cpp = keys; *cpp != NULL; ++cpp, ++count) {
	char *quoted = vnaproperty_quote_key(*cpp);

	if (strcmp(quoted, words[2 * count]) != 0) {
	    (void)printf("71[%d]: key \"%s\" != \"%s\"\n",
		count, quoted, words[2 * count]);
	    free((void *)quoted);
	    result = T_FAIL;
	    goto out;
	}
	free((void *)quoted);
    }
    if (count != (N_WORDS + 1) / 2) {
	(void)printf("72: vnaproperty_keys returned only %d of %d keys\n",
		count, (int)(N_WORDS + 1) / 2);
	result = T_FAIL;
	goto out;
    }
    free((void *)keys);

out:
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
    exit(test_vnaproperty_map());
}
