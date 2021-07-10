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
    "tropaeolaceous"
};
#define N_WORDS		(sizeof(words) / sizeof(char *))

/*
 * test_vnaproperty_map
 */
static libt_result_t test_vnaproperty_map()
{
    vnaproperty_t *map;
    int length;
    vnaproperty_t *first_scalar = NULL;
    libt_result_t result = T_SKIPPED;

    /*
     * Test alloc and get_type.
     */
    map = vnaproperty_map_alloc();
    if (map == NULL) {
	(void)printf("vnaproperty_map_alloc: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_type(map) != VNAPROPERTY_MAP) {
	(void)printf("vnaproperty_type(map) != VNAPROPERTY_MAP\n");
	result = T_FAIL;
	goto out;
    }

    /*
     * Test set.
     */
    for (int i = 0; i < N_WORDS; ++i) {
	vnaproperty_t *scalar;
	char buf[3 * sizeof(int) + 1];

	length = vnaproperty_map_count(map);
	if (length == -1) {
	    (void)printf("vnaproperty_map_count: %s (%d)\n",
		    strerror(errno), i);
	}
	if (length != i) {
	    (void)printf("vnaproperty_map_count mismatch (%d != %d)\n",
		    (int)length, i);
	    result = T_FAIL;
	    goto out;
	}
	(void)sprintf(buf, "%d", i);
	if ((scalar = vnaproperty_scalar_alloc(buf)) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_map_set(map, words[i], scalar) == -1) {
	    (void)printf("vnaproperty_map_set: %s (%d)\n",
		    strerror(errno), i);
	    result = T_FAIL;
	    goto out;
	}
    }
    length = vnaproperty_map_count(map);
    if (length != N_WORDS) {
	(void)printf("vnaproperty_map_count mismatch (%d != %d)\n",
		(int)length, (int)N_WORDS);
	result = T_FAIL;
	goto out;
    }

    /*
     * Test get.
     */
    for (int i = 0; i < N_WORDS; ++i) {
	vnaproperty_t *scalar;
	const char *value;
	char buf[3 * sizeof(int) + 1];

	(void)sprintf(buf, "%d", i);
	if ((scalar = vnaproperty_map_get(map, words[i])) == NULL) {
	    (void)printf("vnaproperty_map_get: %s (%s)\n",
		    strerror(errno), words[i]);
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
	    (void)printf("retrieved list element %d not a scalar\n", i);
	    result = T_FAIL;
	    goto out;
	}
	if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
	    (void)printf("vnaproperty_scalar_get: %s (%d)\n",
		    strerror(errno), i);
	    result = T_FAIL;
	    goto out;
	}
	if (strcmp(value, buf) != 0) {
	    (void)printf("vnaproperty_list_get miscompare \"%s\" != \"%s\"\n",
		value, buf);
	    result = T_FAIL;
	    goto out;
	}
    }
    first_scalar = vnaproperty_map_get(map, words[0]);
    assert(first_scalar != NULL);
    assert(vnaproperty_type(first_scalar) == VNAPROPERTY_SCALAR);

    /*
     * Test get of non-existent.
     */
    errno = 0;
    if (vnaproperty_map_get(map, "NotInList") != NULL || errno != ENOENT) {
	(void)printf("vnaproperty_map_get: %s (NotInList)\n",
		strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Test change via set.
     */
    for (int i = N_WORDS-1; i >= 0; --i) {
	vnaproperty_t *scalar;
	char buf[3 * sizeof(int) + 2];

	(void)sprintf(buf, "%d", -i);
	if ((scalar = vnaproperty_scalar_alloc(buf)) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_map_set(map, words[i], scalar) == -1) {
	    (void)printf("vnaproperty_map_set: %s (%d)\n",
		    strerror(errno), i);
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_map_count(map);
	if (length == -1) {
	    (void)printf("vnaproperty_map_count: %s (%d)\n",
		    strerror(errno), i);
	}
	if (length != N_WORDS) {
	    (void)printf("vnaproperty_map_count mismatch (%d != %d)\n",
		    (int)length, (int)N_WORDS);
	    result = T_FAIL;
	    goto out;
	}
    }
    for (int i = 0; i < N_WORDS; ++i) {
	vnaproperty_t *scalar;
	const char *value;
	char buf[3 * sizeof(int) + 2];

	(void)sprintf(buf, "%d", -i);
	if ((scalar = vnaproperty_map_get(map, words[i])) == NULL) {
	    (void)printf("vnaproperty_map_get: %s (%s)\n",
		    strerror(errno), words[i]);
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
	    (void)printf("retrieved list element %d not a scalar\n", i);
	    result = T_FAIL;
	    goto out;
	}
	if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
	    (void)printf("vnaproperty_scalar_get: %s (%d)\n",
		    strerror(errno), i);
	    result = T_FAIL;
	    goto out;
	}
	if (strcmp(value, buf) != 0) {
	    (void)printf("vnaproperty_list_get miscompare \"%s\" != \"%s\"\n",
		value, buf);
	    result = T_FAIL;
	    goto out;
	}
    }
    if (vnaproperty_type(first_scalar) == VNAPROPERTY_SCALAR) {
	(void)printf("first unheld scalar remained on re-set\n");
	result = T_FAIL;
	goto out;
    }
    first_scalar = vnaproperty_map_get(map, words[0]);
    assert(first_scalar != NULL);
    assert(vnaproperty_type(first_scalar) == VNAPROPERTY_SCALAR);

    /*
     * Test delete by deleting all the odd words.
     */
    for (int i = 0; i < N_WORDS / 2; ++i) {
	if (vnaproperty_map_delete(map, words[2 * i + 1]) == -1) {
	    (void)printf("vnaproperty_map_delete: %s (%s)\n",
		    strerror(errno), words[i]);
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_map_count(map);
	if (length == -1) {
	    (void)printf("vnaproperty_map_count: %s (%d)\n",
		    strerror(errno), i);
	}
	if (length != N_WORDS - i - 1) {
	    (void)printf("vnaproperty_map_count mismatch (%d != %d)\n",
		    (int)length, (int)N_WORDS - i - 1);
	    result = T_FAIL;
	    goto out;
	}
    }
    for (int i = 0; i < N_WORDS; ++i) {
	vnaproperty_t *scalar;
	const char *value;
	char buf[3 * sizeof(int) + 2];

	errno = 0;
	scalar = vnaproperty_map_get(map, words[i]);
	if (i & 1) {	/* odd should no longer be there */
	    if (scalar != NULL || errno != ENOENT) {
		(void)printf("vnaproperty_scalar_get: %s "
			"(still there %s)\n",
			strerror(errno), words[i]);
		result = T_FAIL;
		goto out;
	    }
	    continue;
	}
	if (scalar == NULL) {
	    (void)printf("vnaproperty_map_get: %s (%s)\n",
		    strerror(errno), words[i]);
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
	    (void)printf("retrieved list element %d not a scalar\n", i);
	    result = T_FAIL;
	    goto out;
	}
	if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
	    (void)printf("vnaproperty_scalar_get: %s (%d)\n",
		    strerror(errno), i);
	    result = T_FAIL;
	    goto out;
	}
	(void)sprintf(buf, "%d", -i);
	if (strcmp(value, buf) != 0) {
	    (void)printf("vnaproperty_list_get miscompare \"%s\" != \"%s\"\n",
		value, buf);
	    result = T_FAIL;
	    goto out;
	}
    }

    /*
     * Test delete of non-existent.
     */
    errno = 0;
    if (vnaproperty_map_delete(map, "NotInList") != -1 || errno != ENOENT) {
	(void)printf("vnaproperty_map_delete: %s (NotInList)\n",
		strerror(errno));
	result = T_FAIL;
	goto out;
    }
    length = vnaproperty_map_count(map);
    if (length == -1) {
	(void)printf("vnaproperty_map_count: %s (NotInList)\n",
		strerror(errno));
    }
    if (length != (N_WORDS + 1) / 2) {
	(void)printf("vnaproperty_map_count mismatch (%d != %d)\n",
		(int)length, (int)(N_WORDS + 1) / 2);
	result = T_FAIL;
	goto out;
    }

    /*
     * Test iteration.
     */
    {
	const vnaproperty_map_pair_t *vmprp;
	int count = 0;

	errno = 0;
	if ((vmprp = vnaproperty_map_begin(map)) == NULL) {
	    (void)printf("vnaproperty_map_begin: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	do {
	    if (strcmp(vmprp->vmpr_key, words[2 * count]) != 0) {
		(void)printf("iteration miscompare \"%s\" != \"%s\"\n",
		    vmprp->vmpr_key, words[2 * count]);
		result = T_FAIL;
		goto out;
	    }
	    ++count;
	    vmprp = vnaproperty_map_next(vmprp);
	} while (vmprp != NULL);
	if (count != (N_WORDS + 1) / 2) {
	    (void)printf("iteration length mismatch (%d != %d)\n",
		    (int)count, (int)(N_WORDS + 1) / 2);
	    result = T_FAIL;
	}
    }

    /*
     * Test hold and free.
     */
    {
	vnaproperty_hold(map);
	vnaproperty_free(map);
	if (vnaproperty_type(map) != VNAPROPERTY_MAP) {
	    (void)printf("held map type changed on free\n");
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(first_scalar) != VNAPROPERTY_SCALAR) {
	    (void)printf("first held list scalar element type changed "
		    "on free\n");
	    result = T_FAIL;
	    goto out;
	}
	vnaproperty_free(map);
	if (vnaproperty_type(map) == VNAPROPERTY_LIST) {
	    (void)printf("unheld map type remained on free\n");
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(first_scalar) == VNAPROPERTY_SCALAR) {
	    (void)printf("first unheld list scalar element type "
		    "remained on free\n");
	    result = T_FAIL;
	    goto out;
	}
    }
    result = T_PASS;

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
