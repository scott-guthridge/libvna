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
 * test_vnaproperty_list
 */
static test_result_t test_vnaproperty_list()
{
    vnaproperty_t *list;
    vnaproperty_t *first_scalar = NULL;
    int length;
    test_result_t result = T_SKIPPED;

    /*
     * Test alloc and get_type.
     */
    list = vnaproperty_list_alloc();
    if (list == NULL) {
	(void)printf("vnaproperty_list_alloc: %s\n", strerror(errno));
	result = T_FAIL;
	goto out;
    }
    if (vnaproperty_type(list) != VNAPROPERTY_LIST) {
	(void)printf("vnaproperty_type(list) != VNAPROPERTY_LIST\n");
	result = T_FAIL;
	goto out;
    }

    /*
     * Test append.
     */
    for (int i = 0; i < 100; ++i) {
	char buf[3 * sizeof(int) + 1];
	vnaproperty_t *scalar;

	length = vnaproperty_list_count(list);
	if (length == -1) {
	    (void)printf("vnaproperty_list_count: %s (%d)\n",
		    strerror(errno), i);
	}
	if (length != i) {
	    (void)printf("vnaproperty_list_count mismatch (%d != %d)\n",
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
	if (vnaproperty_list_append(list, scalar) == -1) {
	    (void)printf("vnaproperty_list_append: %s (%d)\n",
		    strerror(errno), i);
	    result = T_FAIL;
	    goto out;
	}
    }
    length = vnaproperty_list_count(list);
    if (length != 100) {
	(void)printf("vnaproperty_list_count mismatch (%d != 100)\n",
		(int)length);
	result = T_FAIL;
	goto out;
    }

    /*
     * Test get.
     */
    for (int i = 0; i < 100; ++i) {
	vnaproperty_t *scalar;
	const char *value;
	char buf[3 * sizeof(int) + 1];

	(void)sprintf(buf, "%d", i);
	if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
	    (void)printf("vnaproperty_list_get: %s (%d)\n",
		    strerror(errno), i);
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
    first_scalar = vnaproperty_list_get(list, 0);

    /*
     * Test set.
     */
    {
	vnaproperty_t *scalar;

	/*
	 * Test bounds check.
	 */
	if ((scalar = vnaproperty_scalar_alloc("out-of-bounds")) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s (fifty)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	errno = 0;
	if (vnaproperty_list_set(list, -1, scalar) != -1 || errno != EINVAL) {
	    (void)printf("vnaproperty_list_set: %s (bounds -1)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	vnaproperty_free(scalar);

	/*
	 * Test set in middle.
	 *   starting state: 0..99
	 *   ending state:   0..49 "fifty" 51..99
	 */
	if ((scalar = vnaproperty_scalar_alloc("fifty")) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s (fifty)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_list_set(list, 50, scalar) == -1) {
	    (void)printf("vnaproperty_list_set: %s (fifty)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_list_count(list);
	if (length != 100) {
	    (void)printf("vnaproperty_list_count mismatch "
		    "(set %d != 100)\n", (int)length);
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < 100; ++i) {
	    vnaproperty_t *scalar;
	    const char *value;
	    char buf[64];

	    (void)sprintf(buf, "%d", i);
	    if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
		(void)printf("vnaproperty_list_get: %s (set %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
		(void)printf("retrieved list element %d not a scalar (set)\n",
			i);
		result = T_FAIL;
		goto out;
	    }
	    if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
		(void)printf("vnaproperty_scalar_get: %s (set %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 50) {
		(void)strcpy(buf, "fifty");
	    }
	    if (strcmp(value, buf) != 0) {
		(void)printf("vnaproperty_list_get miscompare "
		    "\"%s\" != \"%s\" (set)\n", value, buf);
		result = T_FAIL;
		goto out;
	    }
	}

	/*
	 * Test setting past the end.
	 *   starting state: 0..49 "fifty" 51..99
	 *   ending state:   0..49 "fifty" 51..99 ~ ~ "hundred-two"
	 */
	if ((scalar = vnaproperty_scalar_alloc("hundred-two")) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s (fifty)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_list_set(list, 102, scalar) == -1) {
	    (void)printf("vnaproperty_list_set: %s (fifty)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_list_count(list);
	if (length != 103) {
	    (void)printf("vnaproperty_list_count mismatch "
		    "(set %d != 103)\n", (int)length);
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < 103; ++i) {
	    vnaproperty_t *scalar;
	    const char *value;
	    char buf[64];

	    if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
		(void)printf("vnaproperty_list_get: %s (set %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
		(void)printf("retrieved list element %d not a scalar (set)\n",
			i);
		result = T_FAIL;
		goto out;
	    }
	    if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
		(void)printf("vnaproperty_scalar_get: %s (set %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (i == 50) {
		(void)strcpy(buf, "fifty");
	    } else if (i == 100 || i == 101) {
		(void)strcpy(buf, "~");
	    } else if (i == 102) {
		(void)strcpy(buf, "hundred-two");
	    } else {
		(void)sprintf(buf, "%d", i);
	    }
	    if (strcmp(value, buf) != 0) {
		(void)printf("vnaproperty_list_get miscompare "
		    "\"%s\" != \"%s\" (set)\n", value, buf);
		result = T_FAIL;
		goto out;
	    }
	}
    }

    /*
     * Test insert.
     *	 starting state: 0..49 "fifty" 51..99 ~ ~ "hundred-two"
     */
    {
	vnaproperty_t *scalar;

	/*
	 * Test bounds check.
	 */
	if ((scalar = vnaproperty_scalar_alloc("50")) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s (END)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	errno = 0;
	if (vnaproperty_list_insert(list, -1, scalar) != -1 ||
		errno != EINVAL) {
	    (void)printf("vnaproperty_list_insert: %s (bounds -1)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}

	/*
	 * Test insert in the middle.
	 *   ending state: 0..50 [51]="fifty" [52]=51..[100]=99
	 *	[101]=~ [102]=~ [103]="hundred-two"
	 */
	if (vnaproperty_list_insert(list, 50, scalar) == -1) {
	    (void)printf("vnaproperty_list_insert: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_list_count(list);
	if (length != 104) {
	    (void)printf("vnaproperty_list_count mismatch "
		    "(insert %d != 104)\n", (int)length);
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < 102; ++i) {
	    vnaproperty_t *scalar;
	    const char *value;
	    char buf[64];

	    if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
		(void)printf("vnaproperty_list_get: %s (insert %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
		(void)printf("retrieved list element %d not a scalar "
			"(insert)\n", i);
		result = T_FAIL;
		goto out;
	    }
	    if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
		(void)printf("vnaproperty_scalar_get: %s (insert %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (i <= 50) {
		(void)sprintf(buf, "%d", i);
	    } else if (i == 51) {
		(void)strcpy(buf, "fifty");
	    } else if (i == 101 || i == 102) {
		(void)strcpy(buf, "~");
	    } else if (i == 103) {
		(void)strcpy(buf, "hundred-two");
	    } else {	/* 52 <= i <= 100 */
		(void)sprintf(buf, "%d", i - 1);
	    }
	    if (strcmp(value, buf) != 0) {
		(void)printf("vnaproperty_list_get miscompare "
		    "\"%s\" != \"%s\" (insert)\n", value, buf);
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
	if ((scalar = vnaproperty_scalar_alloc("one-o-four")) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s (END)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_list_insert(list, 104, scalar) == -1) {
	    (void)printf("vnaproperty_list_insert: %s\n", strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_list_count(list);
	if (length != 105) {
	    (void)printf("vnaproperty_list_count mismatch "
		    "(insert %d != 105)\n", (int)length);
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < 105; ++i) {
	    vnaproperty_t *scalar;
	    const char *value;
	    char buf[64];

	    if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
		(void)printf("vnaproperty_list_get: %s (insert %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
		(void)printf("retrieved list element %d not a scalar "
			"(insert)\n", i);
		result = T_FAIL;
		goto out;
	    }
	    if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
		(void)printf("vnaproperty_scalar_get: %s (insert %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (i <= 50) {
		(void)sprintf(buf, "%d", i);
	    } else if (i == 51) {
		(void)strcpy(buf, "fifty");
	    } else if (i == 101 || i == 102) {
		(void)strcpy(buf, "~");
	    } else if (i == 103) {
		(void)strcpy(buf, "hundred-two");
	    } else if (i == 104) {
		(void)strcpy(buf, "one-o-four");
	    } else {	/* 52 <= i <= 100 */
		(void)sprintf(buf, "%d", i - 1);
	    }
	    if (strcmp(value, buf) != 0) {
		(void)printf("vnaproperty_list_get miscompare "
		    "\"%s\" != \"%s\" (insert)\n", value, buf);
		result = T_FAIL;
		goto out;
	    }
	}
    }

    /*
     * Test delete.
     *   starting state: 0..50 [51]="fifty" [52]=51..[100]=99
     *	   [101]=~ [102]=~ [103]="hundred-two" [104]="one-o-four"
     */
    {
	/*
	 * Test bounds check.
	 */
	errno = 0;
	if (vnaproperty_list_delete(list, -1) != -1 || errno != EINVAL) {
	    (void)printf("vnaproperty_list_delete: %s (bounds -1)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	errno = 0;
	if (vnaproperty_list_delete(list, 105) != -1 || errno != ENOENT) {
	    (void)printf("vnaproperty_list_delete: %s (bounds 103)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}

	/*
	 * Test delete in the middle.
	 *   ending state:   0..99 [100]=~ [101]=~ [102]="hundred-two"
	 *     [103]="one-o-four"
	 */
	if (vnaproperty_list_delete(list, 51) == -1) {
	    (void)printf("vnaproperty_list_delete: %s (delete 51)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_list_count(list);
	if (length != 104) {
	    (void)printf("vnaproperty_list_count mismatch "
		    "(delete %d != 104)\n", (int)length);
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < 102; ++i) {
	    vnaproperty_t *scalar;
	    const char *value;
	    char buf[64];

	    if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
		(void)printf("vnaproperty_list_get: %s (delete %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
		(void)printf("retrieved list element %d not a scalar "
			"(delete)\n", i);
		result = T_FAIL;
		goto out;
	    }
	    if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
		(void)printf("vnaproperty_scalar_get: %s (delete %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (i <= 99) {
		(void)sprintf(buf, "%d", i);
	    } else if (i == 100 || i == 101) {
		(void)strcpy(buf, "~");
	    } else if (i == 102) {
		(void)strcpy(buf, "hundred-two");
	    } else if (i == 103) {
		(void)strcpy(buf, "one-o-four");
	    } else {
		abort();
	    }
	    if (strcmp(value, buf) != 0) {
		(void)printf("vnaproperty_list_get miscompare "
		    "\"%s\" != \"%s\" (delete 51)\n", value, buf);
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
	if (vnaproperty_list_delete(list, 103) == -1) {
	    (void)printf("vnaproperty_list_delete: %s (delete 103)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	length = vnaproperty_list_count(list);
	if (length != 103) {
	    (void)printf("vnaproperty_list_count mismatch "
		    "(delete %d != 101)\n", (int)length);
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < 103; ++i) {
	    vnaproperty_t *scalar;
	    const char *value;
	    char buf[64];

	    if ((scalar = vnaproperty_list_get(list, i)) == NULL) {
		(void)printf("vnaproperty_list_get: %s (delete %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
		(void)printf("retrieved list element %d not a scalar "
			"(delete)\n", i);
		result = T_FAIL;
		goto out;
	    }
	    if ((value = vnaproperty_scalar_get(scalar)) == NULL) {
		(void)printf("vnaproperty_scalar_get: %s (delete %d)\n",
			strerror(errno), i);
		result = T_FAIL;
		goto out;
	    }
	    if (i <= 99) {
		(void)sprintf(buf, "%d", i);
	    } else if (i == 100 || i == 101) {
		(void)strcpy(buf, "~");
	    } else if (i == 102) {
		(void)strcpy(buf, "hundred-two");
	    } else {
		abort();
	    }
	    if (strcmp(value, buf) != 0) {
		(void)printf("vnaproperty_list_get miscompare "
		    "\"%s\" != \"%s\" (delete 101)\n", value, buf);
		result = T_FAIL;
		goto out;
	    }
	}
    }

    /*
     * Test hold and free.
     */
    {
	vnaproperty_t *scalar;

	if ((scalar = vnaproperty_scalar_alloc("END")) == NULL) {
	    (void)printf("vnaproperty_scalar_alloc: %s (END)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_list_append(list, scalar) == -1) {
	    (void)printf("vnaproperty_list_append: %s (END)\n",
		    strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	vnaproperty_hold(list);
	vnaproperty_free(list);
	if (vnaproperty_type(list) != VNAPROPERTY_LIST) {
	    (void)printf("held list type changed on free\n");
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(scalar) != VNAPROPERTY_SCALAR) {
	    (void)printf("last held list scalar element type changed "
		    "on free\n");
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(first_scalar) != VNAPROPERTY_SCALAR) {
	    (void)printf("first held list scalar element type changed "
		    "on free\n");
	    result = T_FAIL;
	    goto out;
	}
	vnaproperty_free(list);
	if (vnaproperty_type(list) == VNAPROPERTY_LIST) {
	    (void)printf("unheld list type remained on free\n");
	    result = T_FAIL;
	    goto out;
	}
	if (vnaproperty_type(scalar) == VNAPROPERTY_SCALAR) {
	    (void)printf("last unheld list scalar element type "
		    "remained on free\n");
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
    exit(test_vnaproperty_list());
}
