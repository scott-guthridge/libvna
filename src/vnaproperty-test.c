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
#include <unistd.h>
#include "vnaproperty_internal.h"


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
 * test_result_type
 */
typedef enum test_result {
    T_PASS,
    T_FAIL,
    T_SKIPPED
} test_result_type;

/*
 * test counters
 */
static int test_count = 0;
static int fail_count = 0;

/*
 * report_test_result: report a test result
 */
static void report_test_result(const char *test_name, test_result_type result)
{
    const char *result_name;

    switch (result) {
    case T_PASS:
	result_name = "PASS";
	break;
    case T_FAIL:
	result_name = "FAIL";
	break;
    case T_SKIPPED:
	result_name = "SKIPPED";
	break;
    default:
	abort();
	/*NOTREACHED*/
    }
    (void)printf("Test %2d: %-58s %s\n", ++test_count, test_name, result_name);
    (void)fflush(stdout);
    if (result == T_FAIL) {
	++fail_count;
    }
}

/*
 * test_vnaproperty_scalar
 */
void test_vnaproperty_scalar()
{
    vnaproperty_t *scalar;
    const char *value;
    static const char text1[] = "abcdefghijklmnopqrstuvwxyz";
    static const char text2[] = "0123456789";
    char buf[64];
    test_result_type result = T_SKIPPED;

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
    report_test_result("Scalar", result);
}

/*
 * test_vnaproperty_list
 */
void test_vnaproperty_list()
{
    vnaproperty_t *list;
    vnaproperty_t *first_scalar = NULL;
    int length;
    test_result_type result = T_SKIPPED;

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
    report_test_result("List", result);
}

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
void test_vnaproperty_map()
{
    vnaproperty_t *map;
    int length;
    vnaproperty_t *first_scalar = NULL;
    test_result_type result = T_SKIPPED;

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
    report_test_result("Map", result);
}

/*
 * test_vnaproperty_expr
 */
void test_vnaproperty_expr()
{
    vnaproperty_t *root = NULL;
    vnaproperty_type_t type;
    int count;
    const char *value;
    const char **keys;
    test_result_type result = T_SKIPPED;

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
    report_test_result("Expr", result);
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
    /*NOTREACHED*/
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
	    /*NOTREACHED*/
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
	/*NOTREACHED*/
    }

    test_vnaproperty_scalar();
    test_vnaproperty_list();
    test_vnaproperty_map();
    test_vnaproperty_expr();

    exit(fail_count != 0);
    /*NOTREACHED*/
}
