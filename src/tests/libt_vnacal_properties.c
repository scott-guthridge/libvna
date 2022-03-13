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
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "libt.h"
#include "libt_vnacal.h"

/*
 * libt_vnacal_print_properties: print a property list
 */
void libt_vnacal_print_properties(const vnaproperty_t *vprp, int indent)
{
    if (vprp == NULL) {
	for (int i = 0; i < indent; ++i) {
	    (void)printf("    ");
	}
	printf(".\n");
	return;
    }
    switch (vnaproperty_type(vprp, ".")) {
    case 's':
	for (int i = 0; i < indent; ++i) {
	    (void)printf("    ");
	}
	(void)printf("\"%s\"\n", vnaproperty_get(vprp, "."));
	return;

    case 'm':
	{
	    const char **keys;

	    keys = vnaproperty_keys(vprp, "{}");
	    assert(keys != NULL);
	    for (const char **cpp = keys; *cpp != NULL; ++cpp) {
		char *key = NULL;
		vnaproperty_t *subtree;

		if ((key = vnaproperty_quote_key(*cpp)) == NULL) {
		    (void)printf("vnaproperty_quote_keys: %s\n",
			    strerror(errno));
		    continue;
		}
		subtree = vnaproperty_get_subtree(vprp, "%s", key);
		assert(subtree != NULL);
		for (int i = 0; i < indent; ++i) {
		    (void)printf("    ");
		}
		(void)printf(".%s\n", key);
		libt_vnacal_print_properties(subtree, indent + 1);
		free((void *)key);
	    }
	    free((void *)keys);
	}
	return;

    case 'l':
	{
	    int count;

	    count = vnaproperty_count(vprp, "[]");
	    assert(count >= 0);
	    for (int i = 0; i < count; ++i) {
		vnaproperty_t *subtree;

		subtree = vnaproperty_get_subtree(vprp, "[%d]", i);
		assert(subtree != NULL);
		for (int i = 0; i < indent; ++i) {
		    (void)printf("    ");
		}
		(void)printf("[%d]\n", i);
		libt_vnacal_print_properties(subtree, indent + 1);
	    }
	}
	return;

    default:
	abort();
    }
}
