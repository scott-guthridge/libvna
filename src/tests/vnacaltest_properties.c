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
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include "test.h"
#include "vnacaltest.h"

/*
 * test_vnacal_print_properties: print a property list
 */
void test_vnacal_print_properties(const vnaproperty_t *vprp, int indent)
{
    if (vprp == NULL) {
	for (int i = 0; i < indent; ++i) {
	    (void)printf("    ");
	}
	printf(".\n");
	return;
    }
    switch (vnaproperty_type(vprp)) {
    case VNAPROPERTY_SCALAR:
	for (int i = 0; i < indent; ++i) {
	    (void)printf("    ");
	}
	(void)printf("\"%s\"\n", vnaproperty_scalar_get(vprp));
	return;

    case VNAPROPERTY_MAP:
	{
	    const vnaproperty_map_pair_t *vmprp;

	    for (vmprp = vnaproperty_map_begin(vprp); vmprp != NULL;
		    vmprp = vnaproperty_map_next(vmprp)) {
		for (int i = 0; i < indent; ++i) {
		    (void)printf("    ");
		}
		(void)printf(".%s\n", vmprp->vmpr_key);
		test_vnacal_print_properties(vmprp->vmpr_value, indent + 1);
	    }
	}
	return;

    case VNAPROPERTY_LIST:
	{
	    int count = vnaproperty_list_count(vprp);

	    for (int i = 0; i < count; ++i) {
		for (int i = 0; i < indent; ++i) {
		    (void)printf("    ");
		}
		(void)printf("[%d]\n", i);
		test_vnacal_print_properties(vnaproperty_list_get(vprp, i),
			indent + 1);
	    }
	}
	return;

    default:
	abort();
    }
}
