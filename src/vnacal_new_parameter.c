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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "archdep.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/**********************************************************************
 * Parameter Hash Table
 **********************************************************************/

#define INITIAL_HASH_SIZE	8

/*
 * hash_expand: allocate/resize the hash table
 *   @vnphp: new calibration parameter collection
 */
static int hash_expand(vnacal_new_parameter_hash_t *vnphp)
{
    int old_allocation = vnphp->vnph_allocation;
    int new_allocation = MAX(2 * old_allocation, INITIAL_HASH_SIZE);
    vnacal_new_parameter_t **new_table = NULL;

    /*
     * Resize the table, initializing new buckets to NULL.
     */
    if ((new_table = realloc(vnphp->vnph_table, new_allocation *
		    sizeof(vnacal_new_parameter_t *))) == NULL) {
	return -1;
    }
    (void)memset((void *)&new_table[old_allocation], 0,
	(new_allocation - old_allocation) * sizeof(vnacal_new_parameter_t *));
    vnphp->vnph_table = new_table;
    vnphp->vnph_allocation = new_allocation;

    /*
     * Walk the old hash chains and rehash all elements.
     */
    for (int chain = 0; chain < old_allocation; ++chain) {
	vnacal_new_parameter_t *head, *cur;

	head = new_table[chain];
	new_table[chain] = NULL;
	while ((cur = head) != NULL) {
	    int hashval = cur->vnpr_parameter->vpmr_index;
	    vnacal_new_parameter_t *next, **anchor;

	    head = cur->vnpr_hash_next;
	    anchor = &new_table[hashval % new_allocation];
	    for (; (next = *anchor) != NULL; anchor = &next->vnpr_hash_next) {
		if (next->vnpr_parameter->vpmr_index > hashval) {
		    break;
		}
	    }
	    cur->vnpr_hash_next = next;
	    *anchor = cur;
	}
    }
    return 0;
}

/*
 * hash_lookup: find a parameter in the hash table
 *   @vnphp: parameter hash table
 *   @parameter: parameter index
 */
static vnacal_new_parameter_t *hash_lookup(
	const vnacal_new_parameter_hash_t *vnphp, int parameter)
{
    vnacal_new_parameter_t *vnprp;

    vnprp = vnphp->vnph_table[parameter % vnphp->vnph_allocation];
    for (; vnprp != NULL; vnprp = vnprp->vnpr_hash_next) {
	int index = VNACAL_GET_PARAMETER_INDEX(vnprp->vnpr_parameter);

	if (index == parameter) {
	    return vnprp;
	}
	if (index >  parameter) {
	    break;
	}
    }
    return NULL;
}

/*
 * hash_insert: insert a parameter into the hash table
 *   @vnphp: parameter hash table
 *   @vnprp: node to insert
 */
static void hash_insert(vnacal_new_parameter_hash_t *vnphp,
	vnacal_new_parameter_t *vnprp)
{
    int parameter = VNACAL_GET_PARAMETER_INDEX(vnprp->vnpr_parameter);
    vnacal_new_parameter_t **anchor, *next;

    anchor = &vnphp->vnph_table[parameter % vnphp->vnph_allocation];
    for (; (next = *anchor) != NULL; anchor = &next->vnpr_hash_next) {
	if (VNACAL_GET_PARAMETER_INDEX(next->vnpr_parameter) > parameter)
	    break;
    }
    vnprp->vnpr_hash_next = next;
    *anchor = vnprp;
    if (++vnphp->vnph_count >= vnphp->vnph_allocation) {
	(void)hash_expand(vnphp);
    }
}

/*
 * _vnacal_new_init_parameter_hash: set up the parameter hash
 *   @function: name of user-called function
 *   @vnphp: hash table
 *
 * Caller must log errors.
 */
int _vnacal_new_init_parameter_hash(const char *function,
	vnacal_new_parameter_hash_t *vnphp)
{
    (void)memset((void *)vnphp, 0, sizeof(*vnphp));
    return hash_expand(vnphp);
}

/*
 * _vnacal_new_free_parameter_hash: free the parameter hash
 *   @vnphp: hash table
 */
void _vnacal_new_free_parameter_hash(vnacal_new_parameter_hash_t *vnphp)
{
    if (vnphp->vnph_table != NULL) {
	for (int bucket = 0; bucket < vnphp->vnph_allocation; ++bucket) {
	    vnacal_new_parameter_t *vnprp;

	    while ((vnprp = vnphp->vnph_table[bucket]) != NULL) {
		vnphp->vnph_table[bucket] = vnprp->vnpr_hash_next;

		_vnacal_release_parameter(vnprp->vnpr_parameter);
		free((void *)vnprp);
	    }
	}
	free((void *)vnphp->vnph_table);
	(void)memset((void *)vnphp, 0, sizeof(*vnphp));
    }
}

/*
 * check_single_frequency_range: check frequency range of new parameter
 *   @function: name of user-called function
 *   @vnp: pointer to vnacal_new_t structure
 *   @vprmp: new parameter to check
 */
static int check_single_frequency_range(const char *function,
	vnacal_new_t *vnp, double fmin, double fmax,
	vnacal_parameter_t *vpmrp)
{
    vnacal_t *vcp = vnp->vn_vcp;
    double pfmin, pfmax;
    double lower, upper;

    _vnacal_get_parameter_frange(vpmrp, &pfmin, &pfmax);
    lower = (1.0 + VNACAL_F_EXTRAPOLATION) * fmin;
    upper = (1.0 - VNACAL_F_EXTRAPOLATION) * fmin;
    if (pfmin > lower || pfmax < upper) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: frequency range %.3e..%.3e of "
		"parameter %d is outside of calibration range %.3e..%.3e",
		function, pfmin, pfmax, VNACAL_GET_PARAMETER_INDEX(vpmrp),
		fmin, fmax);
	return -1;
    }
    return 0;
}

/*
 * _vnacal_new_check_all_frequency_ranges: check f range of all parameters
 *   @function: name of user-called function
 *   @vnp: pointer to vnacal_new_t structure
 *
 *   Used to check existing parameters in vnacal_new_set_frequency_vector.
 */
int _vnacal_new_check_all_frequency_ranges(const char *function,
	vnacal_new_t *vnp, double fmin, double fmax)
{
    vnacal_new_parameter_hash_t *vnphp = &vnp->vn_parameter_hash;

    for (int bucket = 0; bucket < vnphp->vnph_allocation; ++bucket) {
	vnacal_new_parameter_t *vnprp;

	vnprp = vnphp->vnph_table[bucket];
	for (; vnprp != NULL; vnprp = vnprp->vnpr_hash_next) {
	    if (check_single_frequency_range(function, vnp,
			fmin, fmax, vnprp->vnpr_parameter) == -1) {
		return -1;
	    }
	}
    }
    return 0;
}


/*
 * _vnacal_new_get_parameter: add/find parameter
 *   @function: name of user-called function
 *   @vnp: pointer to vnacal_new_t structure
 *   @parameter: parameter index such as VNACAL_ZERO
 */
vnacal_new_parameter_t *_vnacal_new_get_parameter(const char *function,
	vnacal_new_t *vnp, int parameter)
{
    vnacal_t *vcp = vnp->vn_vcp;
    vnacal_new_parameter_hash_t *vnphp = &vnp->vn_parameter_hash;
    vnacal_new_parameter_t *vnprp;
    vnacal_parameter_t *vpmrp;
    vnacal_parameter_type_t type;
    vnacal_new_parameter_t *ncprp_correlate = NULL;

    /*
     * Search for the parameter in the hash and return if found.
     */
    if ((vnprp = hash_lookup(vnphp, parameter)) != NULL) {
	return vnprp;
    }

    /*
     * Look-up the parameter in the vnacal_t structure.  If not found,
     * parameter is deleted or invalid.
     */
    if ((vpmrp = _vnacal_get_parameter(vcp, parameter)) == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: invalid parameter index %d",
		function, parameter);
	return NULL;
    }
    type = VNACAL_GET_PARAMETER_TYPE(vpmrp);

    /*
     * If the frequency vector has been given, check the frequency range.
     */
    if (vnp->vn_frequencies_valid) {
	if (check_single_frequency_range(function, vnp,
		    vnp->vn_frequency_vector[0],
		    vnp->vn_frequency_vector[vnp->vn_frequencies - 1],
		    vpmrp) == -1) {
	    return NULL;
	}
    }

    /*
     * If we're given a correlated parameter, recurse to get the correlate.
     */
    if (type == VNACAL_CORRELATED) {
	vnacal_parameter_t *vpmrp_correlate = VNACAL_GET_PARAMETER_OTHER(vpmrp);

	if ((ncprp_correlate = _vnacal_new_get_parameter(function, vnp,
			VNACAL_GET_PARAMETER_INDEX(vpmrp_correlate))) == NULL) {
	    free((void *)vnprp);
	    return NULL;
	}
    }

    /*
     * Create a new vnacal_new_parameter_t structure and add to hash table.
     */
    if ((vnprp = malloc(sizeof(vnacal_new_parameter_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)vnprp, 0, sizeof(*vnprp));
    _vnacal_hold_parameter(vpmrp);
    vnprp->vnpr_parameter = vpmrp;
    vnprp->vnpr_cmp = vnp;
    hash_insert(vnphp, vnprp);

    /*
     * If unknown, assign an index and add to the unknown parameter list.
     */
    if (type == VNACAL_UNKNOWN || type == VNACAL_CORRELATED) {
	vnprp->vnpr_unknown = true;
	vnprp->vnpr_unknown_index = vnp->vn_unknown_parameters++;
	vnprp->vnpr_correlate = ncprp_correlate;
	*vnp->vn_unknown_parameter_anchor = vnprp;
	vnp->vn_unknown_parameter_anchor = &vnprp->vnpr_next_unknown;
    }

    return vnprp;
}
