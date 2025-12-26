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
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/*
 * _vnacal_get_parameter: return a pointer to the parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 *
 * Returns NULL if not found.  Caller must report the error.
 */
vnacal_parameter_t *_vnacal_get_parameter(const vnacal_t *vcp, int parameter)
{
    const vnacal_parameter_collection_t *vprmcp = &vcp->vc_parameter_collection;
    vnacal_parameter_t *vpmrp;

    if (parameter < 0 || parameter >= vprmcp->vprmc_allocation ||
	    (vpmrp = vprmcp->vprmc_vector[parameter]) == NULL ||
	    vpmrp->vpmr_deleted) {
	return NULL;
    }
    assert(vpmrp->vpmr_index == parameter);
    assert(vpmrp->vpmr_vcp == vcp);
    return vpmrp;
}

/*
 * _vnacal_get_parameter_frange: get frequency limits for the given parameter
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 *   @fmin: address of double to receive minimum
 *   @fmax: address of double to receive maximum
 */
void _vnacal_get_parameter_frange(vnacal_parameter_t *vpmrp,
	double *fmin, double *fmax)
{
    vnacal_parameter_t *vpmrp_orig = vpmrp;

    for (;;) {
	switch (vpmrp->vpmr_type) {
	case VNACAL_SCALAR:
	    *fmin = 0.0;
	    *fmax = INFINITY;
	    break;

	case VNACAL_VECTOR:
	    *fmin = vpmrp->vpmr_frequency_vector[0];
	    *fmax = vpmrp->vpmr_frequency_vector[vpmrp->vpmr_frequencies - 1];
	    break;

	case VNACAL_UNKNOWN:
	case VNACAL_CORRELATED:
	    vpmrp = vpmrp->vpmr_other;
	    continue;

	default:
	    assert(!"unexpected parameter type");
	}
	break;
    }

    /*
     * If the original object is of type VNA_CORRELATED, then further
     * restrict the range based on the sigma frequencies.
     */
    if (vpmrp_orig->vpmr_type == VNACAL_CORRELATED &&
	    vpmrp_orig->vpmr_sigma_frequency_vector != NULL) {
	int sf = vpmrp_orig->vpmr_sigma_frequencies;
	double smin = vpmrp_orig->vpmr_sigma_frequency_vector[0];
	double smax = vpmrp_orig->vpmr_sigma_frequency_vector[sf - 1];

	if (smin > *fmin) {
	    *fmin = smin;
	}
	if (smax < *fmax) {
	    *fmax = smax;
	}
    }
}

/*
 * _vnacal_get_parameter_value_i: get the value of the parameter at f
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 *   @frequency: frequency at which to evaluate the value
 */
double complex _vnacal_get_parameter_value_i(vnacal_parameter_t *vpmrp,
	double frequency)
{
    for (;;) {
	switch (vpmrp->vpmr_type) {
	case VNACAL_NEW:
	    break;

	case VNACAL_SCALAR:
	    return vpmrp->vpmr_coefficient;

	case VNACAL_VECTOR:
	    return _vnacal_rfi(vpmrp->vpmr_frequency_vector,
		    vpmrp->vpmr_coefficient_vector,
		    vpmrp->vpmr_frequencies,
		    MIN(vpmrp->vpmr_frequencies, VNACAL_MAX_M),
		    &vpmrp->vpmr_segment,
		    frequency);

	case VNACAL_UNKNOWN:
	case VNACAL_CORRELATED:
	    vpmrp = vpmrp->vpmr_other;
	    continue;
	}
	break;
    }
    assert(!"unexpected parameter type");
}

/*
 * _vnacal_alloc_parameter: allocate a vnacal_parameter and return index
 *   @function: name of user-called function
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
vnacal_parameter_t *_vnacal_alloc_parameter(const char *function, vnacal_t *vcp)
{
    vnacal_parameter_collection_t *vprmcp = &vcp->vc_parameter_collection;
    vnacal_parameter_t *vpmrp = NULL;
    int parameter;

    /*
     * Find a free slot in the table, extending the table if necessary.
     */
    if (vprmcp->vprmc_count < vprmcp->vprmc_allocation) {
	parameter = vprmcp->vprmc_first_free;
	while (vprmcp->vprmc_vector[parameter] != NULL) {
	    ++parameter;
	    assert(parameter < vprmcp->vprmc_allocation);
	}
	vprmcp->vprmc_first_free = parameter + 1;

    } else {
	vnacal_parameter_t **vpmrpp;
	int old_allocation = vprmcp->vprmc_allocation;
	int new_allocation;

	if (old_allocation < 3) {
	    new_allocation = 3;
	} else if (old_allocation < 8) {
	    new_allocation = 8;
	} else {
	    new_allocation = 2 * old_allocation;
	}
	if ((vpmrpp = realloc(vprmcp->vprmc_vector, new_allocation *
			sizeof(vnacal_parameter_t *))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "realloc: %s", strerror(errno));
	    return NULL;
	}
	(void)memset((void *)&vpmrpp[old_allocation], 0,
		(new_allocation - old_allocation) *
		sizeof(vnacal_parameter_t *));
	vprmcp->vprmc_vector = vpmrpp;
	vprmcp->vprmc_allocation = new_allocation;
	parameter = vprmcp->vprmc_count;
    }

    /*
     * Allocate and init the new parameter.  Add it to the table.
     */
    if ((vpmrp = malloc(sizeof(vnacal_parameter_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)vpmrp, 0, sizeof(*vpmrp));
    vpmrp->vpmr_type = VNACAL_NEW;
    vpmrp->vpmr_deleted = false;
    vpmrp->vpmr_hold_count = 1;
    vpmrp->vpmr_index = parameter;
    vpmrp->vpmr_segment = 0;
    vpmrp->vpmr_vcp = vcp;
    vprmcp->vprmc_vector[parameter] = vpmrp;
    ++vprmcp->vprmc_count;
    return vpmrp;
}

/*
 * _vnacal_free_parameter: remove a parameter from the table and free
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 */
static void _vnacal_free_parameter(vnacal_parameter_t *vpmrp)
{
    vnacal_t *vcp = vpmrp->vpmr_vcp;
    vnacal_parameter_collection_t *vprmcp = &vcp->vc_parameter_collection;
    int parameter = vpmrp->vpmr_index;

    assert(vpmrp != NULL);
    assert(vprmcp->vprmc_count >= 1);
    vprmcp->vprmc_vector[parameter] = NULL;
    --vprmcp->vprmc_count;
    if (parameter < vprmcp->vprmc_first_free) {
	vprmcp->vprmc_first_free = parameter;
    }
    switch (vpmrp->vpmr_type) {
    case VNACAL_CORRELATED:
	if (vpmrp->vpmr_sigma_frequency_vector !=
		vpmrp->vpmr_other->vpmr_frequency_vector) {
	    free((void *)vpmrp->vpmr_sigma_frequency_vector);
	}
	free((void *)vpmrp->vpmr_sigma_vector);
	free((void *)vpmrp->vpmr_sigma_spline);
	/*FALLTHROUGH*/

    case VNACAL_UNKNOWN:
	if (vpmrp->vpmr_other != NULL) {
	    _vnacal_release_parameter(vpmrp->vpmr_other);
	}
	/*FALLTHROUGH*/

    case VNACAL_VECTOR:
	free((void *)vpmrp->vpmr_frequency_vector);
	free((void *)vpmrp->vpmr_coefficient_vector);
	break;

    default:
	break;
    }
    free((void *)vpmrp);
}

/*
 * _vnacal_hold_parameter: increase the hold count on a parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 */
void _vnacal_hold_parameter(vnacal_parameter_t *vpmrp)
{
    ++vpmrp->vpmr_hold_count;
}

/*
 * _vnacal_release_parameter: decrease the hold count and free if zero
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 */
void _vnacal_release_parameter(vnacal_parameter_t *vpmrp)
{
    assert(vpmrp->vpmr_hold_count > 0);
    if (--vpmrp->vpmr_hold_count == 0) {
	assert(vpmrp->vpmr_deleted);
	_vnacal_free_parameter(vpmrp);
    }
}

/*
 * _vnacal_setup_parameter_collection: allocate the parameter collection
 *   @function: name of user-called function
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
int _vnacal_setup_parameter_collection(const char *function, vnacal_t *vcp)
{
    vnacal_parameter_t *vpmrp;

    /*
     * Zero the parameter collection structure.
     */
    (void)memset((void *)&vcp->vc_parameter_collection, 0,
	    sizeof(vcp->vc_parameter_collection));

    /*
     * Create the match parameter.
     */
    vpmrp = _vnacal_alloc_parameter(function, vcp);
    if (vpmrp == NULL) {
	goto error;
    }
    vpmrp->vpmr_type = VNACAL_SCALAR;
    vpmrp->vpmr_coefficient = 0.0;
    assert(vpmrp->vpmr_index == VNACAL_MATCH);

    /*
     * Create the open parameter.
     */
    vpmrp = _vnacal_alloc_parameter(function, vcp);
    if (vpmrp == NULL) {
	goto error;
    }
    vpmrp->vpmr_type = VNACAL_SCALAR;
    vpmrp->vpmr_coefficient = 1.0;
    assert(vpmrp->vpmr_index == VNACAL_OPEN);

    /*
     * Create the short parameter.
     */
    vpmrp = _vnacal_alloc_parameter(function, vcp);
    if (vpmrp == NULL) {
	goto error;
    }
    vpmrp->vpmr_type = VNACAL_SCALAR;
    vpmrp->vpmr_coefficient = -1.0;
    assert(vpmrp->vpmr_index == VNACAL_SHORT);

#if VNACAL_PREDEFINED_PARAMETERS != 3
#error "missing initializations in _vnacal_setup_parameter_collection"
#endif
    return 0;

error:
    _vnacal_teardown_parameter_collection(vcp);
    return -1;
}

/*
 * _vnacal_teardown_parameter_collection: free the parameter collection
 */
void _vnacal_teardown_parameter_collection(vnacal_t *vcp)
{
    vnacal_parameter_collection_t *vprmcp = &vcp->vc_parameter_collection;

    for (int i = vprmcp->vprmc_allocation - 1; i >= 0; --i) {
	vnacal_parameter_t *vpmrp = vprmcp->vprmc_vector[i];

	if (vpmrp != NULL) {
	    assert(!vpmrp->vpmr_deleted);
	    vpmrp->vpmr_deleted = true;
	    _vnacal_release_parameter(vpmrp);
	    assert(vprmcp->vprmc_vector[i] == NULL);
	}
    }
    free((void *)vprmcp->vprmc_vector);
    (void)memset((void *)&vcp->vc_parameter_collection, 0,
	    sizeof(vcp->vc_parameter_collection));
}
