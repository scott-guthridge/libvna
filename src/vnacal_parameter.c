/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2026 D Scott Guthridge <scott_guthridge@rompromity.net>
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
 * _vnacal_format_sxx: write the name of an S parameter into buffer
 *   @cur: buffer of at least size SXX_BUFFER_ALLOC
 *   @end: one past the end of buffer
 *   @row: zero-based row
 *   @column: zero-based column
 */
char *_vnacal_format_sxx(char *cur, char *end, int row, int column)
{
    if (cur >= end) {
	return end;
    }
    (void)snprintf(cur, end - cur, "S%d%s%d",
	    row + 1,
	    row > 8 || column > 8 ? "_" : "",
	    column + 1);
    end[-1] = '\000';
    cur += strlen(cur);
    return cur;
}

/*
 * _vnacal_get_parameter_name: copy descriptive name for parameter into buffer
 *   @vpmrp: parameter struct
 *   @with_sxx: include "sxx of" in front of standard parameters
 *   @buffer: buffer of at least PARAMETER_BUFFER_ALLOC + 1 chars for result
 */
void _vnacal_get_parameter_name(const vnacal_parameter_t *vpmrp, bool with_sxx,
	char *buffer)
{
    char *cur = buffer;
    char *end = &buffer[PARAMETER_BUFFER_ALLOC];

    switch (vpmrp->vpmr_type) {
    case VNACAL_NEW:
    default:
	abort();

    case VNACAL_SCALAR:
	(void)snprintf(cur, end - cur, "scalar(%f",
		creal(vpmrp->vpmr_coefficient));
	end[-1] = '\000';
	cur += strlen(cur);
	if (cimag(vpmrp->vpmr_coefficient) != 0.0) {
	    (void)snprintf(cur, end - cur, "%+fj",
		    cimag(vpmrp->vpmr_coefficient));
	}
	end[-1] = '\000';
	cur += strlen(cur);
	(void)stpecpy(cur, end, ") parameter");
	return;

    case VNACAL_VECTOR:
	(void)stpecpy(cur, end, "vector parameter");
	return;

    case VNACAL_UNKNOWN:
	(void)stpecpy(cur, end, "unknown parameter");
	return;

    case VNACAL_CORRELATED:
	(void)stpecpy(cur, end, "correlated parameter");
	return;

    case VNACAL_CALKIT:
    case VNACAL_DATA:
    case VNACAL_DEEMBED:
    case VNACAL_EMBED:
	{
	    vnacal_standard_t *stdp = vpmrp->vpmr_stdp;

	    /*
	     * If the standard has more than one port and with_sxx was
	     * given, start with "Sxx of ", where xx describes the S
	     * parameter of the standard.
	     */
	    if (stdp->std_ports > 1 && with_sxx) {
		cur = _vnacal_format_sxx(cur, end,
			vpmrp->vpmr_row, vpmrp->vpmr_column);
		cur = stpecpy(cur, end, " of ");
	    }
	    (void)stpecpy(cur, end, stdp->std_name);
	}
	return;
    }
}

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

	case VNACAL_CALKIT:
	    {
		vnacal_standard_t *stdp = vpmrp->vpmr_stdp;
		vnacal_calkit_standard_t *cstdp =
		    (vnacal_calkit_standard_t *)stdp;

		*fmin = cstdp->cstd_calkit_data.vcd_fmin;
		*fmax = cstdp->cstd_calkit_data.vcd_fmax;
	    }
	    break;

	case VNACAL_DATA:
	    {
		vnacal_standard_t *stdp = vpmrp->vpmr_stdp;
		vnacal_data_standard_t *dstdp = (vnacal_data_standard_t *)stdp;

		*fmin = dstdp->dstd_frequency_vector[0];
		*fmax = dstdp->dstd_frequency_vector[
		    dstdp->dstd_frequencies - 1];
	    }
	    break;

	case VNACAL_EMBED:
	case VNACAL_DEEMBED:
	    {
		vnacal_standard_t *stdp = vpmrp->vpmr_stdp;
		vnacal_embed_standard_t *estdp =
		    (vnacal_embed_standard_t *)stdp;

		*fmin = estdp->estd_fmin;
		*fmax = estdp->estd_fmax;
	    }
	    break;

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
    case VNACAL_NEW:
    case VNACAL_SCALAR:
	break;

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

    case VNACAL_CALKIT:
    case VNACAL_DATA:
    case VNACAL_EMBED:
    case VNACAL_DEEMBED:
	_vnacal_release_standard(&vpmrp->vpmr_stdp);
	break;

    default:
	abort();
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

#ifdef DEBUG
/*
 * print_parameter_collection: show the parameter collection for debug
 *   @vprmcp: parameter collection
 */
void print_parameter_collection(const vnacal_parameter_collection_t *vprmcp)
{
    char name[PARAMETER_BUFFER_ALLOC + 1];

    (void)printf("parameter collection:\n");
    for (int slot = 0; slot < vprmcp->vprmc_allocation; ++slot) {
	vnacal_parameter_t *vpmrp = vprmcp->vprmc_vector[slot];
	vnacal_standard_t *stdp;
	vnacal_embed_standard_t *estdp;

	(void)printf("  slot %d: ", slot);
	if (vpmrp == NULL) {
	    (void)printf("free\n");
	    continue;
	}
	_vnacal_get_parameter_name(vpmrp, /*with_sxx=*/true, name);
	(void)printf("%s\n", name);
	(void)printf("    hold_count %d", vpmrp->vpmr_hold_count);
	if (vpmrp->vpmr_deleted) {
	    (void)printf(" (DELETED)");
	}
	(void)printf("\n");
	assert(vpmrp->vpmr_index == slot);

	switch (vpmrp->vpmr_type) {
	case VNACAL_NEW:
	case VNACAL_VECTOR:
	case VNACAL_UNKNOWN:
	case VNACAL_CORRELATED:
	default:
	    continue;

	case VNACAL_SCALAR:
	    (void)printf("    value %f%+f\n",
		    creal(vpmrp->vpmr_coefficient),
		    cimag(vpmrp->vpmr_coefficient));
	    continue;

	case VNACAL_CALKIT:
	case VNACAL_DATA:
	    stdp = vpmrp->vpmr_stdp;
	    (void)printf("    standard @%p refcount %d\n",
		    (void *)stdp, stdp->std_refcount);
	    continue;

	case VNACAL_EMBED:
	case VNACAL_DEEMBED:
	    stdp = vpmrp->vpmr_stdp;
	    estdp = (vnacal_embed_standard_t *)stdp;
	    (void)printf("    standard @%p refcount %d\n",
		    (void *)stdp, stdp->std_refcount);
	    (void)printf("    target matrix:\n");
	    for (int row = 0; row < stdp->std_ports; ++row) {
		for (int column = 0; column < stdp->std_ports; ++column) {
		    const int cell = stdp->std_ports * row + column;
		    vnacal_parameter_t *vpmrp_temp;

		    (void)printf("        %d %d ", row, column);
		    vpmrp_temp = estdp->estd_target_matrix[cell];
		    if (vpmrp_temp == NULL) {
			(void)printf("null\n");
			continue;
		    }
		    _vnacal_get_parameter_name(vpmrp_temp,
			    /*with_sxx=*/true, name);
		    (void)printf("slot %d \"%s\"\n",
			    vpmrp_temp->vpmr_index, name);
		}
	    }
	    (void)printf("    fixture matrix:\n");
	    for (int row = 0; row < 2 * stdp->std_ports; ++row) {
		for (int column = 0; column < 2 * stdp->std_ports; ++column) {
		    const int cell = 2 * stdp->std_ports * row + column;
		    vnacal_parameter_t *vpmrp_temp;

		    (void)printf("        %d %d ", row, column);
		    vpmrp_temp = estdp->estd_fixture_matrix[cell];
		    if (vpmrp_temp == NULL) {
			(void)printf("null\n");
			continue;
		    }
		    _vnacal_get_parameter_name(vpmrp_temp,
			    /*with_sxx=*/true, name);
		    (void)printf("slot %d \"%s\"\n",
			    vpmrp_temp->vpmr_index, name);
		}
	    }
	    continue;
	}
    }
    (void)fflush(stdout);
}
#endif

/*
 * _vnacal_teardown_parameter_collection: free the parameter collection
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
void _vnacal_teardown_parameter_collection(vnacal_t *vcp)
{
    vnacal_parameter_collection_t *vprmcp = &vcp->vc_parameter_collection;

#ifdef DEBUG
    print_parameter_collection(vprmcp);
#endif /* DEBUG */
    for (int i = vprmcp->vprmc_allocation - 1; i >= 0; --i) {
	vnacal_parameter_t *vpmrp = vprmcp->vprmc_vector[i];

	/*
	 * Delete any parameters the caller didn't delete.
	 */
	if (vpmrp != NULL && !vpmrp->vpmr_deleted) {
	    vpmrp->vpmr_deleted = true;
	    _vnacal_release_parameter(vpmrp);
	}
    }
#if DEBUG > 1
    print_parameter_collection(vprmcp);
#endif
    /* check for refcount/memory leaks */
    for (int i = vprmcp->vprmc_allocation - 1; i >= 0; --i) {
	assert(vprmcp->vprmc_vector[i] == NULL);
    }
    free((void *)vprmcp->vprmc_vector);
    (void)memset((void *)&vcp->vc_parameter_collection, 0,
	    sizeof(vcp->vc_parameter_collection));
}
