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
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _vnacal_calibration_alloc: alloc vnacal_calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @rows: number of VNA ports where signal is detected
 *   @columns: number of VNA ports where signal is generated
 *   @frequencies: number of frequency points
 *   @error_terms: number of error terms
 *
 *   Allocate the internal calibration data structure.
 *
 * Return on NULL on error.
 */
vnacal_calibration_t *_vnacal_calibration_alloc(vnacal_t *vcp,
	vnacal_type_t type, int rows, int columns, int frequencies,
	int error_terms)
{
    vnacal_calibration_t *calp;

    calp = (vnacal_calibration_t *)malloc(sizeof(vnacal_calibration_t));
    if (calp == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)calp, 0, sizeof(vnacal_calibration_t));
    calp->cal_vcp = vcp;
    calp->cal_type = type;
    calp->cal_rows = rows;
    calp->cal_columns = columns;
    calp->cal_frequencies = frequencies;
    calp->cal_frequency_vector = calloc(frequencies, sizeof(double));
    if (calp->cal_frequency_vector == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto error;
    }
    calp->cal_error_term_vector = calloc(error_terms, sizeof(double complex *));
    if (calp->cal_error_term_vector == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto error;
    }
    calp->cal_error_terms = error_terms;
    for (int term = 0; term < error_terms; ++term) {
	if ((calp->cal_error_term_vector[term] = calloc(frequencies,
			sizeof(double complex))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    goto error;
	}
    }
    return calp;

error:
    _vnacal_calibration_free(calp);
    return NULL;
}

/*
 * _vnacal_calibration_get_fmin_bound: get the lower frequency bound
 *   @calp: pointer to calibration structure
 */
double _vnacal_calibration_get_fmin_bound(const vnacal_calibration_t *calp)
{
    return (1.0 - VNACAL_F_EXTRAPOLATION) * calp->cal_frequency_vector[0];
}

/*
 * _vnacal_etermset_get_fmax_bound: get the upper frequency bound
 *   @calp: pointer to calibration structure
 */
double _vnacal_calibration_get_fmax_bound(const vnacal_calibration_t *calp)
{
    return (1.0 + VNACAL_F_EXTRAPOLATION) *
	calp->cal_frequency_vector[calp->cal_frequencies - 1];
}

/*
 * _vnacal_calibration_free: free the memory for a vnacal_calibration_t
 *   @calp: pointer to calibration structure
 */
void _vnacal_calibration_free(vnacal_calibration_t *calp)
{
    if (calp != NULL) {
	(void)vnaproperty_delete(&calp->cal_properties, ".");
	for (int term = 0; term < calp->cal_error_terms; ++term) {
	    free((void *)calp->cal_error_term_vector[term]);
	}
	free((void *)calp->cal_error_term_vector);
	free((void *)calp->cal_frequency_vector);
	free((void *)calp->cal_name);
	free((void *)calp);
    }
}

/*
 * _vnacal_add_calibration_common: add/replace a calibration
 *   @function: name of user-called function
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @calp: pointer to calibration structure
 *   @name: name of new calibration
 *
 * Return the new index.
 */
int _vnacal_add_calibration_common(const char *function, vnacal_t *vcp,
	vnacal_calibration_t *calp, const char *name)
{
    int first_free = -1;
    int cur;

    /*
     * Search the existing calibrations for name.  If a calibration with
     * name exists, we'll replace it.
     */
    for (cur = 0; cur < vcp->vc_calibration_allocation; ++cur) {
	if (vcp->vc_calibration_vector[cur] != NULL) {
	    if (strcmp(vcp->vc_calibration_vector[cur]->cal_name, name) == 0) {
		break;
	    }
	} else {
	    if (first_free == -1) {
		first_free = cur;
	    }
	}
    }

    /*
     * Extend the allocation as needed.  If we found a free slot above,
     * use it.
     */
    if (cur == vcp->vc_calibration_allocation) {
	if (first_free > -1) {
	    cur = first_free;
	} else {
	    int new_allocation;
	    vnacal_calibration_t **calpp = NULL;

	    /*
	     * Most common number of calibrations is expected to be 1, so
	     * initially allocate only one slot.  If a second calibration
	     * is added, though, increase to 8 and double from there on.
	     */
	    switch (vcp->vc_calibration_allocation) {
	    case 0:
		new_allocation = 1;
		break;
	    case 1:
		new_allocation = 8;
		break;
	    default:
		new_allocation = 2 * vcp->vc_calibration_allocation;
		break;
	    }
	    calpp = realloc(vcp->vc_calibration_vector,
		    new_allocation * sizeof(vnacal_calibration_t *));
	    if (calpp == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "realloc: %s",
			strerror(errno));
		return -1;
	    }
	    (void)memset((void *)&calpp[vcp->vc_calibration_allocation], 0,
		    (new_allocation - vcp->vc_calibration_allocation) *
		    sizeof(vnacal_calibration_t *));
	    vcp->vc_calibration_vector = calpp;
	    vcp->vc_calibration_allocation = new_allocation;
	}
    }

    /*
     * Fill in the calibration name.
     */
    assert(calp->cal_name == NULL);
    if ((calp->cal_name = strdup(name)) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"strdup: %s", strerror(errno));
	return -1;
    }

    /*
     * If we're replacing, free the old calibration.
     */
    _vnacal_calibration_free(vcp->vc_calibration_vector[cur]);
    vcp->vc_calibration_vector[cur] = calp;
    return cur;
}
