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
 * vnacal_make_correlated_parameter: create unknown parameter related by sigma
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @other: another parameter close to this one
 *   @sigma_frequency_vector: vector of increasing frequency values
 *   @sigma_frequencies: length of sigma_frequency_vector and sigma_vector
 *   @sigma_vector: frequency dependent parameter deviation from other
 */
int vnacal_make_correlated_parameter(vnacal_t *vcp, int other,
	const double *sigma_frequency_vector, int sigma_frequencies,
	const double *sigma_vector)
{
    vnacal_parameter_t *vpmrp_other;
    vnacal_parameter_t *vpmrp;
    vnacal_parameter_t *vpmrp_end = NULL;
    double *frequency_vector_copy = NULL;
    double *sigma_vector_copy = NULL;
    double (*spline_vector)[3] = NULL;


    /*
     * Validate parameters.
     */
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vpmrp_other = _vnacal_get_parameter(vcp, other);
    if (vpmrp_other == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_correlated_parameter: "
		"other must refer to a valid scalar or vector parameter");
	return -1;
    }
    if (sigma_frequencies < 1) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_correlated_parameter: "
		"at least one frequency must be given");
	return -1;
    }
    if (sigma_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_correlated_parameter: "
		"sigma_vector cannot be NULL");
	return -1;
    }

    /*
     * If sigma_frequencies == 1, ignore sigma_frequency_vector.  If it's
     * greater than 1, handle it here.
     */
    if (sigma_frequencies > 1) {
	/*
	 * Follow the "other" pointers until we get to the VNACAL_SCALAR
	 * or VNACAL_VECTOR parameter serving as the initial guess.
	 */
	vpmrp_end = vpmrp_other;
	while (vpmrp_end->vpmr_type == VNACAL_UNKNOWN ||
		vpmrp_end->vpmr_type == VNACAL_CORRELATED) {
	    vpmrp_end = vpmrp_end->vpmr_other;
	}

	/*
	 * We allow a special case of sigma_frequency_vector == NULL to
	 * mean take the frequencies from the initial guess parameter.
	 * For this to work, the object at the end of the chain has to
	 * have type VNACAL_VECTOR, and the number of frequencies has
	 * to match.
	 *
	 * Otherwise, if sigma_frequency_vector is not NULL, then validate
	 * it and make a copy.
	 */
	if (sigma_frequency_vector == NULL) {
	    if (vpmrp_end->vpmr_type != VNACAL_VECTOR ||
		    sigma_frequencies != vpmrp_end->vpmr_frequencies) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"vnacal_make_correlated_parameter: "
			"sigma_frequency_vector can be NULL only if the "
			"initial guess is a vector parameter and the counts "
			"of frequencies are equal");
		goto error;
	    }

	} else {
	    /*
	     * Validate that the frequencies are non-negative and ascending.
	     */
	    if (sigma_frequency_vector[0] < 0.0) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"vnacal_make_correlated_parameter: "
			"frequencies must be nonnegative");
		goto error;
	    }
	    for (int i = 1; i < sigma_frequencies; ++i) {
		if (sigma_frequency_vector[i] <=
			sigma_frequency_vector[i - 1]) {
		    _vnacal_error(vcp, VNAERR_USAGE,
			    "vnacal_make_correlated_parameter: "
			    "frequencies must be ascending");
		    goto error;
		}
	    }

	    /*
	     * If the initial guess is a vector parameter, make sure the
	     * frequency ranges aren't disjoint.
	     */
	    if (vpmrp_end->vpmr_type == VNACAL_VECTOR) {
		int other_n = vpmrp_end->vpmr_frequencies;
		double other_min = vpmrp_end->vpmr_frequency_vector[0];
		double other_max = vpmrp_end->vpmr_frequency_vector[other_n-1];
		double sigma_min = sigma_frequency_vector[0];
		double sigma_max = sigma_frequency_vector[sigma_frequencies-1];

		if (sigma_min > other_max || sigma_max < other_min) {
		    _vnacal_error(vcp, VNAERR_USAGE,
			    "vnacal_make_correlated_parameter: "
			    "sigma_frequency_vector cannot be disjoint with "
			    "the initial guess");
		    goto error;
		}
	    }

	    /*
	     * Make a copy of sigma_frequency_vector.
	     */
	    frequency_vector_copy = calloc(sigma_frequencies, sizeof(double));
	    if (frequency_vector_copy == NULL) {
		_vnacal_error(vcp, VNAERR_SYSTEM,
			"calloc: %s", strerror(errno));
		goto error;
	    }
	    (void)memcpy((void *)frequency_vector_copy,
		    (void *)sigma_frequency_vector,
		    sigma_frequencies * sizeof(double));
	}
    }

    /*
     * Validate and copy sigma_vector.
     */
    for (int findex = 0; findex < sigma_frequencies; ++findex) {
	if (sigma_vector[findex] <= 0.0) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_make_correlated_parameter: "
		    "sigma values must be positive");
	    goto error;
	}
    }
    if ((sigma_vector_copy = calloc(sigma_frequencies,
		    sizeof(double))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto error;
    }
    (void)memcpy((void *)sigma_vector_copy, (void *)sigma_vector,
	    sigma_frequencies * sizeof(double));

    /*
     * If there are at least two frequencies, generate cubic spline
     * coefficients over sigma_frequency_vector and sigma_vector.
     */
    if (sigma_frequencies > 1) {
	spline_vector = calloc(sigma_frequencies - 1, sizeof(double [3]));
	if (spline_vector == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "calloc: %s", strerror(errno));
	    goto error;
	}
	if (_vnacommon_spline_calc(sigma_frequencies - 1,
		    frequency_vector_copy != NULL ?
		    frequency_vector_copy : vpmrp_end->vpmr_frequency_vector,
		    sigma_vector_copy, spline_vector) == -1) {
	    _vnacal_error(vcp, VNAERR_SYSTEM,
		    "malloc: %s", strerror(errno));
	    goto error;
	}
    }

    /*
     * Create the new parameter.
     */
    vpmrp = _vnacal_alloc_parameter("vnacal_make_correlated_parameter", vcp);
    if (vpmrp == NULL) {
	return -1;
    }
    _vnacal_hold_parameter(vpmrp_other);
    vpmrp->vpmr_type = VNACAL_CORRELATED;
    vpmrp->vpmr_other = vpmrp_other;
    vpmrp->vpmr_sigma_frequencies = sigma_frequencies;
    if (sigma_frequencies == 1) {
	vpmrp->vpmr_sigma_frequency_vector = NULL;
    } else if (frequency_vector_copy != NULL) {
	vpmrp->vpmr_sigma_frequency_vector = frequency_vector_copy;
    } else {
	assert(vpmrp_end->vpmr_type == VNACAL_VECTOR);
	vpmrp->vpmr_sigma_frequency_vector = vpmrp_end->vpmr_frequency_vector;
    }
    vpmrp->vpmr_sigma_vector = sigma_vector_copy;
    vpmrp->vpmr_sigma_spline = spline_vector;
    return vpmrp->vpmr_index;

error:
    free((void *)spline_vector);
    free((void *)sigma_vector_copy);
    free((void *)frequency_vector_copy);
    return -1;
}

/*
 * _vnacal_get_correlated_sigma: return the sigma value for the given frequency
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 *   @frequency: frequency at which to evaluate sigma
 */
double _vnacal_get_correlated_sigma(vnacal_parameter_t *vpmrp, double frequency)
{
    assert(vpmrp->vpmr_type == VNACAL_CORRELATED);

    if (vpmrp->vpmr_sigma_frequencies == 1) {
	return vpmrp->vpmr_sigma_vector[0];
    }
    return _vnacommon_spline_eval(vpmrp->vpmr_sigma_frequencies - 1,
	    vpmrp->vpmr_sigma_frequency_vector, vpmrp->vpmr_sigma_vector,
	    vpmrp->vpmr_sigma_spline, frequency);
}
