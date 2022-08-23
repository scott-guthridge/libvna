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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
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
#include "vnacal_new_internal.h"


/*
 * vnacal_new_set_m_error: set the VNA measurement error
 *   @vnp: pointer to vnacal_new_t structure
 *   @frequency_vector: vector of frequency points
 *   @frequencies: number of frequencies
 *   @noise_error_vector: vector of standard deviation of noise floor
 *   @tracking_error_vector: vector of standard deviation of tracking error
 */
int vnacal_new_set_m_error(vnacal_new_t *vnp,
	const double *frequency_vector, int frequencies,
	const double *noise_error_vector,
	const double *tracking_error_vector)
{
    vnacal_t *vcp;
    vnacal_new_m_error_t *m_error_vector = NULL;

    /*
     * Validate arguments.
     */
    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vcp = vnp->vn_vcp;
    m_error_vector = vnp->vn_m_error_vector;
    if (frequencies < 1) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_m_error: "
		"frequencies must be at least 1");
	return -1;
    }
    if (m_error_vector == NULL && noise_error_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_m_error: "
		"noise_error_vector must be non-NULL");
	return -1;
    }
    if (noise_error_vector == NULL) {
	for (int i = 0; i < frequencies; ++i) {
	    if (noise_error_vector[i] <= 0) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"vnacal_new_set_m_error: noise error "
			"values must be positive");
		return -1;
	    }
	}
    }
    if (tracking_error_vector == NULL) {
	for (int i = 0; i < frequencies; ++i) {
	    if (tracking_error_vector[i] < 0) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"vnacal_new_set_m_error: tracking error "
			"values must be non-negative");
		return -1;
	    }
	}
    }
    if (!vnp->vn_frequencies_valid) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_m_error: "
		"vnacal_new_set_frequency_vector must be called first");
	return -1;
    }
    if (frequency_vector != NULL) {
	double fmin, fmax;
	double lower, upper;

	for (int i = 1; i < frequencies; ++i) {
	    if (frequency_vector[i - 1] >= frequency_vector[i]) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"vnacal_new_set_m_error: "
			"frequencies must be ascending");
		return -1;
	    }
	}
	fmin = vnp->vn_frequency_vector[0];
	fmax = vnp->vn_frequency_vector[vnp->vn_frequencies - 1];
	lower = (1.0 + VNACAL_F_EXTRAPOLATION) * fmin;
	upper = (1.0 - VNACAL_F_EXTRAPOLATION) * fmax;
	if (frequency_vector[0] > lower ||
		frequency_vector[frequencies - 1] < upper) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_new_set_m_error: frequency range "
		    "%.3e..%.3e is outside of calibration range %.3e..%3.e",
		    frequency_vector[0], frequency_vector[frequencies - 1],
		    fmin, fmax);
	    return -1;
	}

    } else if (frequencies != 1 && frequencies != vnp->vn_frequencies) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_m_error: "
		"invalid NULL frequency_vector");
	return -1;
    }
    if (m_error_vector == NULL) {
	if ((m_error_vector = malloc(vnp->vn_frequencies *
			sizeof(vnacal_new_m_error_t))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	    return -1;
	}
	for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
	    m_error_vector[findex].vnme_noise    = 0.0;
	    m_error_vector[findex].vnme_tracking = 0.0;
	}
	vnp->vn_m_error_vector = m_error_vector;
    }

    /*
     * If frequencies is 1, ignore frequency vector and use
     * noise_error_vector[0] and tracking_error_vector[0] for all
     * calibration frequencies.
     */
    if (frequencies == 1) {
	if (noise_error_vector != NULL) {
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_noise =
		    noise_error_vector[0];
	    }
	}
	if (tracking_error_vector != NULL) {
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_tracking =
		    tracking_error_vector[0];
	    }
	}

    /*
     * Else if no frequency_vector was given, use the frequencies of
     * the calibration.
     */
    } else if (frequency_vector == NULL) {
	assert(frequencies == vnp->vn_frequencies);
	if (noise_error_vector != NULL) {
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_noise =
		    noise_error_vector[findex];
	    }
	}
	if (tracking_error_vector != NULL) {
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_tracking =
		    tracking_error_vector[findex];
	    }
	}

    /*
     * Otherwise, use cubic spline interpolation over the frequencies
     * given.
     */
    } else {
	double c_vector[frequencies - 1][3];

	if (noise_error_vector != NULL) {
	    if (_vnacommon_spline_calc(frequencies - 1, frequency_vector,
			noise_error_vector, c_vector) == -1) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s",
			strerror(errno));
		return -1;
	    }
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_noise =
		    _vnacommon_spline_eval(frequencies - 1, frequency_vector,
			    noise_error_vector, c_vector,
			    vnp->vn_frequency_vector[findex]);
	    }
	}
	if (tracking_error_vector != NULL) {
	    if (_vnacommon_spline_calc(frequencies - 1, frequency_vector,
			tracking_error_vector, c_vector) == -1) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s",
			strerror(errno));
		return -1;
	    }
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_tracking =
		    _vnacommon_spline_eval(frequencies - 1, frequency_vector,
			    tracking_error_vector, c_vector,
			    vnp->vn_frequency_vector[findex]);
	    }
	}
    }
    return 0;
}
