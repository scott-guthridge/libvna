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
 *   @sigma_nf_vector: standard deviations of the noise floor for measurements
 *   @sigma_tr_vector: standard deviation of noise proportional to signal level
 */
int vnacal_new_set_m_error(vnacal_new_t *vnp,
	const double *frequency_vector, int frequencies,
	const double *sigma_nf_vector, const double *sigma_tr_vector)
{
    vnacal_t *vcp;
    const vnacal_layout_t *vlp;
    vnacal_new_m_error_t *m_error_vector = NULL;

    /*
     * Validate arguments.
     */
    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vcp = vnp->vn_vcp;
    vlp = &vnp->vn_layout;
    m_error_vector = vnp->vn_m_error_vector;
    if (frequencies < 1) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_m_error: "
		"frequencies must be at least 1");
	return -1;
    }

    /*
     * If both vectors are NULL, clear any previous measurement
     * error setting and return.
     */
    if (sigma_nf_vector == NULL && sigma_tr_vector == NULL) {
	free((void *)vnp->vn_m_error_vector);
	vnp->vn_m_error_vector = NULL;
	return 0;
    }

    /*
     * Continue validating...
     */
    if (sigma_nf_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_m_error: "
		"noise error required if gain error given");
	return -1;
    }
    for (int i = 0; i < frequencies; ++i) {
	if (sigma_nf_vector[i] <= 0) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_new_set_m_error: noise error "
		    "values must be positive");
	    return -1;
	}
    }
    if (sigma_tr_vector != NULL) {
	for (int i = 0; i < frequencies; ++i) {
	    if (sigma_tr_vector[i] < 0) {
		_vnacal_error(vcp, VNAERR_USAGE,
			"vnacal_new_set_m_error: gain error "
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

    /*
     * If the error term type is T16 or U16, validate that all standards
     * given so far fully specify the S matrix.  Full S is needed to
     * generate the V matrices in this case.
     */
    if (VL_TYPE(vlp) == VNACAL_T16 || VL_TYPE(vlp) == VNACAL_U16) {
       const int s_cells = VL_S_ROWS(vlp) * VL_S_COLUMNS(vlp);
       int measurement = 0;
       vnacal_new_measurement_t *vnmp;

       for (vnmp = vnp->vn_measurement_list; vnmp != NULL;
	       vnmp = vnmp->vnm_next) {
	   ++measurement;
	   for (int s_cell = 0; s_cell < s_cells; ++s_cell) {
	       if (vnmp->vnm_s_matrix[s_cell] == NULL) {
		   _vnacal_new_err_need_full_s(vnp, __func__,
			   measurement, s_cell);
		   return -1;
	       }
	   }
       }
    }

    /*
     * Allocate the vector if needed.
     */
    if (m_error_vector == NULL) {
	if ((m_error_vector = malloc(vnp->vn_frequencies *
			sizeof(vnacal_new_m_error_t))) == NULL) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s", strerror(errno));
	    return -1;
	}
	vnp->vn_m_error_vector = m_error_vector;
    }

    /*
     * Always init the vector.
     */
    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
	m_error_vector[findex].vnme_sigma_nf = 0.0;
	m_error_vector[findex].vnme_sigma_tr = 0.0;
    }

    /*
     * If frequencies is 1, ignore frequency vector and use
     * sigma_nf_vector[0] and sigma_tr_vector[0] for all
     * calibration frequencies.
     */
    if (frequencies == 1) {
	for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
	    vnp->vn_m_error_vector[findex].vnme_sigma_nf =
		sigma_nf_vector[0];
	}
	if (sigma_tr_vector != NULL) {
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_sigma_tr =
		    sigma_tr_vector[0];
	    }
	}

    /*
     * Else if no frequency_vector was given, use the frequencies of
     * the calibration.
     */
    } else if (frequency_vector == NULL) {
	assert(frequencies == vnp->vn_frequencies);
	for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
	    vnp->vn_m_error_vector[findex].vnme_sigma_nf =
		sigma_nf_vector[findex];
	}
	if (sigma_tr_vector != NULL) {
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_sigma_tr =
		    sigma_tr_vector[findex];
	    }
	}

    /*
     * Otherwise, use cubic spline interpolation over the frequencies
     * given.
     */
    } else {
	double c_vector[frequencies - 1][3];

	if (_vnacommon_spline_calc(frequencies - 1, frequency_vector,
		    sigma_nf_vector, c_vector) == -1) {
	    _vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s",
		    strerror(errno));
	    return -1;
	}
	for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
	    vnp->vn_m_error_vector[findex].vnme_sigma_nf =
		_vnacommon_spline_eval(frequencies - 1, frequency_vector,
			sigma_nf_vector, c_vector,
			vnp->vn_frequency_vector[findex]);
	}
	if (sigma_tr_vector != NULL) {
	    if (_vnacommon_spline_calc(frequencies - 1, frequency_vector,
			sigma_tr_vector, c_vector) == -1) {
		_vnacal_error(vcp, VNAERR_SYSTEM, "malloc: %s",
			strerror(errno));
		return -1;
	    }
	    for (int findex = 0; findex < vnp->vn_frequencies; ++findex) {
		vnp->vn_m_error_vector[findex].vnme_sigma_tr =
		    _vnacommon_spline_eval(frequencies - 1, frequency_vector,
			    sigma_tr_vector, c_vector,
			    vnp->vn_frequency_vector[findex]);
	    }
	}
    }
    return 0;
}
