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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_make_vector_parameter: create frequency-dependent parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @frequency_vector: vector of increasing frequency values
 *   @frequencies: length of frequency_vector and gamma_vector
 *   @gamma_vector: vector of per-frequency gamma values
 */
int vnacal_make_vector_parameter(vnacal_t *vcp,
	const double *frequency_vector, int frequencies,
	const double complex *gamma_vector)
{
    vnacal_parameter_t *vpmrp;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (frequencies < 1) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_vector_parameter: "
		"at least one frequency must be given");
	return -1;
    }
    if (frequency_vector == NULL || gamma_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_vector_parameter: "
		"frequency_vector and gamma_vector must be non-NULL");
	return -1;
    }
    if (frequency_vector[0] < 0.0) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_vector_parameter: "
		"frequencies must be nonnegative");
	return -1;
    }
    for (int i = 1; i < frequencies; ++i) {
	if (frequency_vector[i - 1] >= frequency_vector[i]) {
	    _vnacal_error(vcp, VNAERR_USAGE, "vnacal_make_vector_parameter: "
		    "frequencies must be ascending");
	    return -1;
	}
    }
    vpmrp = _vnacal_alloc_parameter("vnacal_make_vector_parameter", vcp);
    if (vpmrp == NULL) {
	return -1;
    }
    vpmrp->vpmr_type = VNACAL_VECTOR;
    vpmrp->vpmr_frequencies = frequencies;
    vpmrp->vpmr_frequency_vector = calloc(frequencies, sizeof(double));
    if (vpmrp->vpmr_frequency_vector == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	vpmrp->vpmr_deleted = true;
	_vnacal_release_parameter(vpmrp);
	return -1;
    }
    (void)memcpy((void *)vpmrp->vpmr_frequency_vector,
	    (void *)frequency_vector, frequencies * sizeof(double));
    vpmrp->vpmr_gamma_vector = calloc(frequencies, sizeof(double complex));
    if (vpmrp->vpmr_gamma_vector == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	vpmrp->vpmr_deleted = true;
	_vnacal_release_parameter(vpmrp);
	return -1;
    }
    (void)memcpy((void *)vpmrp->vpmr_gamma_vector,
	(void *)gamma_vector, frequencies * sizeof(double complex));

    return vpmrp->vpmr_index;
}
