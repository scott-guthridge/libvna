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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_get_parameter_value: evaluate a parameter at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 *   @frequency: frequency at which to evaluate parameter
 */
double complex vnacal_get_parameter_value(vnacal_t *vcp, int parameter,
	double frequency)
{
    vnacal_parameter_t *vpmrp;
    double fmin, fmax;
    double lower, upper;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    if ((vpmrp = _vnacal_get_parameter(vcp, parameter)) == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_get_parameter_value: "
		"invalid parameter");
	return HUGE_VAL;
    }
    switch (vpmrp->vpmr_type) {
    case VNACAL_SCALAR:
	return vpmrp->vpmr_coefficient;

    case VNACAL_UNKNOWN:
    case VNACAL_CORRELATED:
	if (vpmrp->vpmr_frequency_vector == NULL) {
	    _vnacal_error(vcp, VNAERR_USAGE, "vnacal_get_parameter_value: "
		    "unknown parameter value");
	    return HUGE_VAL;
	}
	break;

    default:
	break;
    }
    fmin = vpmrp->vpmr_frequency_vector[0];
    fmax = vpmrp->vpmr_frequency_vector[vpmrp->vpmr_frequencies - 1];
    lower = (1.0 - VNACAL_F_EXTRAPOLATION) * fmin;
    upper = (1.0 + VNACAL_F_EXTRAPOLATION) * fmax;
    if (frequency < lower || frequency > upper) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_get_parameter_value: "
		"frequency %e must be between %e and %e\n",
		frequency, fmin, fmax);
	return HUGE_VAL;
    }
    return _vnacal_rfi(vpmrp->vpmr_frequency_vector,
	    vpmrp->vpmr_coefficient_vector,
	    vpmrp->vpmr_frequencies,
	    MIN(vpmrp->vpmr_frequencies, VNACAL_MAX_M),
	    &vpmrp->vpmr_segment,
	    frequency);
}
