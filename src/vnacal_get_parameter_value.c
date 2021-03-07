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
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/*
 * _vnacal_get_parameter_frange: get frequency limits for the given parameter
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 *   @fmin: address of double to receive minimum
 *   @fmax: address of double to receive maximum
 */
void _vnacal_get_parameter_frange(vnacal_parameter_t *vpmrp,
	double *fmin, double *fmax)
{
    for (;;) {
	switch (vpmrp->vpmr_type) {
	case VNACAL_NEW:
	    break;

	case VNACAL_SCALAR:
	    *fmin = 0.0;
	    *fmax = INFINITY;
	    return;

	case VNACAL_VECTOR:
	    *fmin = vpmrp->vpmr_frequency_vector[0];
	    *fmax = vpmrp->vpmr_frequency_vector[vpmrp->vpmr_frequencies - 1];
	    return;

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
 * _vnacal_get_parameter_value: get the value of the parameter at f
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 *   @frequency: frequency at which to evaluate the value
 */
double complex _vnacal_get_parameter_value(vnacal_parameter_t *vpmrp,
	double frequency)
{
    for (;;) {
	switch (vpmrp->vpmr_type) {
	case VNACAL_NEW:
	    break;

	case VNACAL_SCALAR:
	    return vpmrp->vpmr_gamma;

	case VNACAL_VECTOR:
	    return _vnacal_rfi(vpmrp->vpmr_frequency_vector,
		    vpmrp->vpmr_gamma_vector,
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
