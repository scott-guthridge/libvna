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

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_get_z0_vector: return a calibration reference impedance vector
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @vector: caller-provided #ports-long buffer to receive result
 *   @max_entries: number of double complex entries in vector
 *   @f: frequency at which to evaluate
 *
 *   Copies #ports reference impedances into the caller-provided buffer.
 *   The buffer should have space for at least one double complex entry
 *   per VNA port.  The max_entries parameter gives the length of the
 *   provided buffer to guard against buffer overrun.
 *
 *   When the z0 type is VNACAL_Z0_SCALAR, the single reference impedance
 *   is duplicated for each port.  When it's VNACAL_Z0_VECTOR, the entries
 *   are copied into the user's buffer.  When it's VNACAL_Z0_MATRIX,
 *   then the function returns the reference impedances for the given
 *   frequency, interpolating if necessary.  The frequency argument is
 *   ignored if the z0 type is not VNACAL_Z0_MATRIX.
 *
 * Return:
 *   number of VNA ports (number of entries placed into vector), or
 *   -1 on error
 */
int vnacal_get_z0_vector(const vnacal_t *vcp, int ci,
	double complex *vector, int max_entries, double f)
{
    const vnacal_calibration_t *calp;
    int ports;
    int segment = 0;

    if ((calp = _vnacal_get_calibration(__func__, vcp, ci)) == NULL) {
	return -1;
    }
    ports = MAX(calp->cal_rows, calp->cal_columns);
    if (ports > max_entries) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: "
		"result vector must have space for at least %d entries",
		__func__, ports);
	return -1;
    }
    switch (calp->cal_z0_type) {
    case VNACAL_Z0_SCALAR:
	for (int port = 0; port < ports; ++port) {
	    vector[port] = calp->cal_z0;
	}
	break;

    case VNACAL_Z0_VECTOR:
	(void)memcpy((void *)vector, (void *)calp->cal_z0_vector,
		ports * sizeof(double complex));
	break;

    case VNACAL_Z0_MATRIX:
	if (f < _vnacal_calibration_get_fmin_bound(calp)) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: frequency out of bounds %.3e < %.3e", __func__,
		    f, calp->cal_frequency_vector[0]);
	    return -1;
	}
	if (f > _vnacal_calibration_get_fmax_bound(calp)) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "%s: frequency out of bounds %.3e > %.3e", __func__,
		    f, calp->cal_frequency_vector[calp->cal_frequencies - 1]);
	    return -1;
	}
	for (int port = 0; port < ports; ++port) {
	    vector[port] = _vnacal_rfi(calp->cal_frequency_vector,
		    calp->cal_z0_matrix[port], calp->cal_frequencies,
		    MIN(calp->cal_frequencies, VNACAL_MAX_M), &segment, f);
	}
	break;

    default:
	abort();
    }
    return ports;
}
