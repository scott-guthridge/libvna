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
 * vnacal_eval_parameter: evaluate a parameter at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 *   @frequency: frequency at which to evaluate the parameter
 *   @z0: reference impedance for the result
 */
double complex vnacal_eval_parameter(vnacal_t *vcp, int parameter,
        double frequency, double complex z0)
{
    vnacal_parameter_t *vpmrp = NULL;
    vnacal_parameter_matrix_map_t *vpmmp = NULL;
    double complex result;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    if ((vpmrp = _vnacal_get_parameter(vcp, parameter)) == NULL) {
	return HUGE_VAL;
    }
    if ((vpmmp = _vnacal_analyze_parameter_matrix(__func__, vcp,
		    &vpmrp, 1, 1, /*initial=*/false)) == NULL) {
	return HUGE_VAL;
    }
    if (_vnacal_eval_parameter_matrix_i(__func__, vpmmp, frequency,
	    &z0, &result) == -1) {
	result = HUGE_VAL;
	goto out;
    }
out:
    _vnacal_free_parameter_matrix_map(vpmmp);
    return result;
}
