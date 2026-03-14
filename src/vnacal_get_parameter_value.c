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
double complex vnacal_get_parameter_value(vnacal_t *vcp,
	int parameter, double frequency)
{
    vnacal_parameter_t *vpmrp = NULL;
    vnacal_parameter_matrix_map_t *vpmmp = NULL;
    double complex result;
    const char *alt;

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
    if (vpmmp->vpmm_standard_rmap != NULL) {
        if (vpmrp->vpmr_stdp->std_ports == 1) {
            alt = "vnacal_eval_parameter";
        } else {
            alt = "vnacal_eval_parameter_matrix";
        }
        _vnacal_error(vpmrp->vpmr_vcp, VNAERR_USAGE,
                "%s: cannot be used on calkit, data or embedded "
                "parameters; use %s instead",
                __func__, alt);
        result = HUGE_VAL;
	goto out;
    }
    if (_vnacal_eval_parameter_matrix_i(__func__, vpmmp, frequency,
	    NULL, &result) == -1) {
	result = HUGE_VAL;
	goto out;
    }
out:
    _vnacal_free_parameter_matrix_map(vpmmp);
    return result;
}
