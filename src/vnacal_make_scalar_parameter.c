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
 * vnacal_make_scalar_parameter: create frequency-independent parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @gamma: coefficient of reflection or transmission
 */
int vnacal_make_scalar_parameter(vnacal_t *vcp, double complex gamma)
{
    vnacal_parameter_t *vpmrp;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vpmrp = _vnacal_alloc_parameter("vnacal_make_scalar_parameter", vcp);
    if (vpmrp == NULL) {
	return -1;
    }
    vpmrp->vpmr_type = VNACAL_SCALAR;
    vpmrp->vpmr_gamma = gamma;
    return vpmrp->vpmr_index;
}
