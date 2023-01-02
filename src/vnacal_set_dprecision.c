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
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_set_dprecision: set the data precision for vnacal_save
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @precision: precision in decimal places (1..n) or VNACAL_MAX_PRECISION
 */
int vnacal_set_dprecision(vnacal_t *vcp, int precision)
{
    if (precision < 1) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_set_dprecision: precision must be at least 1");
	return -1;
    }
    vcp->vc_dprecision = precision;
    return 0;
}
