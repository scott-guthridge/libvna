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

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"

/*
 * _vnacal_get_calibration: return the calibration at the given index
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
vnacal_calibration_t *_vnacal_get_calibration(const vnacal_t *vcp, int ci)
{
    vnacal_calibration_t *calp;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    if (ci < 0 || ci >= vcp->vc_calibration_allocation) {
	errno = EINVAL;
	return NULL;
    }
    if ((calp = vcp->vc_calibration_vector[ci]) == NULL) {
	errno = EINVAL;
	return NULL;
    }
    return calp;
}

/*
 * vnacal_get_filename: return the calibration file name
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *
 * Return:
 *   NULL if the vnacal_t structure came from vnacal_create and vnacal_save
 *   hasn't get been called.
 */
const char *vnacal_get_filename(const vnacal_t *vcp)
{
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    return vcp->vc_filename;
}

/*
 * vnacal_get_calibration_end: return one past the highest calibration index
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
int vnacal_get_calibration_end(const vnacal_t *vcp)
{
    int n;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    n = vcp->vc_calibration_allocation;
    while (n > 0 && vcp->vc_calibration_vector[n - 1] == NULL) {
	--n;
    }
    return n;
}

/*
 * vnacal_get_name: return the calibration name
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
const char *vnacal_get_name(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	errno = EINVAL;
	return NULL;
    }
    return calp->cal_name;
}

/*
 * vnacal_get_type: return the type of error terms
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
vnacal_type_t vnacal_get_type(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	errno = EINVAL;
	return -1;
    }
    return calp->cal_type;
}

/*
 * vnacal_get_rows: return the number of rows in the calibration matrix
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
int vnacal_get_rows(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	errno = EINVAL;
	return -1;
    }
    return calp->cal_rows;
}

/*
 * vnacal_get_columns: return the number of columns in the calibration matrix
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
int vnacal_get_columns(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	errno = EINVAL;
	return -1;
    }
    return calp->cal_rows;
}

/*
 * vnacal_get_frequencies: return the number of frequency points
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
int vnacal_get_frequencies(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	return -1;
    }
    return calp->cal_frequencies;
}

/*
 * vnacal_get_fmin: return the minimum frequency point
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
double vnacal_get_fmin(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	return HUGE_VAL;
    }
    return calp->cal_frequency_vector[0];
}

/*
 * vnacal_get_fmax: return the maximum frequency point
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
double vnacal_get_fmax(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	return HUGE_VAL;
    }
    return calp->cal_frequency_vector[calp->cal_frequencies - 1];
}

/*
 * vnacal_get_fmax: return a pointer to the calibration frequency vector
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
const double *vnacal_get_frequency_vector(const vnacal_t *vcp, int ci)
{
    const vnacal_calibration_t *calp;

    if ((calp = _vnacal_get_calibration(vcp, ci)) == NULL) {
	return NULL;
    }
    return calp->cal_frequency_vector;
}
