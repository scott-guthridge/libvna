/*
 * Vector Network Analyzer Calibration Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * vnacal_get_filename: return the calibration file name
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *
 * Return:
 *   NULL if the vnacal_t object came from vnacal_create and vnacal_save
 *   hasn't get been called.
 */
const char *vnacal_get_filename(const vnacal_t *vcp)
{
    return vcp->vc_filename;
}

/*
 * vnacal_get_sets: return the number of calibration sets
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 */
int vnacal_get_sets(const vnacal_t *vcp)
{
    return vcp->vc_sets;
}

/*
 * vnacal_get_setname: return the calibration set setname
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
const char *vnacal_get_setname(const vnacal_t *vcp, int set)
{
    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return NULL;
    }
    return vcp->vc_set_vector[set]->ets_setname;
}

/*
 * vnacal_get_rows: return the number of rows in the calibration matrix
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
int vnacal_get_rows(const vnacal_t *vcp, int set)
{
    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return -1;
    }
    return vcp->vc_set_vector[set]->ets_rows;
}

/*
 * vnacal_get_columns: return the number of columns in the calibration matrix
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
int vnacal_get_columns(const vnacal_t *vcp, int set)
{
    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return -1;
    }
    return vcp->vc_set_vector[set]->ets_columns;
}

/*
 * vnacal_get_frequencies: return the number of frequency points
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
int vnacal_get_frequencies(const vnacal_t *vcp, int set)
{
    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return -1;
    }
    return vcp->vc_set_vector[set]->ets_frequencies;
}

/*
 * vnacal_get_fmin: return the minimum frequency point
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
double vnacal_get_fmin(const vnacal_t *vcp, int set)
{
    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    return vcp->vc_set_vector[set]->ets_frequency_vector[0];
}

/*
 * vnacal_get_fmax: return the maximum frequency point
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
double vnacal_get_fmax(const vnacal_t *vcp, int set)
{
    vnacal_etermset_t *etsp;

    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    etsp = vcp->vc_set_vector[set];
    return etsp->ets_frequency_vector[etsp->ets_frequencies - 1];
}

/*
 * vnacal_get_fmax: return a pointer to the calibration frequency vector
 *   @vcp: a pointer to the object returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
const double *vnacal_get_frequency_vector(const vnacal_t *vcp, int set)
{
    if (set < 0 || set >= vcp->vc_sets) {
	errno = EINVAL;
	return NULL;
    }
    return vcp->vc_set_vector[set]->ets_frequency_vector;
}
