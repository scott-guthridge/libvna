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

#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"


/*
 * _vnacal_etermset_alloc: alloc vnacal_etermset
 *   @vcp: a pointer to the vnacal_t structure
 *   @setname: set setname
 *   @rows: number of VNA ports where signal is detected
 *   @columns: number of VNA ports where signal is generated
 *   @frequencies: number of frequency points
 *
 *   Allocate the internal calibration data set data structure.
 *
 * Return on NULL on error.
 */
vnacal_etermset_t *_vnacal_etermset_alloc(vnacal_t *vcp, const char *setname,
	int rows, int columns, int frequencies)
{
    vnacal_etermset_t *etsp;
    int ncells = rows * columns;

    etsp = (vnacal_etermset_t *)malloc(sizeof(vnacal_etermset_t));
    if (etsp == NULL) {
	_vnacal_error(vcp, "malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)etsp, 0, sizeof(vnacal_etermset_t));
    if ((etsp->ets_setname = strdup(setname)) == NULL) {
	_vnacal_error(vcp, "strdup: %s", strerror(errno));
	goto error;
    }
    etsp->ets_vcp = vcp;
    etsp->ets_rows = rows;
    etsp->ets_columns = columns;
    etsp->ets_frequencies = frequencies;
    etsp->ets_frequency_vector = calloc(frequencies, sizeof(double));
    if (etsp->ets_frequency_vector == NULL) {
	_vnacal_error(vcp, "calloc: %s", strerror(errno));
	goto error;
    }
    etsp->ets_error_term_matrix = calloc(ncells, sizeof(vnacal_error_terms_t));
    if (etsp->ets_error_term_matrix == NULL) {
	_vnacal_error(vcp, "calloc: %s", strerror(errno));
	goto error;
    }
    for (int i = 0; i < ncells; ++i) {
	vnacal_error_terms_t *etp = &etsp->ets_error_term_matrix[i];

	for (int j = 0; j < 3; ++j) {
	    etp->et_data_vectors[j] = (double complex *)calloc(frequencies,
		    sizeof(double complex));
	    if (etp->et_data_vectors[j] == NULL) {
		_vnacal_error(vcp, "calloc: %s", strerror(errno));
		goto error;
	    }
	}
    }
    return etsp;

error:
    _vnacal_etermset_free(etsp);
    return NULL;
}

/*
 * _vnacal_etermset_get_fmin_bound: get the lower frequency bound
 *   @etsp: pointer to vnacal_etermset_t
 */
double _vnacal_etermset_get_fmin_bound(const vnacal_etermset_t *etsp)
{
    return (1.0 - VNACAL_F_EXTRAPOLATION) * etsp->ets_frequency_vector[0];
}

/*
 * _vancal_etermset_get_fmax_bound: get the upper frequency bound
 *   @etsp: pointer to vnacal_etermset_t
 */
double _vnacal_etermset_get_fmax_bound(const vnacal_etermset_t *etsp)
{
    return (1.0 + VNACAL_F_EXTRAPOLATION) *
	etsp->ets_frequency_vector[etsp->ets_frequencies - 1];
}

/*
 * _vnacal_etermset_free: free the memory allocated in _vncacal_alloc_data_set
 *   @etsp: pointer to vnacal_etermset_t
 */
void _vnacal_etermset_free(vnacal_etermset_t *etsp)
{
    int ncells = etsp->ets_rows * etsp->ets_columns;

    if (etsp->ets_error_term_matrix != NULL) {
	for (int i = 0; i < ncells; ++i) {
	    vnacal_error_terms_t *etp = &etsp->ets_error_term_matrix[i];

	    for (int j = 0; j < 3; ++j) {
		free((void *)etp->et_data_vectors[j]);
	    }
	}
    }
    free((void *)etsp->ets_error_term_matrix);
    free((void *)etsp->ets_frequency_vector);
    free((void *)etsp->ets_setname);
    free((void *)etsp);
}
