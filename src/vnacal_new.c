/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2022 D Scott Guthridge <scott_guthridge@rompromity.net>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
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
 * vnacal_new_alloc: allocate a new calibration structure
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @type: error term type
 *   @rows: number of VNA ports where signal is detected
 *   @columns: number of VNA ports where signal is generated
 *   @frequencies: number of frequency points
 *
 *   On error, set errno and return NULL.
 */
vnacal_new_t *vnacal_new_alloc(vnacal_t *vcp, vnacal_type_t type,
	int m_rows, int m_columns, int frequencies)
{
    vnacal_new_t *vnp = NULL;
    int systems;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return NULL;
    }
    if (m_rows < 1 || m_columns < 1) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_new_alloc: calibration matrix must be at least 1x1");
	return NULL;
    }
    if (frequencies < 0) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_new_alloc: frequencies cannot be negative");
	return NULL;
    }

    /*
     * In T parameters, fail if the measurement matrix has more rows
     * than columns; in U parameters, fail if the measurement matrix
     * has fewer rows than columns.  Otherwise, we'll construct systems
     * with more equations than measurements.  In principle, these could
     * still be solved, but instead of choosing orthagonal standards and
     * measuring them, one would have to choose orthagonal measurements
     * and then find standards that realize them.  That's not what we're
     * trying to do, so don't allow it.
     */
    switch (type) {
    case VNACAL_T8:
    case VNACAL_TE10:
    case VNACAL_T16:
	if (m_rows > m_columns) {
	    _vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_alloc: "
		    "U parameters must be used when m_rows > m_columns");
	    return NULL;
	}
	break;

    case VNACAL_E12:
	type = _VNACAL_E12_UE14;
	/*FALLTHROUGH*/

    case VNACAL_U8:
    case VNACAL_UE10:
    case VNACAL_UE14:
    case VNACAL_U16:
	if (m_rows < m_columns) {
	    _vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_alloc: "
		    "T parameters must be used when m_rows < m_columns");
	    return NULL;
	}
	break;

    default:
	_vnacal_error(vcp, VNAERR_USAGE,
		"vnacal_new_alloc: invalid calibration type %d", (int)type);
	return NULL;
    }
    systems = VNACAL_IS_UE14(type) ? m_columns : 1;

    /*
     * Allocate and init the vnacal_new_t structure.
     */
    if ((vnp = malloc(sizeof(vnacal_new_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return NULL;
    }
    (void)memset((void *)vnp, 0, sizeof(vnacal_new_t));
    vnp->vn_magic = VN_MAGIC;
    vnp->vn_vcp = vcp;
    _vnacal_layout(&vnp->vn_layout, type, m_rows, m_columns);
    vnp->vn_frequencies = frequencies;
    if ((vnp->vn_frequency_vector = calloc(frequencies,
		    sizeof(double))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	vnacal_new_free(vnp);
	return NULL;
    }
    vnp->vn_frequencies_valid = false;
    if (_vnacal_new_init_parameter_hash(__func__,
		&vnp->vn_parameter_hash) == -1) {
	vnacal_new_free(vnp);
	return NULL;
    }
    if ((vnp->vn_zero = _vnacal_new_get_parameter(__func__, vnp,
		    VNACAL_ZERO)) == NULL) {
	vnacal_new_free(vnp);
	return NULL;
    }
    vnp->vn_unknown_parameters = 0;
    vnp->vn_correlated_parameters = 0;
    vnp->vn_unknown_parameter_list = NULL;
    vnp->vn_unknown_parameter_anchor = &vnp->vn_unknown_parameter_list;
    vnp->vn_z0 = VNADATA_DEFAULT_Z0;
    vnp->vn_m_error_vector = NULL;
    vnp->vn_systems = systems;
    if ((vnp->vn_system_vector = calloc(systems,
		    sizeof(vnacal_new_system_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	vnacal_new_free(vnp);
	return NULL;
    }
    for (int i = 0; i < systems; ++i) {
	vnacal_new_system_t *vnsp = &vnp->vn_system_vector[i];

	vnsp->vns_equation_anchor = &vnsp->vns_equation_list;
    }
    vnp->vn_equations = 0;
    vnp->vn_max_equations = 0;
    vnp->vn_measurement_list = NULL;
    vnp->vn_measurement_anchor = &vnp->vn_measurement_list;
    vnp->vn_calibration = NULL;
    vnp->vn_rms_error_vector = NULL;

    /*
     * Link this structure onto the vnacal_t structure.
     */
    insque((void *)&vnp->vn_next, (void *)&vcp->vc_new_head);

    return vnp;
}

/*
 * vnacal_new_set_frequency_vector: set the frequency vector
 *   @vnp: pointer to vnacal_new_t structure
 *   @frequency_vector: vector of increasing frequencies
 */
int vnacal_new_set_frequency_vector(vnacal_new_t *vnp,
	const double *frequency_vector)
{
    vnacal_t *vcp;

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vcp = vnp->vn_vcp;
    if (frequency_vector == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE, "vnacal_new_set_frequency_vector: "
		"frequency_vector cannot be NULL");
	return -1;
    }
    for (int i = 0; i < vnp->vn_frequencies; ++i) {
	if (isnan(frequency_vector[i]) || frequency_vector[i] < 0.0) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_new_set_frequency_vector: invalid frequency: %f",
		    frequency_vector[i]);
	    return -1;
	}
    }
    for (int i = 0; i < vnp->vn_frequencies - 1; ++i) {
	if (frequency_vector[i] >= frequency_vector[i + 1]) {
	    _vnacal_error(vcp, VNAERR_USAGE,
		    "vnacal_new_set_frequency_vector: "
		    "frequencies must be ascending");
	    return -1;
	}
    }
    if (_vnacal_new_check_all_frequency_ranges(__func__, vnp,
		frequency_vector[0],
		frequency_vector[vnp->vn_frequencies - 1]) == -1) {
	return -1;
    }
    (void)memcpy((void *)vnp->vn_frequency_vector, (void *)frequency_vector,
	    vnp->vn_frequencies * sizeof(double));
    vnp->vn_frequencies_valid = true;
    return 0;
}

/*
 * vnacal_new_set_z0: set the system impedance for all VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @z0: nominal impedance looking into a VNA port
 *
 * Note:
 *   In this implementation, all VNA ports must have the same system
 *   impedance.	 If not set, the default is 50 ohms.
 */
int vnacal_new_set_z0(vnacal_new_t *vnp, double complex z0)
{
    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    vnp->vn_z0 = z0;
    return 0;
}

/*
 * _vnacal_new_free_measurement: free the memory for an vnacal_new_measurement_t
 *   @vnmp: structure to free
 */
void _vnacal_new_free_measurement(vnacal_new_measurement_t *vnmp)
{
    if (vnmp != NULL) {
	vnacal_new_t *vnp = vnmp->vnm_ncp;
	const int m_rows    = vnp->vn_layout.vl_m_rows;
	const int m_columns = vnp->vn_layout.vl_m_columns;

	free((void *)vnmp->vnm_reachability_matrix);
	if (vnmp->vnm_m_matrix != NULL) {
	    for (int m_cell = 0; m_cell < m_rows * m_columns; ++m_cell) {
		free((void *)vnmp->vnm_m_matrix[m_cell]);
	    }
	    free((void *)vnmp->vnm_m_matrix);
	}
	free((void *)vnmp->vnm_s_matrix);
	free((void *)vnmp);
    }
}

/*
 * vnacal_new_free: free a vnacal_new_t structure
 *   @vnp: pointer to vnacal_new_t structure
 */
void vnacal_new_free(vnacal_new_t *vnp)
{
    if (vnp != NULL && vnp->vn_magic == VN_MAGIC) {
	vnacal_new_measurement_t *vnmp;

	remque((void *)&vnp->vn_next);
	_vnacal_calibration_free(vnp->vn_calibration);

	for (int i = 0; i < vnp->vn_systems; ++i) {
	    vnacal_new_system_t *vnsp = &vnp->vn_system_vector[i];

	    while (vnsp->vns_equation_list != NULL) {
		vnacal_new_equation_t *vnep = vnsp->vns_equation_list;

		vnsp->vns_equation_list = vnep->vne_next;
		while (vnep->vne_coefficient_list != NULL) {
		    vnacal_new_coefficient_t *vncp = vnep->vne_coefficient_list;

		    vnep->vne_coefficient_list = vncp->vnc_next;
		    free((void *)vncp);
		}
		free((void *)vnep);
	    }
	}
	free((void *)vnp->vn_system_vector);
	while ((vnmp = vnp->vn_measurement_list) != NULL) {
	    vnp->vn_measurement_list = vnmp->vnm_next;
	    _vnacal_new_free_measurement(vnmp);
	}
	free((void *)vnp->vn_m_error_vector);
	_vnacal_new_free_parameter_hash(&vnp->vn_parameter_hash);
	free((void *)vnp->vn_frequency_vector);
	vnp->vn_magic = -1;
	free((void *)vnp);
    }
}
