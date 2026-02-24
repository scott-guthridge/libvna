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
#include "vnacal_new_internal.h"

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
    const double complex default_z0 = VNADATA_DEFAULT_Z0;

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
	goto error;
    }
    vnp->vn_frequencies_valid = false;
    if (_vnacal_new_init_parameter_hash(__func__,
		&vnp->vn_parameter_hash) == -1) {
	goto error;
    }
    if ((vnp->vn_zero = _vnacal_new_get_parameter(__func__, vnp,
		    VNACAL_ZERO)) == NULL) {
	goto error;
    }
    vnp->vn_unknown_parameters = 0;
    vnp->vn_correlated_parameters = 0;
    vnp->vn_unknown_parameter_list = NULL;
    vnp->vn_unknown_parameter_anchor = &vnp->vn_unknown_parameter_list;
    if (_vnacal_new_set_z0_vector(__func__, vnp, &default_z0, 1) == -1) {
	goto error;
    }
    vnp->vn_m_error_vector = NULL;
    vnp->vn_p_tolerance = VNACAL_NEW_DEFAULT_P_TOLERANCE;
    vnp->vn_et_tolerance = VNACAL_NEW_DEFAULT_ET_TOLERANCE;
    vnp->vn_iteration_limit = VNACAL_NEW_DEFAULT_ITERATION_LIMIT;
    vnp->vn_pvalue_limit = VNACAL_NEW_DEFAULT_PVALUE_LIMIT;
    vnp->vn_systems = systems;
    if ((vnp->vn_system_vector = calloc(systems,
		    sizeof(vnacal_new_system_t))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"calloc: %s", strerror(errno));
	goto error;
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

error:
    vnacal_new_free(vnp);
    return NULL;
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
 * _vnacal_new_set_z0_vector: set the reference impedances
 *   @function: name of user-called function
 *   @vnp: vnacal_new_t structure
 *   @z0_vector: vector of reference impedance values
 *   @length: length of z0_vector: 1, ports or ports * frequencies
 */
int _vnacal_new_set_z0_vector(const char *function, vnacal_new_t *vnp,
	const double complex *z0_vector, int length)
{
    vnacal_t *vcp = vnp->vn_vcp;
    int rows, columns, ports;
    vnacal_z0_type_t z0_type;
    double complex *clfp;

    if (vnp == NULL || vnp->vn_magic != VN_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    rows = VL_M_ROWS(&vnp->vn_layout);
    columns = VL_M_COLUMNS(&vnp->vn_layout);
    ports = MAX(rows, columns);
    if (length == 1) {
	z0_type = VNACAL_Z0_SCALAR;
    } else if (length == ports) {
	z0_type = VNACAL_Z0_VECTOR;
    } else if (length == ports * vnp->vn_frequencies) {
	z0_type = VNACAL_Z0_MATRIX;
    } else {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: z0_vector length must be 1, number of ports, or "
		"number of ports * frequencies", function);
	return -1;
    }

    /*
     * In the single z0 case, allocate a ports long vector and duplicate
     * the z0 value into every cell so that we always have a vector.
     */
    if ((clfp = malloc(MAX(length, ports) * sizeof(double complex))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM,
		"malloc: %s", strerror(errno));
	return -1;
    }
    if (z0_type == VNACAL_Z0_SCALAR) {
	for (int i = 0; i < ports; ++i) {
	    clfp[i] = z0_vector[0];
	}
    } else {
	(void)memcpy((void *)clfp, (void *)z0_vector,
		length * sizeof(double complex));
    }
    free((void *)vnp->vn_z0_vector);
    vnp->vn_z0_type = z0_type;
    vnp->vn_z0_vector = clfp;
    return 0;
}

/*
 * _vnacal_new_get_z0_vector: return the z0 vector for the given freq index
 *   @vnp: vnacal_new_t structure
 *   @findex: frequency index
 */
const double complex *_vnacal_new_get_z0_vector(vnacal_new_t *vnp, int findex)
{
    int rows, columns, ports;

    switch (vnp->vn_z0_type) {
    case VNACAL_Z0_SCALAR:
    case VNACAL_Z0_VECTOR:
	return vnp->vn_z0_vector;

    case VNACAL_Z0_MATRIX:
	rows = VL_M_ROWS(&vnp->vn_layout);
	columns = VL_M_COLUMNS(&vnp->vn_layout);
	ports = MAX(rows, columns);
	return &vnp->vn_z0_vector[findex * ports];

    default:
	abort();
    }
}

/*
 * _vnacal_new_free_measurement: free the memory for an vnacal_new_measurement_t
 *   @vnmp: structure to free
 */
void _vnacal_new_free_measurement(vnacal_new_measurement_t *vnmp)
{
    if (vnmp != NULL) {
	vnacal_new_t *vnp = vnmp->vnm_vnp;
	const int m_rows    = vnp->vn_layout.vl_m_rows;
	const int m_columns = vnp->vn_layout.vl_m_columns;

	free((void *)vnmp->vnm_connectivity_matrix);
	if (vnmp->vnm_m_matrix != NULL) {
	    for (int m_cell = 0; m_cell < m_rows * m_columns; ++m_cell) {
		free((void *)vnmp->vnm_m_matrix[m_cell]);
	    }
	    free((void *)vnmp->vnm_m_matrix);
	}
	_vnacal_free_parameter_matrix_map(vnmp->vnm_parameter_map);
	free((void *)vnmp->vnm_parameter_matrix);
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
		while (vnep->vne_term_list != NULL) {
		    vnacal_new_term_t *vntp = vnep->vne_term_list;

		    vnep->vne_term_list = vntp->vnt_next;
		    free((void *)vntp);
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
