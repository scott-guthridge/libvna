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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _VNACAL_NEW_INTERNAL_H
#define _VNACAL_NEW_INTERNAL_H

#include "vnacal_internal.h"
#include "vnacommon_internal.h"
#include "vnaerr_internal.h"
#include "vnaproperty.h"
#include "vnacal_layout.h"

#ifdef __cplusplus
extern "C" {
#endif

#define VNACAL_NEW_DEFAULT_P_TOLERANCE		1.0e-6

/*
 * vnacal_new_m_error_t: measurement error
 */
typedef struct vnacal_new_m_error {
    /* standard deviation of measurement noise */
    double vnme_noise;

    /* standard deviation of measurement tracking error */
    double vnme_tracking;

} vnacal_new_m_error_t;

/*
 * vnacal_new_parameter_t: a parameter used in a new calibration
 */
typedef struct vnacal_new_parameter {
    /* pointer to corresponding vnacal_parameter_t structure */
    vnacal_parameter_t *vnpr_parameter;

    /* back pointer to vnacal_new_t structure */
    vnacal_new_t *vnpr_cmp;

    /* true if the parameter value is unknown and must be determined */
    bool vnpr_unknown;

    /* union keyed on vnpr_unknown */
    union {
	struct {
	    /* second index to cmprc_unknown_vector */
	    int unknown_index;

	    /* for correlated parameters, the other */
	    struct vnacal_new_parameter *correlate;

	    /* next unknown parameter */
	    struct vnacal_new_parameter *next_unknown;
	} vnpr_unknown;
    } u;

    /* next parameter in hash chain */
    struct vnacal_new_parameter *vnpr_hash_next;

} vnacal_new_parameter_t;

#define vnpr_unknown_index	u.vnpr_unknown.unknown_index
#define vnpr_correlate		u.vnpr_unknown.correlate
#define vnpr_next_unknown	u.vnpr_unknown.next_unknown

/*
 * vnacal_new_parameter_hash_t: collection of new calibration parameters
 */
typedef struct vnacal_new_parameter_hash {
    /* hash table keyed on vnpr_parameter->vpmr_index */
    vnacal_new_parameter_t **vnph_table;

    /* number of hash table buckets */
    int vnph_allocation;

    /* number of elements stored in hash table */
    int vnph_count;

} vnacal_new_parameter_hash_t;

/*
 * vnacal_new_coefficient_t: single coefficient of an equation
 */
typedef struct vnacal_new_coefficient {
    /* column in expanded coefficient matrix, or -1 for "b" vector */
    int vnc_coefficient;

    /* multiply by -1 */
    bool vnc_negative;

    /* index into vnm_m_matrix or -1 if no m */
    int vnc_m_cell;

    /* index into vnm_s_matrix or -1 if no s */
    int vnc_s_cell;

    /* next term in equation */
    struct vnacal_new_coefficient *vnc_next;

} vnacal_new_coefficient_t;

/*
 * vnacal_new_equation_t: equation generated from a measured standard
 */
typedef struct vnacal_new_equation {
    /* associated measured calibration standard */
    struct vnacal_new_measurement *vne_vnmp;

    /* measurement row associated with this equation */
    int vne_row;

    /* measurement column associated with this equation */
    int vne_column;

    /* linked list of terms */
    vnacal_new_coefficient_t *vne_coefficient_list;

    /* next equation in system */
    struct vnacal_new_equation *vne_next;

} vnacal_new_equation_t;

/*
 * vnacal_new_measurement_t: a measured calibration standard
 */
typedef struct vnacal_new_measurement {
    /* index of this struct */
    int vnm_index;

    /* m_rows x m_columns matrix of vectors of per-frequency measurements
       of the standard */
    double complex **vnm_m_matrix;

    /* s_rows x s_columns matrix of the S parameters of the standard */
    vnacal_new_parameter_t **vnm_s_matrix;

    /* transitive closure of vnm_s_matrix */
    bool *vnm_reachability_matrix;

    /* associated vnacal_new_t structure */
    struct vnacal_new *vnm_vnp;

    /* next in list of measured calibration standards */
    struct vnacal_new_measurement *vnm_next;

} vnacal_new_measurement_t;

/*
 * vnacal_new_system_t: a linear system of equations
 */
typedef struct vnacal_new_system {
    /* count of equations in this system */
    int vns_equation_count;

    /* list of equations in this system */
    vnacal_new_equation_t *vns_equation_list;

    /* location where next equation should be linked */
    vnacal_new_equation_t **vns_equation_anchor;

} vnacal_new_system_t;

/*
 * vnacal_new_t: a system of calibration measurements
 */
struct vnacal_new {
    /* magic number used to avoid double-free */
    uint32_t vn_magic;

    /* pointer to vnacal_t structure */
    vnacal_t *vn_vcp;

    /* error parameter type and layout */
    vnacal_layout_t vn_layout;

    /* number of frequencies */
    int vn_frequencies;

    /* vector of frequencies */
    double *vn_frequency_vector;

    /* true if the frequency vector has been set */
    bool vn_frequencies_valid;

    /* hash table of parameters used here */
    vnacal_new_parameter_hash_t vn_parameter_hash;

    /* constant zero parameter */
    vnacal_new_parameter_t *vn_zero;

    /* number of unknown parameters */
    int vn_unknown_parameters;

    /* number of unknown correlated parameters */
    int vn_correlated_parameters;

    /* linked list of unknown parameters */
    vnacal_new_parameter_t *vn_unknown_parameter_list;

    /* point where next unknown parameter should be linked */
    vnacal_new_parameter_t **vn_unknown_parameter_anchor;

    /* system impedance of the VNA ports (currently assumed all the same) */
    double complex vn_z0;

    /* vector of measurement error values */
    vnacal_new_m_error_t *vn_m_error_vector;

    /* consider iterative solve converged when RMS p change less than this */
    double vn_p_tolerance;

    /* number of linear systems */
    int vn_systems;

    /* vector of linear systems of equations */
    vnacal_new_system_t *vn_system_vector;

    /* total number of equations in all systems */
    int vn_equations;

    /* maximum number of equations in any system */
    int vn_max_equations;

    /* list of measured standards */
    vnacal_new_measurement_t *vn_measurement_list;

    /* address where next measured standard should be added */
    vnacal_new_measurement_t **vn_measurement_anchor;

    /* count of measured standards */
    int vn_measurement_count;

    /* solved error parameters */
    vnacal_calibration_t *vn_calibration;

    /* vector of RMS error of the solutions by frequency */
    double *vn_rms_error_vector;

    /* next and previous elements in list of vnacal_new_t structure */
    list_t vn_next;
};

/*
 * vnacal_new_add_arguments_t: common arguments for _vnacal_new_add_common
 */
typedef struct vnacal_new_add_arguments {
    /* name of user called function */
    const char			       *vnaa_function;

    /* associated vnacal_new_t structure */
    vnacal_new_t		       *vnaa_cmp;

    /* matrix of voltages leaving each VNA port */
    const double complex	*const *vnaa_a_matrix;
    int					vnaa_a_rows;
    int					vnaa_a_columns;

    /* matrix of voltages entering each VNA port */
    const double complex	*const *vnaa_b_matrix;
    int					vnaa_b_rows;
    int					vnaa_b_columns;

    /* scattering parameters for the measured standard */
    const int			       *vnaa_s_matrix;
    int					vnaa_s_rows;
    int					vnaa_s_columns;

    /* true if m, s is only the diagonal */
    bool				vnaa_m_is_diagonal;
    bool				vnaa_s_is_diagonal;

    /* measurement type: 'a' for a/b or 'm' for m */
    char				vnaa_m_type;

    /* map from S port to VNA port */
    const int			       *vnaa_s_port_map;

} vnacal_new_add_arguments_t;


/* _vnacal_new_get_parameter: add/find parameter and return held */
extern vnacal_new_parameter_t *_vnacal_new_get_parameter(
	const char *function, vnacal_new_t *vnp, int parameter);

/* _vnacal_new_init_parameter_hash: set up the parameter hash */
extern int _vnacal_new_init_parameter_hash(const char *function,
	vnacal_new_parameter_hash_t *vnphp);

/* _vnacal_new_free_parameter_hash: free the parameter hash */
extern void _vnacal_new_free_parameter_hash(
	vnacal_new_parameter_hash_t *vnphp);

/* free a vnacal_new_measurement_t structure */
extern void _vnacal_new_free_measurement(vnacal_new_measurement_t *vnmp);

/* _vnacal_new_check_all_frequency_ranges: check f range of all parameters */
extern int _vnacal_new_check_all_frequency_ranges(const char *function,
	vnacal_new_t *vnp, double fmin, double fmax);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_NEW_INTERNAL_H */
