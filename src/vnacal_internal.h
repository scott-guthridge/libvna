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

#ifndef _VNACAL_INTERNAL_H
#define _VNACAL_INTERNAL_H

#include "vnacal.h"
#include "vnacommon_internal.h"
#include "vnaerr_internal.h"
#include "vnaproperty.h"
#include "vnacal_layout.h"

#ifdef __cplusplus
extern "C" {
#endif

/* number of predefined parameters */
#define VNACAL_PREDEFINED_PARAMETERS		3

/* default numerical precision for vnacal_save frequencies (in digits) */
#define VNACAL_DEFAULT_FREQUENCY_PRECISION	7

/* default numerical precision for vnacal_save data (in digits) */
#define VNACAL_DEFAULT_DATA_PRECISION		6

/* maximum points to use for rational function interpolation */
#define VNACAL_MAX_M				5

/* factor by which we will extrapolate the frequency */
#define VNACAL_F_EXTRAPOLATION			0.01

/*
 * magic numbers to protect against bad pointers from the user
 */
#define VN_MAGIC				0x564E4557 /* VNEW */
#define VC_MAGIC				0x5643414C /* VCAL */

/*
 * vnacal_parameter_type_t: type of parameter
 */
typedef enum vnacal_parameter_type {
    VNACAL_NEW,
    VNACAL_SCALAR,
    VNACAL_VECTOR,
    VNACAL_UNKNOWN,
    VNACAL_CORRELATED
} vnacal_parameter_type_t;

/*
 * vnacal_parameter_t: internal representation of a parameter
 */
typedef struct vnacal_parameter {
    vnacal_parameter_type_t vpmr_type;		/* parameter type */
    bool vpmr_deleted;				/* true if deleted */
    int vpmr_hold_count;			/* hold count */
    int vpmr_index;				/* position in vector */
    int vpmr_segment;				/* for _vnacal_rfi */
    vnacal_t *vpmr_vcp;				/* back pointer to vnacal_t */
    union {
	double complex scalar;			/* simple gamma value */
	struct {
	    int frequencies;			/* number of frequencies */
	    double *frequency_vector;		/* vector of frequencies */
	    double complex *gamma_vector;	/* vector of gamma values */
	} vector;
	struct {
	    struct vnacal_parameter *other;
	    union {
		struct {
		    double sigma;
		} correlated;
	    } u;
	} unknown;
    } u;
} vnacal_parameter_t;
#define vpmr_gamma		u.scalar
#define vpmr_frequencies	u.vector.frequencies
#define vpmr_frequency_vector	u.vector.frequency_vector
#define vpmr_gamma_vector	u.vector.gamma_vector
#define vpmr_other		u.unknown.other
#define vpmr_sigma		u.unknown.u.correlated.sigma

/*
 * VNACAL_GET_PARAMETER_TYPE: access the parameter type
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 */
#define VNACAL_GET_PARAMETER_TYPE(vpmrp)	((vpmrp)->vpmr_type)

/*
 * VNACAL_GET_PARAMETER_INDEX: access the parameter index
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 */
#define VNACAL_GET_PARAMETER_INDEX(vpmrp)	((vpmrp)->vpmr_index)

/*
 * VNACAL_GET_PARAMETER_OTHER: get the related parameter
 *   @vpmrp: pointer returned from _vnacal_get_parameter
 */
#define VNACAL_GET_PARAMETER_OTHER(vpmrp) \
    ((vpmrp)->vpmr_type == VNACAL_UNKNOWN || \
     (vpmrp)->vpmr_type == VNACAL_CORRELATED ? \
     ((vpmrp)->vpmr_other) : NULL)

/*
 * vnacal_parameter_collection_t: collection of parameters
 */
typedef struct vnacal_parameter_collection {
    /* allocation of vprmc_vector */
    int vprmc_allocation;

    /* number of non-NULL slots in vprmc_vector */
    int vprmc_count;

    /* first slot of vprmc_vector that might be free */
    int vprmc_first_free;

    /* vector of parameters */
    vnacal_parameter_t **vprmc_vector;

} vnacal_parameter_collection_t;

/*
 * vnacal_calibration_t: error terms
 */
typedef struct vnacal_calibration {
    /* calibration name filled in by _vnacal_add_calibration_common */
    char *cal_name;

    /* pointer to associated vnacal_t structure */
    vnacal_t *cal_vcp;

    /* type of error terms */
    vnacal_type_t cal_type;

    /* dimensions of the measurement matrix */
    int cal_rows;
    int cal_columns;

    /* number and vector of frequency values */
    int cal_frequencies;
    double *cal_frequency_vector;

    /* system impedance of the VNA ports (currently assumed all the same) */
    double complex cal_z0;

    /* vector, one per error term, of vectors of values, one per frequency */
    int cal_error_terms;
    double complex **cal_error_term_vector;

    /* per-calibration properties */
    vnaproperty_t *cal_properties;

} vnacal_calibration_t;

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
	/* known value of the parameter at the current frequency */
	double complex known_value;

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
 * vnacal_new_term_t: single term of an equation
 */
typedef struct vnacal_new_term {
    /* column in expanded coefficient matrix, or -1 for "b" vector */
    int vnt_coefficient;

    /* multiply by -1 */
    bool vnt_negative;

    /* index into vnm_m_matrix or -1 if no m */
    int vnt_measurement;

    /* index into vnm_s_matrix or -1 if no s */
    int vnt_sparameter;

    /* next term in equation */
    struct vnacal_new_term *vnt_next;

} vnacal_new_term_t;

/*
 * vnacal_new_equation_t: equation generated from a measured standard
 */
typedef struct vnacal_new_equation {
    /* associated measured calibration standard */
    struct vnacal_new_measurement *vne_vcsp;

    /* equation row and column */
    int vne_row;
    int vne_column;

    /* linked list of terms */
    vnacal_new_term_t *vne_term_list;

    /* next equation in system */
    struct vnacal_new_equation *vne_next;

} vnacal_new_equation_t;

/*
 * vnacal_new_measurement_t: measurement of a calibration standard
 */
typedef struct vnacal_new_measurement {
    /* matrix of vectors of per-frequency measurements of the standard */
    double complex **vnm_m_matrix;

    /* matrix of refrences representing the S parameters of the standard */
    vnacal_new_parameter_t **vnm_s_matrix;

    /* transitive closure of vnm_s_matrix */
    bool *vnm_reachability_matrix;

    /* pointer to the associated vnacal_new_t structure */
    struct vnacal_new *vnm_ncp;

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

    /* constant zero parameter */
    vnacal_new_parameter_t *vn_zero;

    /* hash table of parameters used here */
    vnacal_new_parameter_hash_t vn_parameter_hash;

    /* number of unknown parameters */
    int vn_unknown_parameters;

    /* linked list of unknown parameters */
    vnacal_new_parameter_t *vn_unknown_parameter_list;

    /* point where next unknown parameter should be linked */
    vnacal_new_parameter_t **vn_unknown_parameter_anchor;

    /* vector (entry per frequency) of vectors of unknown parameter values */
    double complex **vn_unknown_parameter_vector;

    /* vector of frequencies */
    double *vn_frequency_vector;

    /* true if the frequency vector has been set */
    bool vn_frequencies_valid;

    /* system impedance of the VNA ports (currently assumed all the same) */
    double complex vn_z0;

    /* number of linear systems */
    int vn_systems;

    /* vector of linear systems of equations */
    vnacal_new_system_t *vn_system_vector;

    /* maximum number of equations in any system */
    int vn_max_equations;

    /* list of measurements of calibration standards */
    vnacal_new_measurement_t *vn_measurement_list;

    /* address where next measurement should be added */
    vnacal_new_measurement_t **vn_measurement_anchor;

    /* solved error parameters */
    vnacal_calibration_t *vn_calibration;

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

/*
 * vnacal_t: structure returned from vnacal_load
 */
struct vnacal {
    /* magic number */
    uint32_t vc_magic;

    /* user-supplied error callback or NULL */
    vnaerr_error_fn_t *vc_error_fn;

    /* user-supplied error callback argument or NULL */
    void *vc_error_arg;

    /* collection of parameters */
    vnacal_parameter_collection_t vc_parameter_collection;

    /* allocation of vc_calibration_vector */
    int vc_calibration_allocation;

    /* vector of pointers to vnacal_calibration_t */
    vnacal_calibration_t **vc_calibration_vector;

    /* calibration filename */
    char *vc_filename;

    /* precision for frequency values */
    int vc_fprecision;

    /* precision for data values */
    int vc_dprecision;

    /* global properties */
    vnaproperty_t *vc_properties;

    /* doubly linked ring list of vnacal_new_t structures */
    list_t vc_new_head;
};

/* report an error */
extern void _vnacal_error(const vnacal_t *vcp, vnaerr_category_t category,
	const char *format, ...)
#ifdef __GNUC__
__attribute__ ((__format__ (__printf__, 3, 4)))
#endif /* __GNUC__ */
;

/* _vnacal_type_to_name: convert type to type name */
extern const char *_vnacal_type_to_name(vnacal_type_t type);

/* _vnacal_name_to_type: convert type name to type */
extern vnacal_type_t _vnacal_name_to_type(const char *name);

/* _vnacal_layout: init the error term layout structure */
extern void _vnacal_layout(vnacal_layout_t *vlp, vnacal_type_t type,
	int m_rows, int m_columns);

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
extern void _vnacal_new_free_standard(vnacal_new_measurement_t *vnmp);

/* _vnacal_new_check_all_frequency_ranges: check f range of all parameters */
extern int _vnacal_new_check_all_frequency_ranges(const char *function,
	vnacal_new_t *vnp, double fmin, double fmax);

/* _vnacal_alloc: allocate a vnacal_new_t structure */
extern vnacal_t *_vnacal_alloc(const char *function,
	vnaerr_error_fn_t *error_fn, void *error_arg);

/* _vnacal_calibration_alloc: alloc vnacal_calibration */
extern vnacal_calibration_t *_vnacal_calibration_alloc(vnacal_t *vcp,
	vnacal_type_t type, int rows, int columns, int frequencies,
	int error_terms);

/* _vnacal_calibration_get_fmin_bound: get the lower frequency bound */
extern double _vnacal_calibration_get_fmin_bound(
	const vnacal_calibration_t *calp);

/* _vnacal_etermset_get_fmax_bound: get the upper frequency bound */
extern double _vnacal_calibration_get_fmax_bound(
	const vnacal_calibration_t *calp);

/* _vnacal_calibration_free: free the memory for a vnacal_calibration_t */
extern void _vnacal_calibration_free(vnacal_calibration_t *calp);

/* _vnacal_get_calibration: return the calibration at the given index */
extern vnacal_calibration_t *_vnacal_get_calibration(const vnacal_t *vcp,
	int ci);

/* _vnacal_add_calibration_common */
extern int _vnacal_add_calibration_common(const char *function, vnacal_t *vcp,
        vnacal_calibration_t *calp, const char *name);

/* _vnacal_get_parameter: return a pointer to the parameter */
extern vnacal_parameter_t *_vnacal_get_parameter(vnacal_t *vcp, int parameter);

/* _vnacal_alloc_parameter: allocate a vnacal_parameter and return index */
extern vnacal_parameter_t *_vnacal_alloc_parameter(const char *function,
	vnacal_t *vcp);

/* _vnacal_hold_parameter: increase the hold count on a parameter */
extern void _vnacal_hold_parameter(vnacal_parameter_t *vpmrp);

/* _vnacal_release_parameter: decrease the hold count and free if zero */
extern void _vnacal_release_parameter(vnacal_parameter_t *vpmrp);

/* _vnacal_get_parameter_frange: get the frequency limits for the parameter */
extern void _vnacal_get_parameter_frange(vnacal_parameter_t *vpmrp,
	double *fmin, double *fmax);

/* _vnacal_get_parameter_value: get the value of the parameter at frequency */
extern double complex _vnacal_get_parameter_value(vnacal_parameter_t *vpmrp,
	double frequency);

/* _vnacal_setup_parameter_collection: allocate the parameter collection */
extern int _vnacal_setup_parameter_collection(const char *function,
	vnacal_t *vcp);

/* _vnacal_teardown_parameter_collection: free the parameter collection */
extern void _vnacal_teardown_parameter_collection(vnacal_t *vcp);

/* _vnacal_rfi: apply rational function interpolation */
extern double complex _vnacal_rfi(const double *xp, double complex *yp,
	int n, int m, int *segment, double x);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_INTERNAL_H */
