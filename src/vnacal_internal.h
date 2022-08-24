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
    /* parameter type */
    vnacal_parameter_type_t vpmr_type;

    /* true if vnacal_delete_parameter called */
    bool vpmr_deleted;

    /* reference count */
    int vpmr_hold_count;

    /* index into vprmc_vector */
    int vpmr_index;

    /* start position for _vnacal_rfi */
    int vpmr_segment;

    /* back pointer to vnacal_t */
    vnacal_t *vpmr_vcp;

    /* type-specific members */
    union {
	double complex scalar;
	struct {
	    /*
	     * For vector parameters, frequencies, frequency_vector and
	     * gamma_vector are set up in vnacal_make_vector_parameter;
	     * for unknown and correlated parameters, these are set up
	     * in vnacal_new_solve.
	     */
	    int frequencies;
	    double *frequency_vector;
	    double complex *gamma_vector;
	    struct {
		/* pointer to related parameter */
		struct vnacal_parameter *other;
		union {
		    struct {
			/* number of freq-dependent sigma values */
			int sigma_frequencies;

			/* sigma frequences; may == other->frequency_vec */
			double *sigma_frequency_vector;

			/* sigma_frequencies long vector of sigma values */
			double *sigma_vector;

			/* cubic spline coefficients for sigma_vector */
			double (*sigma_spline)[3];

		    } correlated;
		} u;
	    } unknown;	/* including correlated */
	} vector;
    } u;
} vnacal_parameter_t;

/*
 * Hide the union
 */
#define vpmr_gamma			u.scalar
#define vpmr_frequencies		u.vector.frequencies
#define vpmr_frequency_vector		u.vector.frequency_vector
#define vpmr_gamma_vector		u.vector.gamma_vector
#define vpmr_other			u.vector.unknown.other
#define vmpr_correlated			u.vector.unknown.u.correlated
#define vpmr_sigma_frequencies		vmpr_correlated.sigma_frequencies
#define vpmr_sigma_frequency_vector	vmpr_correlated.sigma_frequency_vector
#define vpmr_sigma_vector		vmpr_correlated.sigma_vector
#define vpmr_sigma_spline		vmpr_correlated.sigma_spline

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

/* _vnacal_layout: init the error term layout structure */
extern void _vnacal_layout(vnacal_layout_t *vlp, vnacal_type_t type,
	int m_rows, int m_columns);

/* _vnacal_alloc: allocate a vnacal_t structure */
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

/* _vnacal_get_parameter_value_i: get the value of the parameter at frequency */
extern double complex _vnacal_get_parameter_value_i(vnacal_parameter_t *vpmrp,
	double frequency);

/* _vnacal_get_correlated_sigma: return the sigma value for the given f */
extern double _vnacal_get_correlated_sigma(vnacal_parameter_t *vpmrp,
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
