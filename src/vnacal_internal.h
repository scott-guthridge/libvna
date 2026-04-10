/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2026 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#include <assert.h>

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
    VNACAL_CORRELATED,
    VNACAL_CALKIT,
    VNACAL_DATA,
    VNACAL_EMBED,
    VNACAL_DEEMBED
} vnacal_parameter_type_t;

/* forward */
typedef struct vnacal_standard vnacal_standard_t;
typedef struct vnacal_parameter vnacal_parameter_t;
typedef struct vnacal_parameter_matrix_map vnacal_parameter_matrix_map_t;

/*
 * vnacal_standard_ops_t: subclass operations on vnacal_standard_t
 */
typedef struct vnacal_standard_ops {
    /* type of standard */
    vnacal_parameter_type_t stdo_type;

    /* evaluate the standard at the given frequency */
    int (*stdo_eval)(struct vnacal_standard *stdp, const char *function,
	const double complex *z0_vector, double frequency,
	double complex *result_matrix);

    /* destruct the derived struct members */
    void (*stdo_free)(struct vnacal_standard *stdp);

} vnacal_standard_ops_t;

/*
 * vnacal_standard_t: base class for standards
 */
struct vnacal_standard {
    /* class-specific type and operations */
    const vnacal_standard_ops_t *std_ops;

    /* name of standard for error messages */
    char *std_name;

    /* number of ports (assumed square) */
    int std_ports;

    /* reference count */
    int std_refcount;

    /* back pointer to vnacal_t structure */
    vnacal_t *std_vcp;

};

/*
 * vnacal_calkit_standard_t: a calibration kit standard
 */
typedef struct vnacal_calkit_standard {
    /* vnacal_standard_t base class */
    vnacal_standard_t cstd_base;

    /* calkit parameters */
    vnacal_calkit_data_t cstd_calkit_data;

} vnacal_calkit_standard_t;

/*
 * vnacal_data_standard_t: a standard based on network parameter data
 */
typedef struct vnacal_data_standard {
    /* vnacal_standard_t base class */
    vnacal_standard_t dstd_base;

    /* number of frequency points */
    int dstd_frequencies;

    /* vector of frequency values */
    double *dstd_frequency_vector;

    /* indicates per-frequency reference impedances */
    bool dstd_has_fz0;

    /* most recent segment used in _vnacal_rfi */
    int dstd_segment;

    /* reference impedances */
    union {
	/* vector of reference impedances by port */
	double complex *dstd_z0_vector;

	/* vector by port of vector by frequency of reference impedances */
	double complex **dstd_z0_vector_vector;
    } u;

    /* serialized matrix of vectors (by frequency) of S parameters */
    double complex **dstd_data;

} vnacal_data_standard_t;

/*
 * vnacal_embed_standard_t: an embedding or de-embedding
 */
typedef struct vnacal_embed_standard {
    /* vnacal_standard_t base class */
    vnacal_standard_t estd_base;

    /* minimum allowed frequency */
    double estd_fmin;

    /* maximum allowed frequency */
    double estd_fmax;

    /* std_ports x std_ports matrix of parameter with references */
    vnacal_parameter_t **estd_target_matrix;

    /* parameter_matrix map for ves_target_matrix */
    vnacal_parameter_matrix_map_t *estd_target_map;

    /* fixture_ports x fixture_ports matrix of parameter with references */
    vnacal_parameter_t **estd_fixture_matrix;

    /* parameter matrix map for ves_fixture_matrix */
    vnacal_parameter_matrix_map_t *estd_fixture_map;

} vnacal_embed_standard_t;

/*
 * vnacal_parameter_t: internal representation of a parameter
 */
struct vnacal_parameter {
    /* parameter type */
    vnacal_parameter_type_t vpmr_type;

    /* true if vnacal_delete_parameter called */
    bool vpmr_deleted;

    /* reference count */
    int vpmr_hold_count;

    /* index into vprmc_vector */
    int vpmr_index;

    /* back pointer to vnacal_t */
    vnacal_t *vpmr_vcp;

    /* type-specific members */
    union {
	double complex scalar;
	struct {
	    /*
	     * For vector parameters, frequencies, frequency_vector and
	     * coefficient_vector are set up in vnacal_make_vector_parameter;
	     * for unknown and correlated parameters, these are set up
	     * in vnacal_new_solve.
	     */
	    int frequencies;
	    double *frequency_vector;
	    double complex *coefficient_vector;

	    /* start position for _vnacal_rfi */
	    int segment;

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
	struct {
	    /* associated n-port standard */
	    vnacal_standard_t *stdp;

	    /* row and column of standard matrix */
	    int row, column;

	} standard;
    } u;
};

/*
 * Hide the union
 */
#define vpmr_coefficient		u.scalar
#define vpmr_frequencies		u.vector.frequencies
#define vpmr_frequency_vector		u.vector.frequency_vector
#define vpmr_coefficient_vector		u.vector.coefficient_vector
#define vpmr_segment			u.vector.segment
#define vpmr_other			u.vector.unknown.other
#define vmpr_correlated			u.vector.unknown.u.correlated
#define vpmr_sigma_frequencies		vmpr_correlated.sigma_frequencies
#define vpmr_sigma_frequency_vector	vmpr_correlated.sigma_frequency_vector
#define vpmr_sigma_vector		vmpr_correlated.sigma_vector
#define vpmr_sigma_spline		vmpr_correlated.sigma_spline
#define vpmr_stdp			u.standard.stdp
#define vpmr_row			u.standard.row
#define vpmr_column			u.standard.column

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
 * VNACAL_IS_STANDARD_PARAMETER: test if a parameter uses vnacal_standard_t
 */
#define VNACAL_IS_STANDARD_PARAMETER(vpmrp) \
    ((vpmrp)->vpmr_type == VNACAL_CALKIT || \
     (vpmrp)->vpmr_type == VNACAL_DATA || \
     (vpmrp)->vpmr_type == VNACAL_EMBED || \
     (vpmrp)->vpmr_type == VNACAL_DEEMBED)

/*
 * The parameter matrix for a calibration standard that is given to the
 * vnacal_new_add_* functions can consist of regular scalar or vector
 * parameters, unknown or correlated parameters that the library must
 * solve for, or calkit and data standards.  Regular scalar or vector
 * parameters are just numbers to the library -- it's the caller's
 * responsibility to pass values that are consistent with the reference
 * impedance(s) of the calibration.  In order to evaluate calkit or
 * data parameters, however, we need to know the reference impedances
 * of the VNA ports, and this implies that there has to be a consistent
 * mapping between the ports of the standard and the ports of the VNA.
 * For example, suppose we have a four-port VNA and a two-port standard
 * described by a vnadata_t structure is connected to VNA ports 1 and 2,
 * a calkit open parameter is connected to VNA port 3, and a perfect short
 * described by VNACAL_SHORT (a regular scalar parameter) is connected
 * to port 4.  Suppose the S parameters of the data standard are called
 * d11, d12, d21, and d22; the S parameter of the open standard is o11,
 * and the perfect short is called s11.  A valid parameter matrix for
 * this would be:
 *
 *    d22 d21 0   0
 *    d12 d11 0   0
 *    0   0   s11 0
 *    0   0   0   m11
 *
 * Note that the two ports of the data standard are connected in reverse
 * so the port mapping is:
 *
 *    VNA port 1 : Data standard port 2
 *    VNA port 2 : Data standard port 1
 *    VNA port 3 : Open standard (port 1)
 *    VNA port 4 : Short standard (port 1)
 *
 * But an invalid parameter matrix could also have been given, e.g.:
 *
 *   d22 d12 0    0
 *   d21 o11 0    0
 *   0   0   0    0
 *   0   0   s11  0
 *
 * There are several problems here.  First, the elements of the data
 * port standard have an inconsistent mapping where s11 of the parameter
 * matrix suggests that VNA port 1 maps to data standard port 2, but
 * s12 and s21 suggest that VNA port 1 maps to data standard port 1
 * and VNA port 2 maps to data standard port 2.  VNA port 2 appears
 * to be connected to both the data standard and the open standard.
 * And single port standard s11 is not on the major diagonal.
 *
 * The following data structures represent the mapping between the ports
 * of the parameter matrix and those of the standards.
 *
 * Note that the vnacal_new_add_* functions allow another layer of port
 * remapping on top of this, where the ports of the parameter matrix
 * may be only a subset of VNA ports or may be connected out of order.
 */

/*
 * vnacal_standard_rmap_t: map of ports of a standard to those of the
 *	parameter matrix
 */
typedef struct vnacal_standard_rmap {
    /* pointer to standard */
    vnacal_standard_t *vsrm_stdp;

    /* map from port of standard to port of parameter matrix (zero-based) */
    int *vsrm_rmap_vector;

    /* cell of parameter matrix that created each entry in vsrm_rmap_vector */
    int *vsrm_cell_vector;

    /* next in list */
    struct vnacal_standard_rmap *vsrm_next;

} vnacal_standard_rmap_t;

/*
 * vnacal_parameter_rmap_t: location of a regular parameter in the parameter
 *	matrix
 */
typedef struct vnacal_parameter_rmap {
    /* pointer to ordinary parameter */
    vnacal_parameter_t *vprm_parameter;

    /* cell of parameter matrix this parameter fills */
    int vprm_cell;

    /* next in list */
    struct vnacal_parameter_rmap *vprm_next;

} vnacal_parameter_rmap_t;

/*
 * vnacal_parameter_matrix_map_t: how the parameter matrix maps to ordinary
 *     parameters, calkit standards and data standards
 */
struct vnacal_parameter_matrix_map {
    /* associated calibration structure */
    vnacal_t *vpmm_vcp;

    /* number of rows in parameter matrix */
    int vpmm_rows;

    /* number of columns in parameter matrix */
    int vpmm_columns;

    /* linked list of (multi-port) standards */
    vnacal_standard_rmap_t *vpmm_standard_rmap;

    /* linked list of ordinary parameters */
    vnacal_parameter_rmap_t *vpmm_parameter_rmap;

};

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

    /* reference impedances */
    vnacal_z0_type_t cal_z0_type;
    union {
	double complex cal_z0;
	double complex *cal_z0_vector;	/* vector[port] */
	double complex **cal_z0_matrix;	/* matrix[port][findex] */
    } u;

    /* vector, one per error term, of vectors of values, one per frequency */
    int cal_error_terms;
    double complex **cal_error_term_vector;

    /* per-calibration properties */
    vnaproperty_t *cal_properties;

} vnacal_calibration_t;

/*
 * Hide the union.
 */
#define cal_z0		u.cal_z0
#define cal_z0_vector	u.cal_z0_vector
#define cal_z0_matrix	u.cal_z0_matrix

/*
 * vnacal_error_term_matrix_type_t: error term type used in load/save
 */
typedef enum vnacal_error_term_matrix_type {
    VETM_UNDEF,		/* invalid value */
    VETM_VECTOR,	/* vector (rows == 1) */
    VETM_MATRIX_ND,	/* matrix missing the major diagonal */
    VETM_MATRIX		/* rows x columns matrix */
} vnacal_error_term_matrix_type_t;

/*
 * vnacal_error_term_matrix_t: an error term vector/matrix used in load/save
 *
 * This structure is created by vnacal_build_error_term_list and used
 * in vnacal_load and vnacal_save.  The vetm_matrix element is a shallow
 * copy of elements from cal_error_term_vector.
 */
typedef struct vnacal_error_term_matrix {
    vnacal_calibration_t               *vetm_calp;
    vnacal_error_term_matrix_type_t	vetm_type;
    const char                         *vetm_name;
    double complex		      **vetm_matrix; /* matrix of vectors */
    int					vetm_rows;
    int					vetm_columns;
    struct vnacal_error_term_matrix    *vetm_next;
} vnacal_error_term_matrix_t;

/*
 * vnacal_t: structure returned from vnacal_create and vnacal_load
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

    /* calibration filename or NULL */
    char *vc_filename;

    /* precision for frequency values */
    int vc_fprecision;

    /* precision for data values */
    int vc_dprecision;

    /* global properties */
    vnaproperty_t *vc_properties;

    /* doubly linked ring list of associated vnacal_new_t structures */
    list_t vc_new_head;
};

/* _vnacal_error: report an error */
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
	vnacal_z0_type_t z0_type, int error_terms);

/* _vnacal_calibration_get_fmin_bound: get the lower frequency bound */
extern double _vnacal_calibration_get_fmin_bound(
	const vnacal_calibration_t *calp);

/* _vnacal_etermset_get_fmax_bound: get the upper frequency bound */
extern double _vnacal_calibration_get_fmax_bound(
	const vnacal_calibration_t *calp);

/* _vnacal_calibration_free: free the memory for a vnacal_calibration_t */
extern void _vnacal_calibration_free(vnacal_calibration_t *calp);

/* _vnacal_get_calibration: return the calibration at the given index */
extern vnacal_calibration_t *_vnacal_get_calibration(const char *function,
	const vnacal_t *vcp, int ci);

/* _vnacal_add_calibration_common */
extern int _vnacal_add_calibration_common(const char *function, vnacal_t *vcp,
	vnacal_calibration_t *calp, const char *name);

/* _vnacal_build_error_term_list: make list of ET matrices for load/save */
extern int _vnacal_build_error_term_list(vnacal_calibration_t *calp,
	const vnacal_layout_t *vlp, vnacal_error_term_matrix_t **head);

/* _vnacal_free_error_term_matrix: free an error term matrix structure */
extern void _vnacal_free_error_term_matrices(vnacal_error_term_matrix_t **head);

/* _vnacal_get_parameter: return a pointer to the parameter */
extern vnacal_parameter_t *_vnacal_get_parameter(const vnacal_t *vcp,
	int parameter);

/* _vnacal_alloc_parameter: allocate a vnacal_parameter and return index */
extern vnacal_parameter_t *_vnacal_alloc_parameter(const char *function,
	vnacal_t *vcp);

/* _vnacal_hold_parameter: increase the hold count on a parameter */
extern void _vnacal_hold_parameter(vnacal_parameter_t *vpmrp);

/* _vnacal_release_parameter: decrease the hold count and free if zero */
extern void _vnacal_release_parameter(vnacal_parameter_t *vpmrp);

/* _vnacal_get_calkit_name: return calkit name and number of ports */
extern const char *_vnacal_get_calkit_name(const vnacal_calkit_data_t *vcdp,
	int *ip_ports);

/*
 * SXX_BUFFER_ALLOC: size of buffer to hold "s%d_%d" safely
 */
#define SXX_BUFFER_ALLOC	(2 * (1 + 3 * sizeof(int)) + 3)

/*
 * MAX_DATA_STD_NAME: maximum length data standard name to show in error msg
 */
#define MAX_DATA_STD_NAME	31

/*
 * PARAMETER_BUFFER_ALLOC: buffer size to hold maximum length parameter name
 *   correlated parameter\0
 *   s22 of calkit through standard\0
 *   sNNN_NNN of "..............................." standard\0
 */
#define PARAMETER_BUFFER_ALLOC \
    (SXX_BUFFER_ALLOC + 3 + 1 + MAX_DATA_STD_NAME + 1 + 9 + 1)

/* _vnacal_format_sxx: write the name of an S parameter into buffer */
extern char *_vnacal_format_sxx(char *cur, char *end, int row, int column);

/* _vnacal_get_parameter_name: copy descriptive name for parameter */
extern void _vnacal_get_parameter_name(const vnacal_parameter_t *vpmrp,
	bool with_sxx, char *buffer);

/* _vnacal_alloc_standard: allocate a vnacal_standard_t structure */
extern void *_vnacal_alloc_standard(const char *function, vnacal_t *vcp,
	const vnacal_standard_ops_t *ops, int ports, size_t size);

/* _vnacal_release_standard: decrement the reference count on a standard */
extern void _vnacal_release_standard(vnacal_standard_t **stdpp);

/* _vnacal_get_parameter_frange: get the frequency limits for the parameter */
extern void _vnacal_get_parameter_frange(vnacal_parameter_t *vpmrp,
	double *fmin, double *fmax);

/* _vnacal_init_parameter_matrix: initialize parameter matrix to -1's */
static inline void _vnacal_init_parameter_matrix(int *parameter_matrix,
	int rows, int columns)
{
    for (int cell = 0; cell < rows * columns; ++cell) {
	parameter_matrix[cell] = -1;
    }
}

/* _vnacal_fill_standard_parameter_matrix: fill parameter matrix */
extern int _vnacal_fill_standard_parameter_matrix(const char *function,
	vnacal_standard_t *stdp, int *parameter_matrix);

/* _vnacal_eval_parameter_matrix: eval parameter matrix at given frequency */
extern int _vnacal_eval_parameter_matrix_i(const char *function,
	const vnacal_parameter_matrix_map_t *vpmmp, double frequency,
	const double complex *z0_vector, double complex *result_matrix);

/* _vnacal_get_correlated_sigma: return the sigma value for the given f */
extern double _vnacal_get_correlated_sigma(vnacal_parameter_t *vpmrp,
	double frequency);

/* _vnacal_analyze_parameter_matrix: check matrix and build standard maps */
extern vnacal_parameter_matrix_map_t *_vnacal_analyze_parameter_matrix(
	const char *function, vnacal_t *vcp, vnacal_parameter_t **matrix,
	int rows, int columns, bool initial);

/* _vnacal_free_parameter_matrix_map: free a vnacal_parameter_matrix_map_t */
extern void _vnacal_free_parameter_matrix_map(
	vnacal_parameter_matrix_map_t *vpmmp);

/* _vnacal_make_data_parameter_matrix: make a parameter matrix from data */
extern int _vnacal_make_data_parameter_matrix(const char *function,
	vnacal_t *vcp, const vnadata_t *vdp,
	int *parameter_matrix, size_t parameter_matrix_size);

/* _vnacal_embed_parameter_matrix: embed/de-embed a DUT into/from a fixture */
extern const vnacal_standard_ops_t _vnacal_embed_ops, _vnacal_deembed_ops;
extern int _vnacal_embed_parameter_matrix(const char *function,
	const vnacal_standard_ops_t *stdop, vnacal_t *vcp,
	const int *target_matrix, int target_ports, const int *fixture_matrix,
	int *result_matrix, size_t result_matrix_size);

/* _vnacal_setup_parameter_collection: allocate the parameter collection */
extern int _vnacal_setup_parameter_collection(const char *function,
	vnacal_t *vcp);

/* _vnacal_teardown_parameter_collection: free the parameter collection */
extern void _vnacal_teardown_parameter_collection(vnacal_t *vcp);

/* _vnacal_rfi: apply rational function interpolation */
extern double complex _vnacal_rfi(const double *xp, const double complex *yp,
	int n, int m, int *segment, double x);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_INTERNAL_H */
