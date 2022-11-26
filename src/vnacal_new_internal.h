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
#define VNACAL_NEW_DEFAULT_ET_TOLERANCE		1.0e-6
#define VNACAL_NEW_DEFAULT_ITERATION_LIMIT	30
#define VNACAL_NEW_DEFAULT_PVALUE_LIMIT		0.001

/*
 * vnacal_new_m_error_t: measurement error
 */
typedef struct vnacal_new_m_error {
    /* standard deviation of measurement noise */
    double vnme_noise;

    /* standard deviation of measurement gain error */
    double vnme_gain;

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
 * vnacal_new_term_t: an expanded term of an equation
 */
typedef struct vnacal_new_term {
    /* index of associated unknown, or -1 for right hand side "b" vector */
    int vnt_xindex;

    /* multiply by -1 */
    bool vnt_negative;

    /* index into vnm_m_matrix or -1 if no m */
    int vnt_m_cell;

    /* index into vnm_s_matrix or -1 if no s */
    int vnt_s_cell;

    /* index into vnsm_v_matrix or -1 if no v */
    int vnt_v_cell;

    /* thread through the diagonal v_matrix coefficients */
    struct vnacal_new_term *vnt_next_no_v;

    /* next term in equation */
    struct vnacal_new_term *vnt_next;

} vnacal_new_term_t;

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

    /* subset of terms for not using the v matrix */
    vnacal_new_term_t *vne_term_list_no_v;

    /* list of all expanded terms of the equation multiplied by the v matrix */
    vnacal_new_term_t *vne_term_list;

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
    bool *vnm_connectivity_matrix;

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

    /* iterative solve not satified until RMS change in p <= than this */
    double vn_p_tolerance;

    /* iterative solve not satified until RMS chagne in et <= than this */
    double vn_et_tolerance;

    /* maximum number of iterations allowed in iterative solutions */
    int vn_iteration_limit;

    /* vnacal_new_solve pvalue lower than this threshold considered failing */
    double vn_pvalue_limit;

    /* hidden API for test: optional caller supplied vector for pvalue per f. */
    double *vn_pvalue_vector;

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

/*
 * vnacal_new_leakage_term_t: leakage term outside the linear system
 */
typedef struct vnacal_new_leakage_term {
    /* sum of the samples */
    double complex vnlt_sum;

    /* sum of squared magnitudes of the samples */
    double vnlt_sumsq;

    /* count of accumulated samples */
    int vnlt_count;

} vnacal_new_leakage_term_t;

/*
 * vnacal_new_iterator_state_t: equation iterator states
 */
typedef enum {
    VNACAL_NI_INIT,			/* not started */
    VNACAL_NI_SYSTEM,			/* in system */
    VNACAL_NI_EQUATION,			/* in equation */
    VNACAL_NI_TERM,			/* in term list */
    VNACAL_NI_END_TERMS,		/* no remaining terms */
    VNACAL_NI_END_EQUATIONS		/* no remaining equations */
} vnacal_new_iterator_state_t;

/*
 * vnacal_new_msv_matrices_t: a measured standard for solve
 */
typedef struct vnacal_new_msv_matrices {
    /* correponding measured standard */
    struct vnacal_new_measurement *vnmm_vnmp;

    /* matrix of measured values for the current frequency */
    double complex *vnmm_m_matrix;

    /* matrix of values of the standard for the current frequency */
    double complex *vnmm_s_matrix;

    /* vector of matrices (one per system) that transform the residuals
       of the equations into units of the measurement primary associated
       with the equation. */
    double complex **vnsm_v_matrices;

} vnacal_new_msv_matrices_t;

/*
 * vnacal_new_trl_indices_t: indices of simple TRL standards and unknowns
 *   This structure is used only when it's possible to solve TRL using
 *   the simple, analytical method.
 */
typedef struct vnacal_new_trl_indices {
    /* index of the measured T standard */
    int vnti_t_standard;

    /* index of the measured R standard */
    int vnti_r_standard;

    /* index of the measured L standard */
    int vnti_l_standard;

    /* index of the R unknown */
    int vnti_r_unknown;

    /* index of the L unknown */
    int vnti_l_unknown;

} vnacal_new_trl_indices_t;

/*
 * vnacal_new_solve_state: iterate over vnacal_new_term_t's
 */
typedef struct vnacal_new_solve_state {
    /* new calibration structure */
    vnacal_new_t *vnss_vnp;

    /* current frequency index */
    int vnss_findex;

    /* vector of structures correponding to each measured standard */
    vnacal_new_msv_matrices_t *vnss_msv_matrices;

    /* serialized matrix of pointers to leakage term structures */
    vnacal_new_leakage_term_t **vnss_leakage_matrix;

    /* vector of vector of unknown parameter values [index][findex] */
    double complex **vnss_p_vector;

    /* equation iterator state */
    vnacal_new_iterator_state_t vnss_iterator_state;

    /* include coefficients using the v matrix */
    bool vnss_include_v;

    /* current system in iterator */
    int vnss_sindex;

    /* current equation in iterator */
    vnacal_new_equation_t *vnss_vnep;

    /* current term in iterator */
    vnacal_new_term_t *vnss_vntp;

} vnacal_new_solve_state_t;

#define vs_init				_vnacal_new_solve_init
#define vs_start_frequency		_vnacal_new_solve_start_frequency
#define vs_next_equation		_vnacal_new_solve_next_equation
#define vs_next_term			_vnacal_new_solve_next_term
#define vs_update_s_matrices		_vnacal_new_solve_update_s_matrices
#define vs_update_v_matrices		_vnacal_new_solve_update_v_matrices
#define vs_update_all_v_matrices	_vnacal_new_solve_update_all_v_matrices
#define vs_calc_weights			_vnacal_new_solve_calc_weights
#define vs_calc_pvalue                  _vnacal_new_solve_calc_pvalue
#define vs_free				_vnacal_new_solve_free

/*
 * vs_start_system: prepare equation iterator for new system
 *   @vnssp: solve state structure
 *   @vnsp:  vnacal_new_system_t structure
 */
static inline void vs_start_system(vnacal_new_solve_state_t *vnssp, int sindex)
{
    vnssp->vnss_iterator_state = VNACAL_NI_SYSTEM;
    vnssp->vnss_sindex = sindex;
    vnssp->vnss_vnep   = NULL;
    vnssp->vnss_vntp   = NULL;
}

/*
 * vs_get_xindex: return the index of the unknown or -1 for right hand side
 *   @vnssp: solve state structure
 */
static inline int vs_get_xindex(const vnacal_new_solve_state_t *vnssp)
{
    assert(vnssp->vnss_iterator_state == VNACAL_NI_TERM);
    return vnssp->vnss_vntp->vnt_xindex;
}

/*
 * vs_get_negative: test if the current term has a minus sign
 *   @vnssp: solve state structure
 */
static inline bool vs_get_negative(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vntp->vnt_negative;
}

/*
 * vs_have_m: test if the current term has an m factor
 *   @vnssp: solve state structure
 */
static inline bool vs_have_m(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vntp->vnt_m_cell >= 0;
}

/*
 * vs_get_m: return the m value for the current term
 *   @vnssp: solve state structure
 */
static inline double complex vs_get_m(const vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_term_t *vntp = vnssp->vnss_vntp;
    int m_cell;
    vnacal_new_measurement_t *vnmp;
    vnacal_new_msv_matrices_t *vnmmp;

    assert(vnssp->vnss_iterator_state == VNACAL_NI_TERM);
    m_cell = vntp->vnt_m_cell;
    assert(m_cell >= 0);
    vnmp = vnssp->vnss_vnep->vne_vnmp;
    assert(vnmp->vnm_m_matrix[m_cell] != NULL);
    vnmmp = &vnssp->vnss_msv_matrices[vnmp->vnm_index];

    return vnmmp->vnmm_m_matrix[m_cell];
}

/*
 * vs_get_m_cell: get the index in the m matrix for the current term
 *   @vnssp: solve state structure
 */
static inline int vs_get_m_cell(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vntp->vnt_m_cell;
}

/*
 * vs_have_s: test if the current term has an s factor
 *   @vnssp: solve state structure
 */
static inline bool vs_have_s(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vntp->vnt_s_cell >= 0;
}

/*
 * vs_get_s: return the s value for the current term
 *   @vnssp: solve state structure
 */
static inline double complex vs_get_s(const vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_term_t *vntp = vnssp->vnss_vntp;
    int s_cell;
    vnacal_new_measurement_t *vnmp;
    vnacal_new_msv_matrices_t *vnmmp;

    assert(vnssp->vnss_iterator_state == VNACAL_NI_TERM);
    s_cell = vntp->vnt_s_cell;
    assert(s_cell >= 0);
    vnmp = vnssp->vnss_vnep->vne_vnmp;
    assert(vnmp->vnm_s_matrix[s_cell] != NULL);
    vnmmp = &vnssp->vnss_msv_matrices[vnmp->vnm_index];

    return vnmmp->vnmm_s_matrix[s_cell];
}

/*
 * vs_get_s_cell: get the index in the s matrix for the current term
 *   @vnssp: solve state structure
 */
static inline int vs_get_s_cell(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_vntp->vnt_s_cell;
}

/*
 * vs_have_v: test if the current system uses v matrices
 *   @vnssp: solve state structure
 */
static inline bool vs_have_v(const vnacal_new_solve_state_t *vnssp)
{
    return vnssp->vnss_include_v;
}

/*
 * vs_get_v: return the v value for the current coefficient
 *   @vnssp: solve state structure
 */
static inline double complex vs_get_v(vnacal_new_solve_state_t *vnssp)
{
    vnacal_new_term_t *vntp = vnssp->vnss_vntp;
    int v_cell;
    vnacal_new_measurement_t *vnmp;
    vnacal_new_msv_matrices_t *vnmmp;

    assert(vnssp->vnss_iterator_state == VNACAL_NI_TERM);
    v_cell = vntp->vnt_v_cell;
    assert(v_cell >= 0);
    if (!vnssp->vnss_include_v) {
       return 1.0;
    }
    vnmp = vnssp->vnss_vnep->vne_vnmp;
    assert(vnmp->vnm_s_matrix[v_cell] != NULL);
    vnmmp = &vnssp->vnss_msv_matrices[vnmp->vnm_index];
    assert(vnmmp->vnsm_v_matrices != NULL);
    assert(vnmmp->vnsm_v_matrices[vnssp->vnss_sindex] != NULL);

    return vnmmp->vnsm_v_matrices[vnssp->vnss_sindex][v_cell];
}

/* _vnacal_new_get_parameter: add/find parameter and return held */
extern vnacal_new_parameter_t *_vnacal_new_get_parameter(
	const char *function, vnacal_new_t *vnp, int parameter);

/* _vnacal_new_init_parameter_hash: set up the parameter hash */
extern int _vnacal_new_init_parameter_hash(const char *function,
	vnacal_new_parameter_hash_t *vnphp);

/* _vnacal_new_free_parameter_hash: free the parameter hash */
extern void _vnacal_new_free_parameter_hash(
	vnacal_new_parameter_hash_t *vnphp);

/* _vnacal_new_build_equation_terms: build vnacal_new_term_t lists */
extern int _vnacal_new_build_equation_terms(vnacal_new_equation_t *vnep);

/* free a vnacal_new_measurement_t structure */
extern void _vnacal_new_free_measurement(vnacal_new_measurement_t *vnmp);

/* _vnacal_new_err_need_full_s: report error for incomplete S with M errors */
extern void _vnacal_new_err_need_full_s(const vnacal_new_t *vnp,
	const char *function, int measurement, int s_cell);

/* _vnacal_new_check_all_frequency_ranges: check f range of all parameters */
extern int _vnacal_new_check_all_frequency_ranges(const char *function,
	vnacal_new_t *vnp, double fmin, double fmax);

/* _vnacal_new_solve_init: initialize the solve state structure */
extern int _vnacal_new_solve_init(vnacal_new_solve_state_t *vnssp,
       vnacal_new_t *vnp);

/* _vnacal_new_solve_start_frequency: start a new frequency */
extern int _vnacal_new_solve_start_frequency(vnacal_new_solve_state_t *vnssp,
	int findex);

/* _vnacal_new_solve_next_equation: move to the next equation in the system */
extern bool _vnacal_new_solve_next_equation(vnacal_new_solve_state_t *vnssp);

/* _vnacal_new_solve_next_term: move to the next term */
extern bool _vnacal_new_solve_next_term(vnacal_new_solve_state_t *vnssp);

/* _vnacal_new_init_x_vector: initialize the error terms to perfect values */
extern void _vnacal_new_solve_init_x_vector(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length);

/* _vnacal_new_solve_update_s_matrices: update s matrix unknown parameters */
extern void _vnacal_new_solve_update_s_matrices(
	vnacal_new_solve_state_t *vnssp);

/* _vnacal_new_solve_update_v_matrices: update vnsm_v_matrices for 1 system */
extern int _vnacal_new_solve_update_v_matrices(const char *function,
       vnacal_new_solve_state_t *vnssp, int sindex,
       const double complex *x_vector, int x_length);

/* _vnacal_new_solve_update_all_v_matrices: update v_matrices for all systems */
extern int _vnacal_new_solve_update_all_v_matrices(const char *function,
       vnacal_new_solve_state_t *vnssp, const double complex *x_vector,
       int x_length);

/* _vnacal_new_solve_calc_weights: calculate weights from measurement errors */
extern double *_vnacal_new_solve_calc_weights(vnacal_new_solve_state_t *vnssp);

/* _vnacal_new_solve_simple: solve error terms when all s-parameters known */
extern int _vnacal_new_solve_simple(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length);

/* _vnacal_new_solve_auto: solve error terms and unknown s-parameters */
extern int _vnacal_new_solve_auto(vnacal_new_solve_state_t *vnssp,
	double complex *x_vector, int x_length);

/* _vnacal_new_solve_calc_pvalue: calc p-value that system is consistent */
extern double _vnacal_new_solve_calc_pvalue(vnacal_new_solve_state_t *vnssp,
       const double complex *x_vector, int x_length);

/* _vnacal_new_solve_free: free resources held by the solve state structure */
extern void _vnacal_new_solve_free(vnacal_new_solve_state_t *vnssp);

/* _vnacal_new_solve_is_trl: test if we can solve using simple TRL */
extern bool _vnacal_new_solve_is_trl(const vnacal_new_t *vnp,
	vnacal_new_trl_indices_t *vntip);

/* _vnacal_new_solve_trl: solve when all s-parameters are known */
extern int _vnacal_new_solve_trl(vnacal_new_solve_state_t *vnssp,
	const vnacal_new_trl_indices_t *vntip,
	double complex *x_vector, int x_length);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_NEW_INTERNAL_H */
