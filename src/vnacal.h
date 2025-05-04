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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _VNACAL_H
#define _VNACAL_H

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <vnadata.h>
#include <vnaerr.h>
#include <vnaproperty.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * VNACAL_MAX_PRECISION: argument to vnacal_set_fprecision and
 *	vnacal_set_dprecision for hexadecimal floating point
 */
#define VNACAL_MAX_PRECISION	1000

/*
 * Predefined Parameters
 */
#define VNACAL_MATCH	0	/* perfect reflections */
#define VNACAL_OPEN	1
#define VNACAL_SHORT	2
#define VNACAL_ZERO	0	/* aliases for transmission */
#define VNACAL_ONE	1

/*
 * vnacal_type_t: calibration type
 */
typedef enum vnacal_type {
    VNACAL_NOTYPE = -1, /* an invalid type value */
    VNACAL_T8,		/* 8-term T parameters */
    VNACAL_U8,		/* 8-term U (inverse T) parameters */
    VNACAL_TE10,	/* 8-term T plus off-diagonal E leakage terms */
    VNACAL_UE10,	/* 8-term U plus off-diagonal E leakage terms */
    VNACAL_T16,		/* 16-term T parameters */
    VNACAL_U16,		/* 16-term U (inverse T) parameters */
    VNACAL_UE14,	/* 14-term columns (rows x 1) U7 systems */
    _VNACAL_E12_UE14,	/* internal only -- used to compute E12 */
    VNACAL_E12,		/* 12-term generalized classic SOLT */
} vnacal_type_t;


/*
 * VNACAL_CK_MAGIC: magic number for validating vnacal_calkit_data_t
 */
#define VNACAL_CK_MAGIC	0x636b0000

/*
 * vnacal_calkit_type_t: type of parameterized calibration kit standard
 */
typedef enum vnacal_calkit_type {
    VNACAL_CALKIT_SHORT		= VNACAL_CK_MAGIC + 0,
    VNACAL_CALKIT_OPEN		= VNACAL_CK_MAGIC + 1,
    VNACAL_CALKIT_LOAD		= VNACAL_CK_MAGIC + 2,
    VNACAL_CALKIT_THROUGH	= VNACAL_CK_MAGIC + 3,
} vnacal_calkit_type_t;

/*
 * vnacal_calkit_data_t: parameters describing a calibration kit standard
 *
 * Note: we may add additional unions to this structure with defines
 * to hide them for source-level compatiblity, and we may extend the
 * struct to support new types of standards, but we must maintain binary
 * compatibilty with existing compiled code.
 */
typedef struct vnacal_calkit_data {
    vnacal_calkit_type_t vcd_type;
    uint32_t vcd_flags;         /* flags: see below */
    double vcd_offset_delay;	/* delay (s) */
    double vcd_offset_loss;	/* loss (Ω/s) */
    double vcd_offset_z0;	/* lossless characteristic Z (Ω) */
    double vcd_fmin;		/* minimum allowed frequency */
    double vcd_fmax;		/* maximum allowed frequency */
    union {
        double _vcd_coefficients[4];
        double complex _vcd_zl;
    } u;
} vnacal_calkit_data_t;
#define vcd_l_coefficients	u._vcd_coefficients
#define vcd_c_coefficients	u._vcd_coefficients
#define vcd_l0			u._vcd_coefficients[0]
#define vcd_l1			u._vcd_coefficients[1]
#define vcd_l2			u._vcd_coefficients[2]
#define vcd_l3			u._vcd_coefficients[3]
#define vcd_c0			u._vcd_coefficients[0]
#define vcd_c1			u._vcd_coefficients[1]
#define vcd_c2			u._vcd_coefficients[2]
#define vcd_c3			u._vcd_coefficients[3]
#define vcd_zl			u._vcd_zl

/*
 * Calkit Flags: vcd_flags is a bitwise OR of these flags
 *
 * VNACAL_CKF_TRADITIONAL:
 *     Use the traditional transmission line model described in Agilent
 *     note AN-1287-11 that uses an approximation to avoid the need for
 *     complex square root.  Otherwise, use the Keysight revised version.
 */
#define VNACAL_CKF_TRADITIONAL	0x0001U

/*
 * vnacal_t: opaque type returned from vnacal_load and vnacal_create
 */
typedef struct vnacal vnacal_t;

/*
 * vnacal_new_t: opaque type containing calibration measurements
 */
typedef struct vnacal_new vnacal_new_t;


/*
 * vnacal_name_to_type: convert error-term type name to enum
 *   @name: name of error term type (case insensitive)
 */
extern vnacal_type_t vnacal_name_to_type(const char *name);

/*
 * vnacal_type_to_name: convert error term type to name
 *   @type: type of error terms
 */
extern const char *vnacal_type_to_name(vnacal_type_t type);

/*
 * vnacal_new_alloc: allocate a vnacal_new_t structure
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @type: error term type
 *   @rows: number of VNA ports where signal is detected
 *   @columns: number of VNA ports where signal is generated
 *   @frequencies: number of frequency points
 *
 *   On error, set errno and return NULL.
 */
extern vnacal_new_t *vnacal_new_alloc(vnacal_t *vcp, vnacal_type_t type,
	int rows, int columns, int frequencies);

/*
 * vnacal_new_set_frequency_vector: set the frequency vector
 *   @vnp: pointer to vnacal_new_t structure
 *   @frequency_vector: vector of increasing frequencies
 */
extern int vnacal_new_set_frequency_vector(vnacal_new_t *vnp,
	const double *frequency_vector);

/*
 * vnacal_new_set_z0: set a common reference impedance for all VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @z0: reference impedance for all ports
 *
 * Note: if not set, the default is 50 ohms for all ports.
 */
static inline int vnacal_new_set_z0(vnacal_new_t *vnp, double complex z0)
{
    extern int _vnacal_new_set_z0_vector(const char *function,
	    vnacal_new_t *vnp, const double complex *z0_vector, int length);

    return _vnacal_new_set_z0_vector(__func__, vnp, &z0, 1);
}

/*
 * vnacal_new_set_z0_vector: set port/frequency specific reference impedances
 *   @vnp: pointer to vnacal_new_t structure
 *   @z0_vector: vector of z0 by port, or matrix of z0 by frequency, port
 *   @length: 1, #ports, or #frequencies
 *
 * Note: if not set, the default is 50 ohms for all ports.
 */
static inline int vnacal_new_set_z0_vector(vnacal_new_t *vnp,
	const double complex *z0_vector, int length)
{
    extern int _vnacal_new_set_z0_vector(const char *function,
	    vnacal_new_t *vnp, const double complex *z0_vector, int length);

    return _vnacal_new_set_z0_vector(__func__, vnp, z0_vector, length);
}

/*
 * vnacal_new_set_m_error: set VNA measurement error by frequency
 *   @vnp: pointer to vnacal_new_t structure
 *   @frequency_vector: vector of frequency points
 *   @frequencies: number of frequencies
 *   @sigma_nf_vector: standard deviations of the noise floor for measurements
 *   @sigma_tr_vector: standard deviation of noise proportional to signal level
 */
extern int vnacal_new_set_m_error(vnacal_new_t *vnp,
	const double *frequency_vector, int frequencies,
	const double *sigma_nf_vector, const double *sigma_tr_vector);

/*
 * vnacal_new_set_p_tolerance: set vnacal_new_solve iteration tolerance
 *   @vnp: pointer to vnacal_new_t structure
 *   @tolerance: new tolerance
 */
extern int vnacal_new_set_p_tolerance(vnacal_new_t *vnp,
	double tolerance);

/*
 * vnacal_new_set_et_tolerance: set the error term iteration tolerance
 *   @vnp: pointer to vnacal_new_t structure
 *   @tolerance: new tolerance
 */
extern int vnacal_new_set_et_tolerance(vnacal_new_t *vnp, double tolerance);

/*
 * vnacal_new_set_iteration_limit: set iteration limit for iterative solutions
 *   @vnp: pointer to vnacal_new_t structure
 *   @iterations: maximum number of iterations before failing
 */
extern int vnacal_new_set_iteration_limit(vnacal_new_t *vnp, int iterations);

/*
 * vnacal_new_set_pvalue_limit: set the pvalue under which we reject the soln.
 *   @vnp: pointer to vnacal_new_t structure
 *   @significance: reject if pvalue less than this value
 */
extern int vnacal_new_set_pvalue_limit(vnacal_new_t *vnp, double significance);

/*
 * vnacal_new_add_single_reflect: add a single reflect on the given port
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s11: s11 parameter index, e.g. VNACAL_MATCH
 *   @port: VNA port on which the measurement is made
 */
extern int vnacal_new_add_single_reflect(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	int s11, int port);

/*
 * vnacal_new_add_single_reflect_m: add a single reflect on the given port
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s11: s11 parameter index, e.g. VNACAL_MATCH
 *   @port: VNA port on which the measurement is made
 */
extern int vnacal_new_add_single_reflect_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	int s11, int port);

/*
 * vnacal_new_add_double_reflect: add a pair of reflects
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s11: s11 parameter index
 *   @s22: s22 parameter index
 *   @port1: first reflect is on this VNA port
 *   @port2: second reflect is on this VNA port
 */
extern int vnacal_new_add_double_reflect(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	int s11, int s22, int port1, int port2);

/*
 * vnacal_new_add_double_reflect_m: add a pair of reflects
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s11: s11 parameter index, e.g. VNACAL_MATCH
 *   @port1: first reflect is on this VNA port
 *   @port2: second reflect is on this VNA port
 */
extern int vnacal_new_add_double_reflect_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	int s11, int s22, int port1, int port2);

/*
 * vnacal_new_add_line: add an arbitrary two-port standard
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s_2x2: 2x2 matrix of parameter indices (by rows) describing the standard
 *   @port1: first VNA port attached to standard
 *   @port2: second VNA port attached to standard
 */
extern int vnacal_new_add_line(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	const int *s_2x2, int port1, int port2);

/*
 * vnacal_new_add_line_m: add an arbitrary two-port standard
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s_2x2: 2x2 matrix of parameter indices (by rows) describing the standard
 *   @port1: first VNA port attached to standard
 *   @port2: second VNA port attached to standard
 */
extern int vnacal_new_add_line_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	const int *s_2x2, int port1, int port2);

/*
 * vnacal_new_add_through: add a perfect through between two VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @port1: first VNA port attached to through
 *   @port2: second VNA port attached to through
 */
extern int vnacal_new_add_through(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	int port1, int port2);

/*
 * vnacal_new_add_through_m: add a perfect through between two VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @port1: first VNA port attached to through
 *   @port2: second VNA port attached to through
 */
extern int vnacal_new_add_through_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	int port1, int port2);

/*
 * vnacal_new_add_mapped_matrix: add a matrix of measurements with port map
 *   @vnp: pointer to vnacal_new_t structure
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s: matrix of parameter indices (by rows)
 *   @s_rows: number of rows in s
 *   @s_columns: number of columns in s
 *   @port_map: vector of VNA port numbers corresponding to these ports
 */
extern int vnacal_new_add_mapped_matrix(vnacal_new_t *vnp,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	const int *s, int s_rows, int s_columns,
	const int *port_map);

/*
 * vnacal_new_add_mapped_matrix_m: add a matrix of measurements with port map
 *   @vnp: pointer to vnacal_new_t structure
 *   @m: matrix of measured voltages entering the VNA
 *   @m_rows: number of rows in m
 *   @m_columns: number of columns in m
 *   @s: matrix of parameter indices (by rows)
 *   @s_rows: number of rows in s
 *   @s_columns: number of columns in s
 *   @port_map: vector of VNA port numbers corresponding to these ports
 */
extern int vnacal_new_add_mapped_matrix_m(vnacal_new_t *vnp,
	double complex *const *m, int m_rows, int m_columns,
	const int *s, int s_rows, int s_columns,
	const int *port_map);

/*
 * vnacal_new_solve: solve for the error parameters
 *   @vnp: pointer to vnacal_new_t structure
 */
extern int vnacal_new_solve(vnacal_new_t *vnp);

/*
 * vnacal_new_free: free a vnacal_new_t structure
 *   @vnp: pointer to vnacal_new_t structure
 */
extern void vnacal_new_free(vnacal_new_t *vnp);

/*
 * vnacal_make_scalar_parameter: create frequency-independent parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @coefficient: coefficient of reflection or transmission
 */
extern int vnacal_make_scalar_parameter(vnacal_t *vcp,
	double complex coefficient);

/*
 * vnacal_make_vector_parameter: create frequency-dependent parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @frequencies: length of frequency_vector and coefficient_vector
 *   @frequency_vector: vector of increasing frequency values
 *   @coefficient_vector: vector of per-frequency parameter values
 */
extern int vnacal_make_vector_parameter(vnacal_t *vcp,
	const double *frequency_vector, int frequencies,
	const double complex *coefficient_vector);

/*
 * vnacal_make_unknown_parameter: create an unknown parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @initial_guess: index of scalar or vector parameter serving as guess
 */
extern int vnacal_make_unknown_parameter(vnacal_t *vcp, int initial_guess);

/*
 * vnacal_make_correlated_parameter: create unknown parameter related by sigma
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @other: another parameter close to this one
 *   @sigma_frequency_vector: vector of increasing frequency values
 *   @sigma_frequencies: length of frequency_vector and sigma_vector
 *   @sigma_vector: frequency dependent parameter deviation from other
 */
extern int vnacal_make_correlated_parameter(vnacal_t *vcp, int other,
	const double *sigma_frequency_vector, int sigma_frequencies,
	const double *sigma_vector);

/*
 * vnacal_make_calkit_parameter: make a parameter for a one-port kit standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vcdp: data describing the calkit standard
 */
extern int vnacal_make_calkit_parameter(vnacal_t *vcp,
	const vnacal_calkit_data_t *vcdp);

/*
 * vnacal_make_calkit_parameter_matrix: make parameter matrix for kit standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vcdp: data describing the calkit standard
 *   @parameter_matrix: caller-supplied result matrix
 *   @parameter_matrix_size: size in bytes of the result matrix
 *
 * Fill parameter_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_single_reflect* or vnadata_new_add_line* functions.
 * The parameter_matrix_size parameter is the allocation in bytes of the
 * result matrix, used to protect against buffer overrun.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
extern int vnacal_make_calkit_parameter_matrix(vnacal_t *vcp,
	const vnacal_calkit_data_t *vcdp, int *parameter_matrix,
	size_t parameter_matrix_size);

/*
 * vnacal_make_data_parameter: make a parameter from 1x1 network parameter data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vdp: network parameter data for a calibration standard
 *
 * Dimensions of vdp must be 1x1. Data must be convertable to
 * S-parameters.  Automatically handles parameter
 * conversion, interpolation and renormalization as needed.
 */
extern int vnacal_make_data_parameter(vnacal_t *vcp,
	const vnadata_t *vdp);

/*
 * vnacal_make_data_parameter_matrix: make a parameter matrix from data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vdp: network parameter data for a calibration standard
 *   @parameter_matrix: caller-allocated matrix to receive result
 *   @parameter_matrix_size: size in bytes of the result matrix
 *
 * Fill parameter_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_* functions.  Automatically handles parameter
 * conversion, interpolation and renormalization.  Data must be
 * convertable to S-parameters.  The parameter_matrix_size parameter is
 * the allocation in bytes of the result matrix, used to protect against
 * buffer overrun.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
extern int vnacal_make_data_parameter_matrix(vnacal_t *vcp,
	const vnadata_t *vdp, int *parameter_matrix,
	size_t parameter_matrix_size);

/*
 * vnacal_get_parameter_value: evaluate a parameter at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 *   @frequency: frequency at which to evaluate the parameter
 *
 * Note: this function works only on scalar, vector, unknown and
 * correlated parameters where the reference impedance is implicit.
 * For calkit and data parameters, use the more general functions
 * vnacal_eval_parameter, vnacal_eval_parameter_matrix and
 * vnacal_parameter_matrix_to_data.
 */
extern double complex vnacal_get_parameter_value(vnacal_t *vcp,
	int parameter, double frequency);

/*
 * vnacal_eval_parameter: evaluate a parameter at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 *   @frequency: frequency at which to evaluate the parameter
 *   @z0: reference impedance for the result
 */
extern double complex vnacal_eval_parameter(vnacal_t *vcp, int parameter,
	double frequency, double complex z0);

/*
 * vnacal_eval_parameter_matrix: evaluate parameter matrix at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter_matrix: index of parameter
 *   @rows: number of rows in parameter matrix
 *   @columns: number of columns in parameter matrix
 *   @frequency: frequency at which to evaluate the parameter
 *   @z0_vector: reference impedance for each port of the result
 *   @result_matrix: caller-supplied matrix to hold the result
 */
extern int vnacal_eval_parameter_matrix(vnacal_t *vcp,
	const int *parameter_matrix, int rows, int columns, double frequency,
	const double complex *z0_vector, double complex *result_matrix);

/*
 * vnacal_parameter_matrix_to_data: convert parameter matrix to parameter data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter_matrix: index of parameter
 *   @rows: number of rows in parameter matrix
 *   @columns: number of columns in parameter matrix
 *   @vdp: takes frequency vector, z0 and type as input; returns data
 */
extern int vnacal_parameter_matrix_to_data(vnacal_t *vcp,
	const int *parameter_matrix, int rows, int columns, vnadata_t *vdp);

/*
 * vnacal_delete_parameter: delete the parameter with given index
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: parameter to delete
 */
extern int vnacal_delete_parameter(vnacal_t *vcp, int parameter);

/*
 * vnacal_delete_parameter_matrix: delete the parameters in the given matrix
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter_matrix: matrix of parameter indices
 *   @rows: rows in parameter_matrix
 *   @columns: columns in parameter_matrix
 *
 * Note: this function does not free the parameter matrix itself.  Use
 * free to return the memory.
 */
extern void vnacal_delete_parameter_matrix(vnacal_t *vcp,
	const int *parameter_matrix, int rows, int columns);

/*
 * vnacal_create: create the main structure for a new calibration
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
extern vnacal_t *vnacal_create(vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnacal_load: load an existing calibration from a file
 *   @pathname: calibration file name
 *   @error_fn: error reporting callback or NULL
 *   @error_arg: arbitrary argument passed through to error_fn or NULL
 *
 *   If error_fn is non-NULL, then vnacal_load and subsequent functions report
 *   error messages using error_fn before returning failure to the caller.
 */
extern vnacal_t *vnacal_load(const char *pathname,
	vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnacal_save: create or overwrite a calibration file with new data
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @pathname: calibration file name
 *
 *   The pathname and basename parameters work as in vnacal_load except
 *   that the $HOME/{pathname} directory path is created if necessary.
 */
extern int vnacal_save(vnacal_t *vcp, const char *pathname);

/*
 * vnacal_get_filename: return the calibration file name, if established
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
extern const char *vnacal_get_filename(const vnacal_t *vcp);

/*
 * vnacal_add_calibration: add a new calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @name: name of new calibration
 *   @vnp: pointer to vnacal_new_t structure
 */
extern int vnacal_add_calibration(vnacal_t *vcp, const char *name,
	vnacal_new_t *vnp);

/*
 * vnacal_find_calibration: find a calibration by name
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @name: name of calibration to find
 */
extern int vnacal_find_calibration(const vnacal_t *vcp, const char *name);

/*
 * vnacal_delete_calibration: delete a calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern int vnacal_delete_calibration(vnacal_t *vcp, int ci);

/*
 * vnacal_get_calibration_end: return one past the highest calibration index
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
extern int vnacal_get_calibration_end(const vnacal_t *vcp);

/*
 * vnacal_get_name: return the name of the given calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern const char *vnacal_get_name(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_type: return the type of error terms
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern vnacal_type_t vnacal_get_type(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_rows: return the number of rows in the calibration matrix
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern int vnacal_get_rows(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_columns: return the number of columns in the calibration matrix
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern int vnacal_get_columns(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_frequencies: return the number of frequency points
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern int vnacal_get_frequencies(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_fmin: return the minimum frequency point
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern double vnacal_get_fmin(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_fmax: return the maximum frequency point
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern double vnacal_get_fmax(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_frequency_vector: return a pointer to the calibration frequencies
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern const double *vnacal_get_frequency_vector(const vnacal_t *vcp,
	int ci);

/*
 * vnacal_z0_type_t: type of reference impedances in a calibration
 */
typedef enum vnacal_z0_type {
    VNACAL_Z0_INVALID = -1,
    VNACAL_Z0_SCALAR,	/* all ports are referenced to the same impedance */
    VNACAL_Z0_VECTOR,	/* VNA ports have per-port reference impedances */
    VNACAL_Z0_MATRIX	/* reference impedances vary by port and frequency */
} vnacal_z0_type_t;

/*
 * vnacal_get_z0_type: return the reference impedance type
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *
 * Returns:
 *   VNACAL_Z0_SCALAR  if all ports have the same reference impedance
 *   VNACAL_Z0_VECTOR  if reference impedances vary by port
 *   VNACAL_Z0_MATRIX  if reference impedances vary by port & frequency
 *   -1                if the arguments are invalid
 */
extern vnacal_z0_type_t vnacal_get_z0_type(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_z0: return the reference impedance for the given calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *
 *   This function can be used only when the z0 type is VNACAL_Z0_SCALAR,
 *   i.e. the reference impedances of all ports are the same.
 */
extern double complex vnacal_get_z0(const vnacal_t *vcp, int ci);

/*
 * vnacal_get_z0_vector: return a calibration reference impedance vector
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @vector: caller-provided #ports-long buffer to receive result
 *   @max_entries: number of double complex entries in vector
 *   @f: frequency at which to evaluate
 *
 *   Copies #ports reference impedances into the caller-provided buffer.
 *   The buffer should have space for at least one double complex entry
 *   per VNA port.  The max_entries parameter gives the length of the
 *   provided buffer to guard against buffer overrun.
 *
 *   When the z0 type is VNACAL_Z0_SCALAR, the single reference impedance
 *   is duplicated for each port.  When it's VNACAL_Z0_VECTOR, the entries
 *   are copied into the user's buffer.  When it's VNACAL_Z0_MATRIX,
 *   then the function returns the reference impedances for the given
 *   frequency, interpolating if necessary.  The frequency argument is
 *   ignored if the z0 type is not VNACAL_Z0_MATRIX.
 *
 * Return:
 *   number of VNA ports (number of entries placed into vector), or
 *   -1 on error
 */
extern int vnacal_get_z0_vector(const vnacal_t *vcp, int ci,
	double complex *vector, int max_entries, double f);

/*
 * vnacal_set_fprecision: set the frequency value precision for vnacal_save
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @precision: precision in decimal places (1..n) or VNACAL_MAX_PRECISION
 */
extern int vnacal_set_fprecision(vnacal_t *vcp, int precision);

/*
 * vnacal_set_dprecision: set the data value precision for vnacal_save
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @precision: precision in decimal places (1..n) or VNACAL_MAX_PRECISION
 */
extern int vnacal_set_dprecision(vnacal_t *vcp, int precision);

/*
 * vnacal_property_type: get the type of the given property expression
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_type(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
    ;

/*
 * vnacal_property_count: return count of elements in given collection
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_count(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
    ;

/*
 * vnacal_property_keys: return a vector of keys for the given map expr
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
extern const char **vnacal_property_keys(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
    ;

/*
 * vnacal_property_get: get a property value from a property expression
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
extern const char *vnacal_property_get(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
    ;

/*
 * vnacal_property_set: set a property value from a property expression
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_set(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
    ;

/*
 * vnacal_property_delete: delete the value described by format
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format: printf format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_delete(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
    ;

/*
 * vnacal_property_get_subtree: get the subtree described by format
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
extern vnaproperty_t *vnacal_property_get_subtree(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
;

/*
 * vnacal_property_set_subtree: for subtree and return address
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @format:  printf-like format string forming the property expression
 *   @...:     optional variable arguments
 */
extern vnaproperty_t **vnacal_property_set_subtree(vnacal_t *vcp, int ci,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)))
#endif
;

/*
 * vnacal_free: free a vnacal_t structure
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 */
extern void vnacal_free(vnacal_t *vcp);

/*
 * vnacal_apply: apply the calibration to measured values
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @frequency_vector: vector of increasing frequency points
 *   @frequencies: number of frequencies in matrix and s_parameters
 *   @a: matrix of measured voltages leaving the VNA
 *   @a_rows: number of rows in a
 *   @a_columns: number of columns in a
 *   @b: matrix of measured voltages entering the VNA
 *   @b_rows: number of rows in b
 *   @b_columns: number of columns in b
 *   @s_parameters: caller-allocated vnadata_t structure
 */
extern int vnacal_apply(vnacal_t *vcp, int ci,
	const double *frequency_vector, int frequencies,
	double complex *const *a, int a_rows, int a_columns,
	double complex *const *b, int b_rows, int b_columns,
	vnadata_t *s_parameters);

/*
 * vnacal_apply_m: apply the calibration to measured values
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 *   @frequency_vector: vector of increasing frequency points
 *   @frequencies: number of frequencies in matrix and s_parameters
 *   @m: matrix of measured voltages
 *   @m_rows: number of rows in b
 *   @m_columns: number of columns in b
 *   @s_parameters: caller-allocated vnadata_t structure
 */
extern int vnacal_apply_m(vnacal_t *vcp, int ci,
	const double *frequency_vector, int frequencies,
	double complex *const *m, int m_rows, int m_columns,
	vnadata_t *s_parameters);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_H */
