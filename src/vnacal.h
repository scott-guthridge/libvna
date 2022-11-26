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
 * vnacal_new_set_z0: set the system impedance for all VNA ports
 *   @vnp: pointer to vnacal_new_t structure
 *   @z0: nominal impedance looking into a VNA port
 *
 * Note:
 *   In this implementation, all VNA ports must have the same system
 *   impedance.  If not set, the default is 50 ohms.
 */
extern int vnacal_new_set_z0(vnacal_new_t *vnp, double complex z0);

/*
 * vnacal_new_set_m_error: set VNA measurement error by frequency
 *   @vnp: pointer to vnacal_new_t structure
 *   @frequency_vector: vector of frequency points
 *   @frequencies: number of frequencies
 *   @noise_error_vector: vector of standard deviation of noise floor
 *   @tracking_error_vector: vector of standard deviation of tracking error
 */
extern int vnacal_new_set_m_error(vnacal_new_t *vnp,
	const double *frequency_vector, int frequencies,
	const double *noise_error_vector,
	const double *gain_error_vector);

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
 *   @limit: maximum number of iterations before failing
 */
extern int vnacal_new_set_iteration_limit(vnacal_new_t *vnp, int limit);

/*
 * vnacal_new_set_pvalue_limit: set the pvalue under which we reject the soln.
 *   @vnp: pointer to vnacal_new_t structure
 *   @limit: probability limit
 */
extern int vnacal_new_set_pvalue_limit(vnacal_new_t *vnp, double limit);

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
 *   @gamma: coefficient of reflection or transmission
 */
extern int vnacal_make_scalar_parameter(vnacal_t *vcp, double complex gamma);

/*
 * vnacal_make_vector_parameter: create frequency-dependent parameter
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @frequencies: length of frequency_vector and gamma_vector
 *   @frequency_vector: vector of increasing frequency values
 *   @gamma_vector: vector of per-frequency gamma values
 */
extern int vnacal_make_vector_parameter(vnacal_t *vcp,
	const double *frequency_vector, int frequencies,
	const double complex *gamma_vector);

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
 * vnacal_get_parameter_value: evaluate a parameter at a given frequency
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: index of parameter
 *   @frequency: frequency at which to evaluate parameter
 */
extern double complex vnacal_get_parameter_value(vnacal_t *vcp,
	int parameter, double frequency);

/*
 * vnacal_delete_parameter: delete the parameter with given index
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @parameter: parameter to delete
 */
extern int vnacal_delete_parameter(vnacal_t *vcp, int parameter);

/*
 * vnacal_create: construct a calibration structure from measured data
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
extern vnacal_t *vnacal_create(vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnacal_load: load the calibration from a file
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
 * vnacal_get_z0: return the system impedance for the given calibration
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @ci: calibration index
 */
extern double complex vnacal_get_z0(const vnacal_t *vcp, int ci);

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
