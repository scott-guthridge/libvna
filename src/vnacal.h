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

#ifndef _VNACAL_H
#define _VNACAL_H

#include <complex.h>
#include <stdbool.h>
#include <stdint.h>
#include <vnadata.h>
#include <vnaerr.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * VNACAL_MAX_PRECISION: argument to vnacal_set_fprecision and
 *      vnacal_set_dprecision for hexadecimal floating point
 */
#define VNACAL_MAX_PRECISION	1000

/*
 * vnacal_type_t: calibration type
 */
typedef enum vnacal_type {
    VNACAL_E12	= 1,	/* traditional "12" error term model */
} vnacal_type_t;

/*
 * vnacal_calset_t: opaque type containing calibration measurements
 */
typedef struct vnacal_calset vnacal_calset_t;

/*
 * vnacal_t: opaque type returned from vnacal_load and vnacal_create
 */
typedef struct vnacal vnacal_t;

/*
 * vnacal_apply_t: opaque type returned from vnacal_apply_alloc
 */
typedef struct vnacal_apply vnacal_apply_t;

/*
 * Optional parameter tag used to sanity check that diagnonal and
 * off-diagonal calibration parameters aren't mixed.  If the tag
 * is isn't given then a parameter in 0-2 is assumed correct.
 */
#define _VNACAL_DIAGONAL	0x40
#define _VNACAL_OFF_DIAGONAL	0x80

/*
 * Calibration Term
 */
#define VNACAL_Sii_REF0		(_VNACAL_DIAGONAL     | 0)
#define VNACAL_Sii_REF1		(_VNACAL_DIAGONAL     | 1)
#define VNACAL_Sii_REF2		(_VNACAL_DIAGONAL     | 2)
#define VNACAL_Sjj_THROUGH	(_VNACAL_OFF_DIAGONAL | 0)
#define VNACAL_Sij_THROUGH	(_VNACAL_OFF_DIAGONAL | 1)
#define VNACAL_Sij_LEAKAGE	(_VNACAL_OFF_DIAGONAL | 2)

/*
 * vnacal_calset_alloc: allocate a vnacal_calset
 *   @type: calibration type
 *   @setname: name of this set
 *   @rows: number of VNA ports where signal is detected
 *   @columns: number of VNA ports where signal is generated
 *   @frequencies: number of frequency points
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 *
 *   Allocate a vnacal_calset_t and contained vcs_frequency_vector,
 *   vcs_matrix, vcd_sii_reference_vector, vcd_sii_through_vector,
 *   vcd_sji_through_vector and vcd_sdj_leakage vectors.
 *
 *   Set errno and return NULL on error.
 */
extern vnacal_calset_t *vnacal_calset_alloc(vnacal_type_t type,
	const char *setname, int rows, int columns, int frequencies,
	vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnacal_calset_set_frequency_vector: set the frequency vector
 *   @vcsp: pointer to vnacal_calset_t
 *   @frequency_vector: vector of increasing frequencies
 */
extern int vnacal_calset_set_frequency_vector(vnacal_calset_t *vcsp,
	const double *frequency_vector);

/*
 * vnacal_calset_set_z0: set the system impedance for all VNA ports
 *   @vcsp: pointer to vnacal_calset_t
 *   @z0:   nominal impedance looking into a VNA port
 *
 * Note:
 *   We currently assume all VNA ports have the same system impedance.
 *   To change this, we'd probably first want to be able to set the
 *   reference gamma values on a per port basis.  Also, some changes
 *   would be needed in the calibration calculations for the "through"
 *   impedance tests which would see an impedance mismatch.
 *
 *   If not set, the default is 50 ohms.
 */
extern int vnacal_calset_set_z0(vnacal_calset_t *vcsp, double complex z0);

/*
 * vnacal_calset_add_vector: add a data vector to the calibration set
 *   @vcsp: pointer to vnacal_calset_t
 *   @row: VNA detector port (zero-based)
 *   @column: VNA driving port (zero-based)
 *   @term: calibration term
 *   @vector: vector of measured complex voltages
 */
extern int vnacal_calset_add_vector(vnacal_calset_t *vcsp, int row,
	int column, int term, const double complex *vector);

/*
 * vnacal_calset_set_reference: set a scalar reference gamma value
 *   @vcsp: pointer to vnacal_calset_t
 *   @reference: reference index 0-3
 *   @gamma: gamma value, e.g. -1 short, 1 open, 0 load
 */
extern int vnacal_calset_set_reference(vnacal_calset_t *vcsp,
	int reference, double complex gamma);

/*
 * vnacal_calset_set_reference_vector: set a vector of reference gamma values
 *   @vcsp: pointer to vnacal_calset_t
 *   @reference: reference index 0-3
 *   @frequencies: length of frequency_vector and gamma_vector
 *   @frequency_vector: vector of increasing frequency values
 *   @gamma_vector: vector of gamma values
 *
 *   The frequency vector given to this function must span the full
 *   range of that given to vnacal_calset_set_frequency_vector, but doesn't
 *   have to be identical.
 */
extern int vnacal_calset_set_reference_vector(vnacal_calset_t *vcsp,
	int reference, int frequencies, const double *frequency_vector,
	const double complex *gamma_vector);

/*
 * vnacal_calset_free: free a vnacal_calset_t structure
 *   @vcsp: pointer to vnacal_calset_t
 *
 *   Free vcs_frequency_vector, vcs_matrix, vcd_sii_reference_vector,
 *   vcd_sii_through_vector, vcd_sji_through_vector and vcd_sdj_leakage.
 */
extern void vnacal_calset_free(vnacal_calset_t *vcsp);

/*
 * vnacal_create: construct a calibration structure from measured data
 *   @sets: number of calibration sets
 *   @vcspp: vector of pointers to vnacal_calset_t structures (sets long)
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
extern vnacal_t *vnacal_create(int sets, vnacal_calset_t **vcspp,
	vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnacal_load: load the calibration from a file
 *   @pathname: calibration file name
 *   @dotdir: directory under $HOME or NULL
 *   @error_fn: error reporting callback or NULL
 *   @error_arg: arbitrary argument passed through to error_fn or NULL
 *
 *   If error_fn is non-NULL, then vnacal_load and subsequent functions report
 *   error messages using error_fn before returning failure to the caller.
 */
extern vnacal_t *vnacal_load(const char *pathname, const char *dotdir,
	vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnacal_save: create or overwrite a calibration file with new data
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @pathname: calibration file name
 *   @dotdir: directory under $HOME or NULL
 *
 *   The pathname and basename parameters work as in vnacal_load except
 *   that the $HOME/{pathname} directory path is created if necessary.
 */
extern int vnacal_save(vnacal_t *vcp, const char *pathname, const char *dotdir);

/*
 * vnacal_get_filename: return the calibration file name
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *
 * Return:
 *   NULL if the vnacal_t structure came from vnacal_create and vnacal_save
 *   hasn't get been called.
 */
extern const char *vnacal_get_filename(const vnacal_t *vcp);

/*
 * vnacal_get_sets: return the number of calibration sets
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 */
extern int vnacal_get_sets(const vnacal_t *vcp);

/*
 * vnacal_get_setname: return the name of the given calibration set
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern const char *vnacal_get_setname(const vnacal_t *vcp, int set);

/*
 * vnacal_get_rows: return the number of rows in the calibration matrix
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern int vnacal_get_rows(const vnacal_t *vcp, int set);

/*
 * vnacal_get_columns: return the number of columns in the calibration matrix
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern int vnacal_get_columns(const vnacal_t *vcp, int set);

/*
 * vnacal_get_frequencies: return the number of frequency points
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern int vnacal_get_frequencies(const vnacal_t *vcp, int set);

/*
 * vnacal_get_fmin: return the minimum frequency point
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern double vnacal_get_fmin(const vnacal_t *vcp, int set);

/*
 * vnacal_get_fmax: return the maximum frequency point
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern double vnacal_get_fmax(const vnacal_t *vcp, int set);

/*
 * vnacal_get_frequency_vector: return a pointer to the calibration frequencies
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 */
extern const double *vnacal_get_frequency_vector(const vnacal_t *vcp, int set);

/*
 * vnacal_set_fprecision: set the frequency value precision for vnacal_save
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @precision: precision in decimal places (1..n) or VNACAL_MAX_PRECISION
 */
extern int vnacal_set_fprecision(vnacal_t *vcp, int precision);

/*
 * vnacal_set_dprecision: set the data value precision for vnacal_save
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @precision: precision in decimal places (1..n) or VNACAL_MAX_PRECISION
 */
extern int vnacal_set_dprecision(vnacal_t *vcp, int precision);

/*
 * vnacal_property_type: get the type of the given property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_type(vnacal_t *vcp, int set,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/*
 * vnacal_property_count: return count of structures in given collection
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_count(vnacal_t *vcp, int set,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/*
 * vnacal_property_keys: return a vector of keys for the given map expr
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 *
 * Caller can free the vector by a call to free.
 */
extern const char **vnacal_property_keys(vnacal_t *vcp, int set,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/*
 * vnacal_property_get: get a property value from a property expression
 *   @root:   property data root (can be NULL)
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern const char *vnacal_property_get(vnacal_t *vcp, int set,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/*
 * vnacal_property_set: set a property value from a property expression
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_set(vnacal_t *vcp, int set,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/*
 * vnacal_property_delete: delete the value described by format
 *   @anchor: address of root property pointer
 *   @format: printf-like format string forming the property expression
 *   @...:    optional variable arguments
 */
extern int vnacal_property_delete(vnacal_t *vcp, int set,
	const char *format, ...)
#ifdef __GNUC__
    __attribute__((__format__(__printf__, 3, 4)));
#else
    ;
#endif

/*
 * vnacal_free: free a vnacal_t structure
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 */
extern void vnacal_free(vnacal_t *vcp);

/*
 * vnacal_apply: apply the calibration to measured values (simple form)
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 *   @frequencies: number of frequencies in matrix and s_parameters
 *   @frequency_vector: vector of increasing frequency points
 *   @matrix: vrows x vcolumns matrix of measured voltage ratios
 *   @s_parameters: caller-allocated vnadata_t structure
 */
extern int vnacal_apply(vnacal_t *vcp, int set, int frequencies,
        const double *frequency_vector, const double complex *const *matrix,
        vnadata_t *s_parameters);

/*
 * vnacal_apply_alloc: allocate a vnacal_apply_t structure
 *   @vcp: a pointer to the structure returned by vnacal_create or vnacal_load
 *   @set: set index (beginning with zero)
 *   @drows: number of rows in the DUT S-parameter matrix
 *   @dcolumns: number of columns in the DUT S-parameter matrix
 *   @dfrequencies: number of measured frequencies
 */
extern vnacal_apply_t *vnacal_apply_alloc(vnacal_t *vcp, int set,
	int drows, int dcolumns, int dfrequencies);

/*
 * vnacal_apply_set_frequency_vector: set the frequency vector
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @frequency_vector: vector of dfrequencies measured frequencies
 *
 * The frequencies don't have to be the same as those in the calibration;
 * however, they may not extend outside of the calibration range -- the
 * the library will interpolate, but won't extrapolate.
 */
extern int vnacal_apply_set_frequency_vector(vnacal_apply_t *vap,
	const double *frequency_vector);

/*
 * vnacal_apply_add_column: add a column of measurements with map
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @vcolumn: column of the calibration matrix being added
 *   @column_vector: pointer to vector, vrows long, of pointers
 *      to vectors of voltage measurements, dfrequencies long
 *   @map: map of VNA port to DUT port, MAX(vrows, vcolumns) long
 *
 * When the DUT has more ports than the VNA, this function gives the
 * mapping for this set of measurements.  The map must be the greater
 * dimension of the calibration matrix.  DUT port indices start at zero.
 * The special value -1 indicates that the VNA port and any unused DUT
 * ports are connected to terminators.
 */
extern int vnacal_apply_add_column(vnacal_apply_t *vap,
	int vcolumn, const double complex *const *column_vector,
	const int *map);

/*
 * vnacal_apply_add_matrix: add matrix of vector with port map
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @matrix: vrows x vcolumns matrix of dfrequencies long vectors of voltages
 *   @map: map of VNA port to DUT port, MAX(vrows, vcolumns) long
 */
extern int vnacal_apply_add_matrix(vnacal_apply_t *vap,
	const double complex *const *matrix, const int *map);

/*
 * vnacal_apply_get_data: create the s-parameter matrix
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 *   @s_parameters: caller-allocated vnadata_t structure
 */
extern int vnacal_apply_get_data(const vnacal_apply_t *vap,
	vnadata_t *s_parameters);

/*
 * vnacal_apply_free: free the apply structure
 *   @vap: pointer to the structure returned from vnacal_apply_alloc
 */
extern void vnacal_apply_free(vnacal_apply_t *vap);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_H */
