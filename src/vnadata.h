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

#ifndef VNADATA_H
#define VNADATA_H

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <vnaerr.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * VNADATA_DEFAULT_Z0: default reference impedance
 */
#define VNADATA_DEFAULT_Z0	50.0

/*
 * vnadata_parameter_type_t: parameter type
 *   When changing, also update vnadata_get_type_name and conversion_table
 *   in vnadata_convert.c.
 */
typedef enum vnadata_parameter_type {
    VPT_UNDEF	=  0,
    VPT_S	=  1,
    VPT_T	=  2,
    VPT_U	=  3,
    VPT_Z	=  4,
    VPT_Y	=  5,
    VPT_H	=  6,
    VPT_G	=  7,
    VPT_A	=  8,
    VPT_B	=  9,
    VPT_ZIN	= 10,
    VPT_NTYPES	= 11,
} vnadata_parameter_type_t;

/*
 * VNADATA_MAX_PRECISION: argument to vnadata_set_fprecision and
 *      vnadata_set_dprecision for hexadecimal floating point
 *
 * Note: must be the same as VNACAL_MAX_PRECISION
 */
#define VNADATA_MAX_PRECISION	1000

/* vnadata_filetype_t: file type */
typedef enum vnadata_filetype {
	/* automatically determine format from the filename */
	VNADATA_FILETYPE_AUTO		= 0,
	/* touchstone v1 format */
	VNADATA_FILETYPE_TOUCHSTONE1	= 1,
	/* touchstone v2 format */
	VNADATA_FILETYPE_TOUCHSTONE2	= 2,
	/* network parameter data format */
	VNADATA_FILETYPE_NPD		= 3,
} vnadata_filetype_t;

/*
 * vnadata_t: network parameter data
 *
 * Note: The members of this structure should be treated as opaque
 * by users of the library.  Accessing these directly will expose
 * you to future compatibility breaks.
 */
typedef struct vnadata {
    vnadata_parameter_type_t vd_type;
    int vd_rows;
    int vd_columns;
    int vd_frequencies;
    double *vd_frequency_vector;
    double complex **vd_data;
} vnadata_t;

/*
 * vnadata_alloc: allocate an empty vnadata_t structure
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 */
extern vnadata_t *vnadata_alloc(vnaerr_error_fn_t *error_fn, void *error_arg);

/*
 * vnadata_free: free a vnadata_t structure
 *   @vdp: pointer to vnacal_data_t structure to free
 */
extern void vnadata_free(vnadata_t *vdp);

/*
 * vnadata_init: resize and initialize a vnadata_t structure
 *   @type: parameter type (see above)
 *   @rows: number of matrix rows
 *   @columns: number of matrix columns
 *   @frequencies: number of frequency points
 */
extern int vnadata_init(vnadata_t *vdp, vnadata_parameter_type_t type,
	int rows, int columns, int frequencies);

/*
 * _vnadata_bounds_error: internal function to report bounds error
 *   @function: name of calling function
 *   @vdp: pointer to vnacal_data_t structure
 *   @what: name of parameter
 *   @value: value of parameter
 */
extern void _vnadata_bounds_error(const char *function, const vnadata_t *vdp,
	const char *what, int value);

/*
 * vnadata_alloc_and_init: allocate a vnadata_t structure and initialize to zero
 *   @error_fn: optional error reporting function (NULL if not used)
 *   @error_arg: user data passed through to the error function (or NULL)
 *   @type: parameter type (see above)
 *   @rows: number of matrix rows
 *   @columns: number of matrix columns
 *   @frequencies: number of frequency points
 */
static inline vnadata_t *vnadata_alloc_and_init(vnaerr_error_fn_t *error_fn,
	void *error_arg, vnadata_parameter_type_t type,
	int rows, int columns, int frequencies)
{
    vnadata_t *vdp;

    if ((vdp = vnadata_alloc(error_fn, error_arg)) == NULL) {
	return NULL;
    }
    if (vnadata_init(vdp, type, rows, columns, frequencies) == -1) {
	vnadata_free(vdp);
	return NULL;
    }
    return vdp;
}

/*
 * vnadata_resize: resize a vnadata_t structure without clearing values
 *   @vdp: pointer to vnacal_data_t structure
 *   @type: new parameter type
 *   @rows: new number of rows
 *   @columns: new number of columns
 *   @frequencies: new number of frequencies
 *
 * Notes:
 *   When changing dimensions, this function increases or decreases
 *   the allocation as necessary and initializes any newly allocated
 *   or vacated cells, but doesn't otherwise reorganize existing data.
 *   Because of the way the the data are stored, changing the number
 *   of frequencies or the number of rows is value preserving, as is
 *   changing between a row vector and a column vector; however, changing
 *   the number of columns, in general, results in values appearing
 *   in the wrong cells.  This can actually be desirable, for example,
 *   when transforming the first row of a 4x4 matrix incorrectly stored
 *   as a 2x2 matrix.
 *
 *   Changing the dimensions invalidates earlier pointers returned
 *   from vnadata_get_frequency_vector, vnadata_get_matrix and
 *   vnadata_get_z0_vector.
 */
extern int vnadata_resize(vnadata_t *vdp, vnadata_parameter_type_t type,
	int rows, int columns, int frequencies);

/*
 * vnadata_get_frequencies: return the number of frequencies
 *   @vdp: a pointer to the vnadata_t structure
 */
static inline int vnadata_get_frequencies(const vnadata_t *vdp)
{
    return vdp->vd_frequencies;
}

/*
 * vnadata_get_rows: return the number of rows
 *   @vdp: a pointer to the vnadata_t structure
 */
static inline int vnadata_get_rows(const vnadata_t *vdp)
{
    return vdp->vd_rows;
}

/*
 * vnadata_get_columns: return the number of columns
 *   @vdp: a pointer to the vnadata_t structure
 */
static inline int vnadata_get_columns(const vnadata_t *vdp)
{
    return vdp->vd_columns;
}

/*
 * vnadata_get_type: return the parameter type of the vnadata_t structure
 *   @vdp: a pointer to the vnadata_t structure
 */
static inline vnadata_parameter_type_t vnadata_get_type(const vnadata_t *vdp)
{
    return vdp->vd_type;
}

/*
 * vnadata_set_type: change the parameter type without conversion
 *   @vdp: a pointer to the vnadata_t structure
 *   @type: new parameter type
 */
extern int vnadata_set_type(vnadata_t *vdp, vnadata_parameter_type_t type);

/*
 * vnadata_get_fmin: get the minimum frequency
 *   @vdp:    a pointer to the vnadata_t structure
 */
static inline double vnadata_get_fmin(const vnadata_t *vdp)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    if (vdp->vd_frequencies == 0) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", 0);
	return HUGE_VAL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_frequency_vector[0];
}

/*
 * vnadata_get_fmax: get the maximum frequency
 *   @vdp:    a pointer to the vnadata_t structure
 */
static inline double vnadata_get_fmax(const vnadata_t *vdp)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    if (vdp->vd_frequencies == 0) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", 0);
	return HUGE_VAL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_frequency_vector[vdp->vd_frequencies - 1];
}

/*
 * vnadata_get_frequency: get the indexed frequency
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 */
static inline double vnadata_get_frequency(const vnadata_t *vdp, int findex)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", findex);
	return HUGE_VAL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_frequency_vector[findex];
}

/*
 * vnadata_set_frequency: get the indexed frequency
 *   @vdp:       a pointer to the vnadata_t structure
 *   @findex:    frequency index
 *   @frequency: new frequency value
 */
static inline int vnadata_set_frequency(vnadata_t *vdp, int findex,
	double frequency)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", findex);
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    vdp->vd_frequency_vector[findex] = frequency;
    return 0;
}

/*
 * vnadata_get_frequency_vector: get the frequency vector
 *   @vdp: a pointer to the vnadata_t structure
 */
static inline const double *vnadata_get_frequency_vector(const vnadata_t *vdp)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return NULL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_frequency_vector;
}

/*
 * vnadata_set_frequency_vector: set the frequency vector
 *   @vdp: a pointer to the vnadata_t structure
 *   @findex: frequency index
 *   @frequency_vector: new frequency value
 */
static inline int vnadata_set_frequency_vector(vnadata_t *vdp,
	const double *frequency_vector)
{
    if (vdp == NULL || frequency_vector == NULL) {
	errno = EINVAL;
	return -1;
    }
    (void)memcpy((void *)vdp->vd_frequency_vector, (void *)frequency_vector,
	vdp->vd_frequencies * sizeof(double));
    return 0;
}

/*
 * vnadata_get_cell: get a value from the matrix
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 *   @row:    matrix row
 *   @column: matrix column
 */
static inline double complex vnadata_get_cell(const vnadata_t *vdp,
	int findex, int row, int column)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return HUGE_VAL;
    }
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", findex);
	return HUGE_VAL;
    }
    if (row    < 0 || row    >= vdp->vd_rows) {
	_vnadata_bounds_error(__func__, vdp, "row", row);
	return HUGE_VAL;
    }
    if (column < 0 || column >= vdp->vd_columns) {
	_vnadata_bounds_error(__func__, vdp, "column", column);
	return HUGE_VAL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_data[findex][row * vdp->vd_columns + column];
}

/*
 * vnadata_set_cell: set a matrix value
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 *   @row:    matrix row
 *   @column: matrix column
 *   @value:  value to write
 */
static inline int vnadata_set_cell(vnadata_t *vdp, int findex, int row,
	int column, double complex value)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", findex);
	return -1;
    }
    if (row    < 0 || row    >= vdp->vd_rows) {
	_vnadata_bounds_error(__func__, vdp, "row", row);
	return -1;
    }
    if (column < 0 || column >= vdp->vd_columns) {
	_vnadata_bounds_error(__func__, vdp, "column", column);
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    vdp->vd_data[findex][row * vdp->vd_columns + column] = value;
    return 0;
}

/*
 * vnadata_get_matrix: return the serialized matrix at the given freq. index
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 */
static inline double complex *vnadata_get_matrix(const vnadata_t *vdp,
	int findex)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return NULL;
    }
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", findex);
	return NULL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_data[findex];
}

/*
 * vnadata_set_matrix: set the matrix at the given frequency index
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 *   @matrix: serialized matrix (concatenation of rows)
 */
static inline int vnadata_set_matrix(vnadata_t *vdp, int findex,
	const double complex *matrix)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	_vnadata_bounds_error(__func__, vdp, "frequency index", findex);
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    (void)memcpy((void *)vdp->vd_data[findex], (void *)matrix,
	vdp->vd_rows * vdp->vd_columns * sizeof(double complex));
    return 0;
}

/*
 * vnadata_get_from_vector: copy a matrix cell into a by-frequency vector
 *   @vdp:    a pointer to the vnadata_t structure
 *   @row:    matrix row
 *   @column: matrix column
 *   @vector: vector of data values by frequency
 *
 * Vector must by frequencies entries long.
 */
static inline int vnadata_get_to_vector(const vnadata_t *vdp,
	int row, int column, double complex *vector)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (row    < 0 || row    >= vdp->vd_rows) {
	_vnadata_bounds_error(__func__, vdp, "row", row);
	return -1;
    }
    if (column < 0 || column >= vdp->vd_columns) {
	_vnadata_bounds_error(__func__, vdp, "column", column);
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    for (int findex = 0; findex < vdp->vd_frequencies; ++findex) {
	vector[findex] = vdp->vd_data[findex][row * vdp->vd_columns + column];
    }
    return 0;
}

/*
 * vnadata_set_from_vector: set a matrix cell from by-frequency vector
 *   @vdp:    a pointer to the vnadata_t structure
 *   @row:    matrix row
 *   @column: matrix column
 *   @vector: vector of data values by frequency
 *
 * Vector must by frequencies entries long.
 */
static inline int vnadata_set_from_vector(vnadata_t *vdp, int row, int column,
	const double complex *vector)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (vdp == NULL) {
	errno = EINVAL;
	return -1;
    }
    if (row    < 0 || row    >= vdp->vd_rows) {
	_vnadata_bounds_error(__func__, vdp, "row", row);
	return -1;
    }
    if (column < 0 || column >= vdp->vd_columns) {
	_vnadata_bounds_error(__func__, vdp, "column", column);
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    for (int findex = 0; findex < vdp->vd_frequencies; ++findex) {
	vdp->vd_data[findex][row * vdp->vd_columns + column] = vector[findex];
    }
    return 0;
}

/*
 * vnadata_get_z0: return the z0 value for the given port
 *   @vdp:  a pointer to the vnadata_t structure
 *   @port: port number (zero-based)
 *
 * Note: fails if vdp contains per-frequency z0 vectors.
 */
extern double complex vnadata_get_z0(const vnadata_t *vdp, int port);

/*
 * vnadata_set_z0: set the z0 value for the given port
 *   @vdp:  a pointer to the vnadata_t structure
 *   @port: port number (zero-based)
 *   @z0:   new value
 */
extern int vnadata_set_z0(vnadata_t *vdp, int port, double complex z0);

/*
 * vnadata_set_all_z0: set all z0's to the same value
 *   @vdp:  a pointer to the vnadata_t structure
 *   @z0:   new value
 */
extern int vnadata_set_all_z0(vnadata_t *vdp, double complex z0);

/*
 * vnadata_get_z0_vector: return the z0 vector
 *   @vdp:  a pointer to the vnadata_t structure
 *
 * Note: fails if vdp contains per-frequency z0 vectors.
 */
extern const double complex *vnadata_get_z0_vector(const vnadata_t *vdp);

/*
 * vnadata_set_z0_vector: set the z0 vector
 *   @vdp:       a pointer to the vnadata_t structure
 *   @z0_vector: new values (length is max of rows and columns)
 */
extern int vnadata_set_z0_vector(vnadata_t *vdp,
	const double complex *z0_vector);

/*
 * vnadata_has_fz0: return true if reference impedances are per-frequency
 *   @vdp:       a pointer to the vnadata_t structure
 */
extern bool vnadata_has_fz0(const vnadata_t *vdp);

/*
 * vnadata_get_fz0: return the z0 value for the given frequency and port
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 *   @port:   port number (zero-based)
 */
extern double complex vnadata_get_fz0(const vnadata_t *vdp, int findex,
	int port);

/*
 * vnadata_set_fz0: set the z0 value for the given frequency and port
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 *   @port:   port number (zero-based)
 *   @z0:     new value
 */
extern int vnadata_set_fz0(vnadata_t *vdp, int findex, int port,
	double complex z0);

/*
 * vnadata_get_fz0_vector: return the z0 vector for the given frequency
 *   @vdp:    a pointer to the vnadata_t structure
 *   @findex: frequency index
 */
extern const double complex *vnadata_get_fz0_vector(const vnadata_t *vdp,
	int findex);

/*
 * vnadata_set_fz0_vector: set the z0 vector for the given frequency
 *   @vdp:       a pointer to the vnadata_t structure
 *   @findex:    frequency index
 *   @z0_vector: new values (length is max of rows and columns)
 */
extern int vnadata_set_fz0_vector(vnadata_t *vdp, int findex,
	const double complex *z0_vector);

/*
 * vnadata_convert: convert the matrix from one parameter type to another
 *   @vdp_in:  input network parameter data
 *   @vdp_out: output network parameter data (can be same as vdp_in)
 *   @newtype: new parameter type (can be the same as old)
 */
extern int vnadata_convert(const vnadata_t *vdp_in, vnadata_t *vdp_out,
	vnadata_parameter_type_t new_type);

/*
 * vnadata_rconvert: renormalizing conversion
 *   @vdp_in:  input network parameter data
 *   @vdp_out: output network parameter data (can be same as vdp_in)
 *   @newtype: new parameter type (can be the same as old)
 *   @new_z0:  new reference impedances
 *   @new_z0_length:  length of new_z0 (1, #ports, or #frequencies * #ports)
 *
 * Note: vdp_out and vdp_in may be the same.  The length field
 * gives the number of elements in new_z0.  It can be 1, #ports,
 * or #frequencies x #ports.
 */
extern int vnadata_rconvert(const vnadata_t *vdp_in, vnadata_t *vdp_out,
	vnadata_parameter_type_t newtype, const double complex *new_z0,
	int new_z0_length);

/*
 * vnadata_add_frequency: add a new frequency entry
 *   @vdp: a pointer to the vnadata_t structure
 *   @frequency: new frequency value
 *
 *   Increase the number of frequencies in the data set by one creating
 *   zero-filled data elements.  This is useful when parsing Touchstone
 *   version 1 files where we don't know the number of frequencies
 *   up front.
 */
extern int vnadata_add_frequency(vnadata_t *vdp, double frequency);

/*
 * vnadata_get_type_name: convert parameter type to name
 *   @type: parameter type
 */
extern const char *vnadata_get_type_name(vnadata_parameter_type_t type);

/*
 * vnadata_get_filetype: return the current file type
 *   @vdp: a pointer to the vnadata_t structure
 */
extern vnadata_filetype_t vnadata_get_filetype(const vnadata_t *vdp);

/*
 * vnadata_set_filetype: set the file type
 *   @vdp: a pointer to the vnadata_t structure
 *   @filetype: file type
 *
 *   The default type is VNADATA_FILETYPE_AUTO where the library tries to
 *   intuit the type from the filename.
 */
extern int vnadata_set_filetype(vnadata_t *vdp, vnadata_filetype_t filetype);

/*
 * vnadata_get_format: get the load/save format string
 *   @vdp: pointer to the structure returned from vnadata_alloc
 */
extern const char *vnadata_get_format(const vnadata_t *vdp);

/*
 * vnadata_set_format: set the load/save format string
 *   @vdp: a pointer to the vnadata_t structure
 *   @format: a comma-separated case-insensitive list of the following:
 *     {S,Z,Y,T,H,G,A,B}[{ri,ma,dB}]
 *     {il,rl}
 *     zin[{ri,ma}]
 *     {prc,prl,src,srl}
 *     vswr
 *
 *   If not set, a suitable default will be provided.
 */
extern int vnadata_set_format(vnadata_t *vdp, const char *format);

/*
 * vnadata_get_fprecision: get the frequency value precision
 *   @vdp: a pointer to the vnadata_t structure
 */
extern int vnadata_get_fprecision(const vnadata_t *vdp);

/*
 * vnadata_set_fprecision: set the frequency value precision
 *   @vdp: a pointer to the vnadata_t structure
 *   @precision: precision in decimal places (1..n) or VNADATA_MAX_PRECISION
 */
extern int vnadata_set_fprecision(vnadata_t *vdp, int precision);

/*
 * vnadata_get_dprecision: set the data value precision
 *   @vdp: a pointer to the vnadata_t structure
 */
extern int vnadata_get_dprecision(const vnadata_t *vdp);

/*
 * vnadata_set_dprecision: set the data value precision for vnadata_save
 *   @vdp: a pointer to the vnadata_t structure
 *   @precision: precision in decimal places (1..n) or VNADATA_MAX_PRECISION
 */
extern int vnadata_set_dprecision(vnadata_t *vdp, int precision);

/*
 * vnadata_load: load network parameters from filename
 *   @vdp: a pointer to the vnadata_t structure (reshaped as needed)
 *   @filename: file to load
 */
extern int vnadata_load(vnadata_t *vdp, const char *filename);

/*
 * vnadata_load: load network parameters from a file pointer
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 */
extern int vnadata_fload(vnadata_t *vdp, FILE *fp, const char *filename);

/*
 * vnadata_cksave: check that the given parameters and format are valid for save
 *   @vdp: a pointer to the vnadata_t structure
 *   @filename: file to save
 */
extern int vnadata_cksave(vnadata_t *vdp, const char *filename);

/*
 * vnadata_save: save network parameters to filename
 *   @vdp: a pointer to the vnadata_t structure
 *   @filename: file to save
 */
extern int vnadata_save(vnadata_t *vdp, const char *filename);

/*
 * vnadata_fsave: save network parameters to a file pointer
 *   @vdp: a pointer to the vnadata_t structure
 *   @fp: file pointer
 *   @filename: filename used in error messages and to intuit the file type
 */
extern int vnadata_fsave(vnadata_t *vdp, FILE *fp, const char *filename);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNADATA_H */
