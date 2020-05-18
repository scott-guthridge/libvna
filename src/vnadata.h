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

#ifndef VNADATA_H
#define VNADATA_H

#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * VNADATA_DEFAULT_Z0: default system impedance
 */
#define VNADATA_DEFAULT_Z0	50.0

/*
 * vnadata_parameter_type_t: network parameter data type
 *   When updating, also update vnadata_get_typename.
 */
typedef enum vnadata_parameter_type {
    VPT_UNDEF	=  0,
    VPT_S	=  1,
    VPT_Z	=  2,
    VPT_Y	=  3,
    VPT_T	=  4,
    VPT_H	=  5,
    VPT_G	=  6,
    VPT_A	=  7,
    VPT_B	=  8,
    VPT_ZIN	=  9,
    VPT_NTYPES	= 10,
} vnadata_parameter_type_t;

/*
 * vnadata_t: network parameter data
 *
 * Note: The members of this structure should be treated as opaque
 * by users of the library.  Accessing these directly will expose
 * you to future source-level compatibility breaks.
 */
typedef struct vnadata {
    vnadata_parameter_type_t vd_type;
    int vd_frequencies;
    int vd_rows;
    int vd_columns;
    double *vd_frequency_vector;
    double complex **vd_data;
} vnadata_t;

/*
 * vnadata_alloc: allocate an empty vnadata_t object
 */
extern vnadata_t *vnadata_alloc();

/*
 * vnadata_free: free a vnadata_t object
 *   @vdp: object to free
 */
extern void vnadata_free(vnadata_t *vdp);

/*
 * vnadata_init: resize and initialize a vnadata_t object
 *   @frequencies: number of frequency points
 *   @rows: number of matrix rows
 *   @columns: number of matrix columns
 *   @type: matrix type (see above)
 */
extern int vnadata_init(vnadata_t *vdp, int frequencies, int rows,
	int columns, vnadata_parameter_type_t type);

/*
 * vnadata_alloc_and_init: allocate a vnadata_t object and initialize to zero
 *   @frequencies: number of frequency points
 *   @rows: number of matrix rows
 *   @columns: number of matrix columns
 *   @type: matrix type (see above)
 */
static inline vnadata_t *vnadata_alloc_and_init(int frequencies, int rows,
	int columns, vnadata_parameter_type_t type)
{
    vnadata_t *vdp;

    if ((vdp = vnadata_alloc()) == NULL) {
	return NULL;
    }
    if (vnadata_init(vdp, frequencies, rows, columns, type) == -1) {
	vnadata_free(vdp);
	return NULL;
    }
    return vdp;
}

/*
 * vnadata_resize: resize a vnadata_t object without clearing values
 *   @vdp: object pointer
 *   @frequencies: new number of frequencies
 *   @rows: new number of rows
 *   @columns: new number of columns
 *   @type: new parameter type
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
extern int vnadata_resize(vnadata_t *vdp, int frequencies,
	int rows, int columns, vnadata_parameter_type_t type);

/*
 * vnadata_get_frequencies: return the number of frequencies
 *   @vdp: vnadata object pointer
 */
static inline int vnadata_get_frequencies(const vnadata_t *vdp)
{
    return vdp->vd_frequencies;
}

/*
 * vnadata_get_rows: return the number of rows
 *   @vdp: vnadata object pointer
 */
static inline int vnadata_get_rows(const vnadata_t *vdp)
{
    return vdp->vd_rows;
}

/*
 * vnadata_get_columns: return the number of columns
 *   @vdp: vnadata object pointer
 */
static inline int vnadata_get_columns(const vnadata_t *vdp)
{
    return vdp->vd_columns;
}

/*
 * vnadata_get_type: return the type of data in the vnadata_t object
 *   @vdp: vnadata object pointer
 */
static inline vnadata_parameter_type_t vnadata_get_type(const vnadata_t *vdp)
{
    return vdp->vd_type;
}

/*
 * vnadata_set_type: change the parameter type without conversion
 *   @vdp: vnadata object pointer
 *   @type: new parameter type
 */
extern int vnadata_set_type(vnadata_t *vdp, vnadata_parameter_type_t type);

/*
 * vnadata_get_frequency: get the indexed frequency
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 */
static inline double vnadata_get_frequency(const vnadata_t *vdp, int findex)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	errno = EINVAL;
	return HUGE_VAL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_frequency_vector[findex];
}

/*
 * vnadata_set_frequency: get the indexed frequency
 *   @vdp:       vnadata object pointer
 *   @findex:    frequency index
 *   @frequency: new frequency value
 */
static inline int vnadata_set_frequency(vnadata_t *vdp, int findex,
	double frequency)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	errno = EINVAL;
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    vdp->vd_frequency_vector[findex] = frequency;
    return 0;
}

/*
 * vnadata_get_frequency_vector: get the frequency vector
 *   @vdp: vnadata object pointer
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
 *   @vdp: vnadata object pointer
 *   @findex: frequency index
 *   @frequency_vector: new frequency value
 */
static inline int vnadata_set_frequency_vector(vnadata_t *vdp,
	const double *frequency_vector)
{
    (void)memcpy((void *)vdp->vd_frequency_vector, (void *)frequency_vector,
	vdp->vd_frequencies * sizeof(double));
    return 0;
}

/*
 * vnadata_get_cell: get a value from the matrix
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 *   @row:    matrix row
 *   @column: matrix column
 */
static inline double complex vnadata_get_cell(const vnadata_t *vdp,
	int findex, int row, int column)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies ||
	row    < 0 || row    >= vdp->vd_rows ||
	column < 0 || column >= vdp->vd_columns) {
	errno = EINVAL;
	return HUGE_VAL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_data[findex][row * vdp->vd_columns + column];
}

/*
 * vnadata_set_cell: set a matrix value
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 *   @row:    matrix row
 *   @column: matrix column
 *   @value:  value to write
 */
static inline int vnadata_set_cell(vnadata_t *vdp, int findex, int row,
	int column, double complex value)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies ||
	row    < 0 || row    >= vdp->vd_rows ||
	column < 0 || column >= vdp->vd_columns) {
	errno = EINVAL;
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    vdp->vd_data[findex][row * vdp->vd_columns + column] = value;
    return 0;
}

/*
 * vnadata_get_matrix: return the serialized matrix at the given freq. index
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 */
static inline const double complex *vnadata_get_matrix(const vnadata_t *vdp,
	int findex)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	errno = EINVAL;
	return NULL;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    return vdp->vd_data[findex];
}

/*
 * vnadata_set_matrix: set the matrix at the given frequency index
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 *   @matrix: serialized matrix (concatenation of rows)
 */
static inline int vnadata_set_matrix(vnadata_t *vdp, int findex,
	const double complex *matrix)
{
#ifndef VNADATA_NO_BOUNDS_CHECK
    if (findex < 0 || findex >= vdp->vd_frequencies) {
	errno = EINVAL;
	return -1;
    }
#endif /* VNADATA_NO_BOUNDS_CHECK */
    (void)memcpy((void *)vdp->vd_data[findex], (void *)matrix,
	vdp->vd_rows * vdp->vd_columns * sizeof(double complex));
    return 0;
}

/*
 * vnadata_get_from_vector: copy a matrix cell into a by-frequency vector
 *   @vdp:    vnadata object pointer
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
    if (row    < 0 || row    >= vdp->vd_rows ||
	column < 0 || column >= vdp->vd_columns) {
	errno = EINVAL;
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
 *   @vdp:    vnadata object pointer
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
    if (row    < 0 || row    >= vdp->vd_rows ||
	column < 0 || column >= vdp->vd_columns) {
	errno = EINVAL;
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
 *   @vdp:  vnadata object pointer
 *   @port: port number (zero-based)
 *
 * Note: fails if vdp contains per-frequency z0 vectors.
 */
extern double complex vnadata_get_z0(const vnadata_t *vdp, int port);

/*
 * vnadata_set_z0: set the z0 value for the given port
 *   @vdp:  vnadata object pointer
 *   @port: port number (zero-based)
 *   @z0:   new value
 */
extern int vnadata_set_z0(vnadata_t *vdp, int port, double complex z0);

/*
 * vnadata_set_all_z0: set all z0's to the same value
 *   @vdp:  vnadata object pointer
 *   @z0:   new value
 */
extern int vnadata_set_all_z0(vnadata_t *vdp, double complex z0);

/*
 * vnadata_get_z0_vector: return the z0 vector
 *   @vdp:  vnadata object pointer
 *
 * Note: fails if vdp contains per-frequency z0 vectors.
 */
extern const double complex *vnadata_get_z0_vector(const vnadata_t *vdp);

/*
 * vnadata_set_z0_vector: set the z0 vector
 *   @vdp:       vnadata object pointer
 *   @z0_vector: new values (length is max of rows and columns)
 */
extern int vnadata_set_z0_vector(vnadata_t *vdp,
	const double complex *z0_vector);

/*
 * vnadata_get_fz0: return the z0 value for the given frequency and port
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 *   @port:   port number (zero-based)
 */
extern double complex vnadata_get_fz0(const vnadata_t *vdp, int findex,
	int port);

/*
 * vnadata_set_fz0: set the z0 value for the given frequency and port
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 *   @port:   port number (zero-based)
 *   @z0:     new value
 */
extern int vnadata_set_fz0(vnadata_t *vdp, int findex, int port,
	double complex z0);

/*
 * vnadata_get_fz0_vector: return the z0 vector for the given frequency
 *   @vdp:    vnadata object pointer
 *   @findex: frequency index
 */
extern const double complex *vnadata_get_fz0_vector(const vnadata_t *vdp,
	int findex);

/*
 * vnadata_set_fz0_vector: set the z0 vector for the given frequency
 *   @vdp:       vnadata object pointer
 *   @findex:    frequency index
 *   @z0_vector: new values (length is max of rows and columns)
 */
extern int vnadata_set_fz0_vector(vnadata_t *vdp, int findex,
	const double complex *z0_vector);

/*
 * vnadata_convert: convert the matrix from one parameter type to another
 *   @vdp_in:  input vnadata object
 *   @vdp_out: output vnadata object (may be same as vdp_in)
 *   @newtype: new type
 */
extern int vnadata_convert(const vnadata_t *vdp_in,
	vnadata_t *vdp_out, vnadata_parameter_type_t newtype);

/*
 * vnadata_add_frequency: add a new frequency entry
 *   @vdp: object pointer
 *   @frequency: new frequency value
 *
 *   Increase the number of frequencies in the data set by one creating
 *   zero-filled data elements.  This is useful when parsing Touchstone
 *   version 1 files where we don't know the number of frequencies
 *   up front.
 */
extern int vnadata_add_frequency(vnadata_t *vdp, double frequency);

/*
 * vnadata_get_typename: convert parameter type to name
 *   @type: parameter type
 */
extern const char *vnadata_get_typename(vnadata_parameter_type_t type);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* VNADATA_H */
