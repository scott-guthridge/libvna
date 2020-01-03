/*
 * Vector Network Analyzer Calibration Library
 * Copyright © 2020 D Scott Guthridge <scott_guthridge@rompromity.net>
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

#include <stdio.h>
#include "vnaproperty.h"
#include "vnacal.h"

#ifdef __cplusplus
extern "C" {
#endif

/* default numerical precision for vnacal_save frequencies (in digits) */
#define VNACAL_DEFAULT_FREQUENCY_PRECISION	7

/* default numerical precision for vnacal_save data (in digits) */
#define VNACAL_DEFAULT_DATA_PRECISION		6

/* maximum points to use for rational function interpolation */
#define VNACAL_MAX_M				5

/* factor by which we will extrapolate the frequency */
#define VNACAL_F_EXTRAPOLATION			0.01

/*
 * vnacal_calset_reference: a reference value or reference vector
 */
typedef struct vnacal_calset_reference {
    bool vcdsr_is_vector;
    union {
	double complex vcdsr_gamma;
	struct {
	    int vcdsr_frequencies;
	    double *vcdsr_frequency_vector;
	    double complex *vcdsr_gamma_vector;
	} v;
    } u ;
} vnacal_calset_reference_t;

/*
 * vnacal_cdata_t: a cell of the matrix of measured calibration data
 *
 *   Each column of the matrix represents the driving VNA port for a
 *   set of measurements.  Each row of the matrix represents the receiving
 *   VNA port.  The diagonal and off diagonal elements of the matrix
 *   contain different information and use different union members.
 *   Diagonal:
 *     vcd_data_vectors[0]	reflected voltage from reference 0
 *     vcd_data_vectors[1]	reflected voltage from reference 1
 *     vcd_data_vectors[2]	reflected voltage from reference 2
 *   Off-diagnonal:
 *     vcd_data_vectors[0]	reflected signal back to driving port
 *     vcd_data_vectors[1]	transmitted signal to receiving port
 *     vcd_data_vectors[2]	leakage from driving port to isolated
 *     					receiving port
 *
 *   Each member is a pointer to a vector of complex voltage ratios,
 *   one entry for each frequency.
 *
 */
typedef struct vnacal_cdata {

    /* three column vectors of data (see aliases below) */
    double complex *vcd_data_vectors[3];

    /* how many vectors have been added into each slot */
    int vcd_counts[3];

} vnacal_cdata_t;

/*
 * vnacal_calset_t: measured calibration data for vnacal_create
 *   This structure represents one set of VNA calibration data.  A VNA
 *   calibration file may contain more than one set of calibration
 *   data, for example, if a relay-controlled attenuator is used and
 *   it's desired to make a separate calibration for each setting.
 *   Another example is where the reflection bridge can be switched
 *   between returning forward / reflected signal and forward+reflected /
 *   forward-reflected signal.
 */
struct vnacal_calset {

    /* name of this set */
    const char *vcs_setname;

    /* number of ports where signal is measured (must be >= columns) */
    int vcs_rows;

    /* number of ports where signal is generated */
    int vcs_columns;

    /* number of frequencies */
    int vcs_frequencies;

    /* reference gamma values */
    vnacal_calset_reference_t vcs_references[3];

    /* vector of frequencies */
    double *vcs_frequency_vector;

    /* true if the frequency vector has been set */
    bool vcs_frequencies_valid;

    /* system impedance of the VNA ports (currently assumed all
       the same) */
    double complex vcs_z0;

    /* serialized matrix of vnacal_cdata_t objects */
    vnacal_cdata_t *vcs_matrix;

    /* user-supplied error callback or NULL */
    vnacal_error_fn_t *vcs_error_fn;

    /* user-supplied error callback argument or NULL */
    void *vcs_error_arg;

};

/*
 * VNACAL_CALIBRATION_DATA: return a pointer to vnacal_cdata_t *
 * 	given row and column
 */
#define VNACAL_CALIBRATION_DATA(vcsp, row, column) \
	(&(vcsp)->vcs_matrix[(row) * (vcsp)->vcs_columns + (column)])

/*
 * vnacal_error_terms_t: a cell of the error term matrix
 */
typedef struct vnacal_error_terms {

    /* three column vectors of data (see aliases below) */
    double complex *et_data_vectors[3];

} vnacal_error_terms_t;

/* directivity error (diagonal entry only) */
#define et_e00		et_data_vectors[0]

/* reflection tracking error (diagonal entry only) */
#define et_e10e01	et_data_vectors[1]

/* port 1 match error (diagnonal entry only) */
#define et_e11		et_data_vectors[2]

/* leakage error (off-diagonal entry only) */
#define et_e30		et_data_vectors[0]

/* transmission tracking error (off-diagonal entry only) */
#define et_e10e32	et_data_vectors[1]

/* port 2 match error (off-diagonal entry only) */
#define et_e22		et_data_vectors[2]

/*
 * vnacal_etermset_t: set of calibration data
 *   This structure represents one complete set of VNA calibration
 *   measurements.
 *
 *   A VNA calibration file may contain more than one set of calibration data,
 *   for example, if a relay-controlled attenuator is used and it's desired to
 *   make a separate calibration for each setting.  Another example is where
 *   the reflection bridge can be switched between returning forward /
 *   reflected signal and forward+reflected / forward-reflected signal.
 */
typedef struct vnacal_etermset {
    /* back pointer to vnacal_t object */
    vnacal_t *ets_vcp;

    /* name of this set */
    const char *ets_setname;

    /* number of ports where signal is measured (must be >= columns) */
    int ets_rows;

    /* number of ports where signal is generated */
    int ets_columns;

    /* number of frequencies */
    int ets_frequencies;

    /* vector of frequencies */
    double *ets_frequency_vector;

    /* system impedance of the VNA ports (currently assumed all
       the same) */
    double complex ets_z0;

    /* per-set property list */
    vnaproperty_t *ets_properties;

    /* serialized matrix of error terms */
    vnacal_error_terms_t *ets_error_term_matrix;

} vnacal_etermset_t;

/*
 * VNACAL_ERROR_TERMS: return a pointer to vnacal_error_terms given row
 *   and column
 */
#define VNACAL_ERROR_TERMS(etsp, row, column) \
	(&(etsp)->ets_error_term_matrix[(row) * (etsp)->ets_columns + (column)])

/*
 * vnacal_t: structure returned from vnacal_load
 */
struct vnacal {
    /* calibration filename */
    char *vc_filename;

    /* number of independent sets of calibration data */
    int vc_sets;

    /* precision for frequency values */
    int vc_fprecision;

    /* precision for data values */
    int vc_dprecision;

    /* vector of pointers to vnacal_etermset_t structures, vc_sets long */
    vnacal_etermset_t **vc_set_vector;

    /* user-supplied error callback or NULL */
    vnacal_error_fn_t *vc_error_fn;

    /* user-supplied error callback argument or NULL */
    void *vc_error_arg;

    /* global properties */
    vnaproperty_t *vc_properties;
};

/*
 * vnacal_input_t: measured data from the device under test
 */
struct vnacal_input {
    /* associated vnacal_t object */
    vnacal_t *vi_vcp;

    /* associated calibration set */
    int vi_set;

    /* number of rows */
    int vi_rows;

    /* number of columns */
    int vi_columns;

    /* number of frequencies */
    int vi_frequencies;

    /* vector of frequencies, ascending */
    double *vi_frequency_vector;

    /* true if the frequency vector has been set */
    bool vi_frequencies_valid;

    /* serialized matrix of pointer to vector of values (sum of all adds) */
    double complex **vi_matrix;

    /* serialized matrix of number of vectors added to each cell */
    int *vi_counts;

    /* serialized matrix of VNA matrix cell used to measure each DUT cell */
    int *vi_map;
};

/* report an error */
extern void _vnacal_error(const vnacal_t *vcp, const char *format, ...)
#ifdef __GNUC__
__attribute__ ((__format__ (__printf__, 2, 3)))
#endif /* __GNUC__ */
;

/* _vnacommon_lu: find replace A with its LU decomposition */
extern double complex _vnacommon_lu(complex double *a, int *row_index, int n);

/* _vnacommon_mldivide: find X = A^-1 * B, X mxn, A mxm, B mxn */
double complex _vnacommon_mldivide(complex double *x, complex double *a,
	const double complex *b, int m, int n);

/* _vnacommon_mrdivide: find X = B * A^-1, X mxn, B mxn, A nxn */
extern double complex _vnacommon_mrdivide(double complex *x,
	const double complex *b, double complex *a, int m, int n);

/* _vnacal_calset_get_value: return the requested calibration value */
extern double complex _vnacal_calset_get_value(const vnacal_cdata_t *vcdp,
	int term, int findex);

/* _vnacal_calset_get_reference: get the given reference value */
extern double complex _vnacal_calset_get_reference(
	const vnacal_calset_t *vcsp, int reference, int findex);

/* _vnacal_etermset_alloc: alloc vnacal_etermset_t and internal vectors */
extern vnacal_etermset_t *_vnacal_etermset_alloc(vnacal_t *vcp,
	const char *name, int rows, int columns, int frequencies);

/* _vnacal_etermset_get_fmin_bound: get the lower frequency bound */
extern double _vnacal_etermset_get_fmin_bound(const vnacal_etermset_t *etsp);

/* _vnacal_etermset_get_fmax_bound: get the upper frequency bound */
extern double _vnacal_etermset_get_fmax_bound(const vnacal_etermset_t *etsp);

/* _vnacal_etermset_free: free vnacal_etermset_t and internal vectors */
extern void _vnacal_etermset_free(vnacal_etermset_t *etsp);

/* open a calibration file with given fopen mode */
extern FILE *_vnacal_open(vnacal_t *vcp, const char *pathname,
	const char *dotdir, const char *mode);

/* _vnacal_rfi: apply rational function interpolation */
extern double complex _vnacal_rfi(const double *xp, double complex *yp,
	int n, int m, int *segment, double x);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* _VNACAL_INTERNAL_H */
