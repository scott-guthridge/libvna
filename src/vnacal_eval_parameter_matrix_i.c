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

#include "archdep.h"

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vnacal_internal.h"
#include "vnaconv.h"

/*
 * calc_tline_coefficients0: calc Zc, gl (classic version)
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zc: address of complex to receive the characteristic impedance
 *
 * Returns the transmission coefficient times electrical length (gamma
 * el).  This is the original version described in Keysight note 1287-11:
 * https://people.ece.ubc.ca/robertor/Links_files/Files/AN-1287-11.pdf
 * This form uses an approximation to avoid the need for complex square
 * root.
 */
static double complex calc_tline_coefficients0(const vnacal_calkit_data_t *vcdp,
	double f, double complex *Zc)
{
    double w = 2.0 * M_PI * f;				/* rad/s */
    double fGrt = sqrt(f / 1.0e+9/*Hz*/);		/* unitless */
    double offset_delay = vcdp->vcd_offset_delay;	/* s */
    double offset_loss = vcdp->vcd_offset_loss;		/* Ω/s */
    double offset_z0 = vcdp->vcd_offset_z0;		/* Ω */
    double alpha_l = offset_loss * offset_delay * fGrt /
		     (2.0 * offset_z0);
    double beta_l = w * offset_delay + alpha_l;
    double complex gamma_l = alpha_l + I * beta_l;
    *Zc = offset_z0 + (f != 0.0 ?
	(1.0 - I) * offset_loss * fGrt / (2.0 * w) : 0.0);

    return gamma_l;
}

/*
 * calc_tline_coefficients: calc Z, gl (revised version)
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zc: address of complex to receive the characteristic impedance
 *
 * Returns the transmission coefficient times electrical length
 * (gamma el).  This is the revised version described here:
 * https://www.keysight.com/us/en/assets/7018-01375/application-notes/
 * 5989-4840.pdf
 */
static double complex calc_tline_coefficients(const vnacal_calkit_data_t *vcdp,
	double f, double complex *Zc)
{
    double complex temp;
    double offset_delay = vcdp->vcd_offset_delay;	/* s */
    double offset_loss = vcdp->vcd_offset_loss;		/* Ω/s */
    double offset_z0 = vcdp->vcd_offset_z0;		/* Ω */

    if (f != 0.0) {
	temp = csqrt(1.0 +
		     (1.0 - I) * offset_loss /
		     (2.0 * M_PI * sqrt(1.0e+9 * f) * offset_z0));
    } else {
	temp = 1.0;
    }
    *Zc = offset_z0 * temp;
    return I * 2.0 * M_PI * f * offset_delay * temp;
}

/*
 * eval_calkit_short: evaluate a calkit short standard at given frequency
 *   @vcdp: vnacal_calkit_data_t structure
 *   @z0: the reference impedance
 *   @f: frequency in Hz
 */
static double complex eval_calkit_short(const vnacal_calkit_data_t *vcdp,
    double complex z0, double f)
{
    double L = vcdp->vcd_l_coefficients[0] +
          f * (vcdp->vcd_l_coefficients[1] +
	  f * (vcdp->vcd_l_coefficients[2] +
	  f *  vcdp->vcd_l_coefficients[3]));
    double complex Zl = I * 2.0 * M_PI * f * L;
    double complex Zc, gl, ht, Zi;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = calc_tline_coefficients0(vcdp, f, &Zc);
    } else {
	gl = calc_tline_coefficients(vcdp, f, &Zc);
    }
    ht = ctanh(gl);
    Zi = Zc * (Zl + Zc * ht) / (Zc + Zl * ht);

    return (Zi - conj(z0)) / (Zi + z0);
}

/*
 * eval_calkit_open: evaluate a calkit open standard at given frequency
 *   @vcdp: vnacal_calkit_data_t structure
 *   @z0: the reference impedance
 *   @f: frequency in Hz
 */
static double complex eval_calkit_open(const vnacal_calkit_data_t *vcdp,
    double complex z0, double f)
{
    double C = vcdp->vcd_c_coefficients[0] +
          f * (vcdp->vcd_c_coefficients[1] +
	  f * (vcdp->vcd_c_coefficients[2] +
	  f *  vcdp->vcd_c_coefficients[3]));
    double complex Zl;
    double complex Zc, gl, ht, Zi;

    /*
     * Special-case zero frequency.  The result is 1.0 in the limit
     * regardless of z0.
     */
    if (f == 0.0) {
	return 1.0;
    }

    /*
     * Handle the normal case.
     */
    Zl = 1.0 / (I * 2.0 * M_PI * f * C);
    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = calc_tline_coefficients0(vcdp, f, &Zc);
    } else {
	gl = calc_tline_coefficients(vcdp, f, &Zc);
    }
    ht = ctanh(gl);
    Zi = Zc * (Zl + Zc * ht) / (Zc + Zl * ht);

    return (Zi - conj(z0)) / (Zi + z0);
}

/*
 * eval_calkit_load: evaluate a calkit load standard at given frequency
 *   @vcdp: vnacal_calkit_data_t structure
 *   @z0: the reference impedance
 *   @f: frequency in Hz
 */
static double complex eval_calkit_load(const vnacal_calkit_data_t *vcdp,
    double complex z0, double f)
{
    double complex Zl = vcdp->vcd_zl;
    double complex Zc, gl, ht, Zi;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = calc_tline_coefficients0(vcdp, f, &Zc);
    } else {
	gl = calc_tline_coefficients(vcdp, f, &Zc);
    }
    ht = ctanh(gl);
    Zi = Zc * (Zl + Zc * ht) / (Zc + Zl * ht);

    return (Zi - conj(z0)) / (Zi + z0);
}

/*
 * eval_calkit_through: evaluate a calkit through standard at given frequency
 *   @vcdp: vnacal_calkit_data_t structure
 *   @z0_vector: the reference impedances
 *   @f: frequency in Hz
 *   @result_matrix: 2x2 complex matrix to receive the result
 */
static void eval_calkit_through(const vnacal_calkit_data_t *vcdp,
    const double complex *z0_vector, double f, double complex *result_matrix)
{
    double complex zc, gl, z1, z2, p, p2, mp, pp, c, d;
    double z1r, z2r, rt;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = calc_tline_coefficients0(vcdp, f, &zc);
    } else {
	gl = calc_tline_coefficients(vcdp, f, &zc);
    }

    /*
     * Effectively, we find the ABCD parameters of the transmission line
     * and convert them to S parameters, e.g.:
     *   double complex a[2][2];
     *
     *   a[0][0] = ccosh(gl);
     *   a[0][1] = csinh(gl) * zc;
     *   a[1][0] = csinh(gl) / zc;
     *   a[1][1] = ccosh(gl);
     *   vnaconv_atos(a, (double complex (*)[2])result_matrix, z0_vector);
     *
     * If we convert the trig functions to exponential form, expand the
     * conversion and then refactor, however; we get the more numerically
     * stable form below.
     */
    p = cexp(-gl);
    p2 = p * p;
    pp = 1.0 + p2;
    mp = 1.0 - p2;
    z1 = z0_vector[0];
    z2 = z0_vector[1];
    z1r = creal(z1);
    z2r = creal(z2);
    rt = sqrt(fabs(z1r / z2r));
    d = pp * (z1 + z2) * zc + mp * (z1 * z2 + zc * zc);
    c = 4.0 * p * zc / d;
    result_matrix[0] = ((pp * z2 + mp * zc) * zc -
                        (mp * z2 + pp * zc) * conj(z1)) / d;
    result_matrix[1] = c * z1r / rt;
    result_matrix[2] = c * z2r * rt;
    result_matrix[3] = ((pp * z1 + mp * zc) * zc -
                        (mp * z1 + pp * zc) * conj(z2)) / d;
}

/*
 * eval_data_standard: evaluate a data standard at a given frequency
 *   @function: name of user-called function
 *   @stdp: vnacal_standard_t structure
 *   @segment_ptr: most recent index to optimize _vnacal_rfi (init to 0)
 *   @zr_vector: reference impedances
 *   @frequency: frequency at which to evaluate
 *   @result_matrix: caller-provided buffer to receive the result
 */
static int eval_data_standard(const char *function,
	vnacal_standard_t *stdp, const double complex *zr_vector,
	double frequency, double complex *result_matrix)
{
    vnacal_t *vcp = stdp->std_vcp;
    const int ports = stdp->std_ports;
    vnacal_data_standard_t *vdsp = &stdp->std_data_standard;
    const int frequencies = vdsp->vds_frequencies;
    const double *frequency_vector = vdsp->vds_frequency_vector;
    const double fmin = frequency_vector[0];
    const double fmax = frequency_vector[frequencies - 1];
    double f_lower, f_upper;
    double complex *zd_vector;
    double complex zd_temp[ports];
    int segment = vdsp->vds_segment;

    /*
     * Test if frequency is in bounds.
     */
    f_lower = (1.0 - VNACAL_F_EXTRAPOLATION) * fmin;
    f_upper = (1.0 + VNACAL_F_EXTRAPOLATION) * fmax;
    if (frequency < f_lower || frequency > f_upper) {
	_vnacal_error(vcp, VNAERR_USAGE, "%s: "
		"frequency %e must be between %e and %e for %s standard\n",
		function, frequency, fmin, fmax, stdp->std_name);
	return -1;
    }

    /*
     * Copy the data matrix, interpolating as necessary.
     */
    for (int cell = 0; cell < ports * ports; ++cell) {
	result_matrix[cell] = _vnacal_rfi(frequency_vector,
		vdsp->vds_data[cell], frequencies,
		MIN(frequencies, VNACAL_MAX_M), &segment, frequency);
    }

    /*
     * Find the reference impedances of the data.  If different,
     * renormalize the result matrix.
     */
    if (!vdsp->vds_has_fz0) {
	zd_vector = vdsp->u.vds_z0_vector;
    } else {
	for (int port = 0; port < ports; ++port) {
	    zd_temp[port] = _vnacal_rfi(frequency_vector,
		    vdsp->u.vds_z0_vector_vector[port], frequencies,
		    MIN(frequencies, VNACAL_MAX_M), &segment, frequency);
	}
	zd_vector = zd_temp;
    }
    for (int port = 0; port < ports; ++port) {
	if (cabs(zd_vector[port] - zr_vector[port]) > 1.0e-5) {
	    vnaconv_stosrn(result_matrix, result_matrix,
		    zd_vector, zr_vector, ports);
	    break;
	}
    }
    vdsp->vds_segment = segment;
    return 0;
}

/*
 * vnacal_eval_parameter_matrix: evaluate parameter matrix at given frequency
 *   @function: name of function user called
 *   @vpmmp: vnacal_parameter_matrix_map structure
 *   @frequency: frequency at which to evaluate
 *   @z0_vector: reference impedances results should be returned in
 *   @result_matrix: caller-allocated matrix to hold result
 */
int _vnacal_eval_parameter_matrix_i(const char *function,
        const vnacal_parameter_matrix_map_t *vpmmp, double frequency,
	const double complex *z0_vector, double complex *result_matrix)
{
    int rows = vpmmp->vpmm_rows;
    int columns = vpmmp->vpmm_columns;

    /*
     * Init result matrix to all zeros.
     */
#if BINARY_ZERO_IS_DOUBLE_ZERO
    (void)memset((void *)result_matrix, 0,
	    rows * columns * sizeof(double complex));
#else
    for (int cell = 0; cell < rows * columns; ++cell) {
	result_matrix[cell] = 0.0;
    }
#endif

    /*
     * Evaluate standards.
     */
    assert(vpmmp->vpmm_standard_rmap == NULL || z0_vector != NULL);
    for (const vnacal_standard_rmap_t *vsrmp = vpmmp->vpmm_standard_rmap;
	    vsrmp != NULL; vsrmp = vsrmp->vsrm_next) {
	vnacal_standard_t *stdp = vsrmp->vsrm_stdp;
	const int *port_map = vsrmp->vsrm_rmap_vector;
	const int std_ports = stdp->std_ports;
	double complex std_z0_vector[std_ports];
	double complex std_result_matrix[std_ports * std_ports];

	/*
	 * Fill std_z0_vector.
	 */
	for (int port = 0; port < std_ports; ++port) {
	    std_z0_vector[port] = z0_vector[port_map[port]];
	}

	/*
	 * Evaluate the standard into std_result_matrix.
	 */
	switch (stdp->std_type) {
	case VNACAL_NEW:
	case VNACAL_SCALAR:
	case VNACAL_VECTOR:
	case VNACAL_UNKNOWN:
	case VNACAL_CORRELATED:
	default:
	    abort();

	case VNACAL_CALKIT:
	    switch (stdp->std_calkit_data.vcd_type) {
	    case VNACAL_CALKIT_SHORT:
		std_result_matrix[0] = eval_calkit_short(&stdp->std_calkit_data,
			std_z0_vector[0], frequency);
		break;

	    case VNACAL_CALKIT_OPEN:
		std_result_matrix[0] = eval_calkit_open(&stdp->std_calkit_data,
			std_z0_vector[0], frequency);
		break;

	    case VNACAL_CALKIT_LOAD:
		std_result_matrix[0] = eval_calkit_load(&stdp->std_calkit_data,
			std_z0_vector[0], frequency);
		break;

	    case VNACAL_CALKIT_THROUGH:
		eval_calkit_through(&stdp->std_calkit_data, std_z0_vector,
			frequency, std_result_matrix);
		break;

	    default:
		abort();
	    }
	    break;

	case VNACAL_DATA:
	    if (eval_data_standard(function, stdp, std_z0_vector,
			frequency, std_result_matrix) == -1) {
		return -1;
	    }
	    break;
	}

	/*
	 * Copy the result back to the full matrix being careful to
	 * skip rows or columns of the standard that are missing because
	 * the result matrix is rectangular.
	 */
	for (int std_row = 0; std_row < std_ports; ++std_row) {
	    int row = port_map[std_row];

	    assert(row >= 0);
	    if (row >= rows)
		continue;

	    for (int std_column = 0; std_column < std_ports; ++std_column) {
		int column = port_map[std_column];
		int std_cell;
		int cell;

		assert(column >= 0);
		if (column >= columns)
		    continue;

		std_cell = std_row * std_ports + std_column;
		cell = row * columns + column;
		result_matrix[cell] = std_result_matrix[std_cell];
	    }
	}
    }

    /*
     * Evaluate regular parameters.
     */
    for (const vnacal_parameter_rmap_t *vprmp = vpmmp->vpmm_parameter_rmap;
	    vprmp != NULL; vprmp = vprmp->vprm_next) {
	vnacal_parameter_t *vpmrp = vprmp->vprm_parameter;;
	double complex value;

	switch (vpmrp->vpmr_type) {
	case VNACAL_SCALAR:
	    value = vpmrp->vpmr_coefficient;
	    break;

	case VNACAL_VECTOR:
	case VNACAL_UNKNOWN:	/* always solved values here */
	case VNACAL_CORRELATED:
	    {
		double fmin, fmax;
		double lower, upper;

		assert(vpmrp->vpmr_frequency_vector != NULL);
		fmin = vpmrp->vpmr_frequency_vector[0];
		fmax = vpmrp->vpmr_frequency_vector[vpmrp->vpmr_frequencies-1];
		lower = (1.0 - VNACAL_F_EXTRAPOLATION) * fmin;
		upper = (1.0 + VNACAL_F_EXTRAPOLATION) * fmax;
		if (frequency < lower || frequency > upper) {
		    _vnacal_error(vpmmp->vpmm_vcp, VNAERR_USAGE,
			    "%s: frequency %e must be between %e and %e\n",
			    function, frequency, fmin, fmax);
		    return -1;
		}
		value = _vnacal_rfi(vpmrp->vpmr_frequency_vector,
			vpmrp->vpmr_coefficient_vector,
			vpmrp->vpmr_frequencies,
			MIN(vpmrp->vpmr_frequencies, VNACAL_MAX_M),
			&vpmrp->vpmr_segment,
			frequency);
	    }
	    break;

	default:
	    abort();
	}
	result_matrix[vprmp->vprm_cell] = value;
    }
    return 0;
}
