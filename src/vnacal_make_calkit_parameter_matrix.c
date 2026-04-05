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
 * add_tline_from_zl: return s11 from impedance at end of transmission line
 *   @vcdp: vnacal_calkit_data_t structure
 *   @z0: the reference impedance
 *   @f: frequency in Hz
 *   @zl: load impedance
 */
static double complex add_tline_from_zl(const vnacal_calkit_data_t *vcdp,
	double complex z0, double f, double complex zl)
{
    double complex zc, gl;
    double complex e, tanh_num, tanh_den;
    double complex num, den;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = calc_tline_coefficients0(vcdp, f, &zc);
    } else {
	gl = calc_tline_coefficients(vcdp, f, &zc);
    }

    /*
     * The input impedance of a transmission line terminated in zl is:
     *   zi = zc * (zl + zc * tanh(gl)) / (zc + zl * tanh(gl)).
     *
     * But tanh is infinite at quarter wavelength delays.  To avoid
     * infinity, use the substitution:
     *
     *     tanh(gl) = (exp(2 gl) - 1) / (exp(2 gl) + 1)
     */
    e = cexp(2.0 * gl);
    tanh_num = e - 1.0;
    tanh_den = e + 1.0;

    /*
     * Compute s11 = (zi - conj(z0)) / (zi + z0) with zi expressed
     * in terms of tanh_num and tanh_den with inner fractions removed.
     */
    num = zc * (zc * tanh_num + zl * tanh_den)
          - conj(z0) * (zc * tanh_den + zl * tanh_num);
    den = zc * (z0 + zl) * tanh_den
          + (zc * zc + z0 * zl) * tanh_num;

    return num / den;
}

/*
 * add_tline_from_yl: return s11 from admittance at end of transmission line
 *   @vcdp: vnacal_calkit_data_t structure
 *   @z0: the reference impedance
 *   @f: frequency in Hz
 *   @yl: load admittance
 */
static double complex add_tline_from_yl(const vnacal_calkit_data_t *vcdp,
	double complex z0, double f, double complex yl)
{
    double complex zc, gl;
    double complex e, tanh_num, tanh_den;
    double complex num, den;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = calc_tline_coefficients0(vcdp, f, &zc);
    } else {
	gl = calc_tline_coefficients(vcdp, f, &zc);
    }

    /*
     * Use tanh(gl) = (exp(2 gl) - 1) / (exp(2 gl) + 1) to avoid
     * infinity at quarter wavelengths.
     */
    e = cexp(2.0 * gl);
    tanh_num = e - 1.0;
    tanh_den = e + 1.0;

    /*
     * Compute s11 = (zi - conj(z0)) / (zi + z0) with zi expressed
     * in terms of yl, tanh_num and tanh_den with inner fractions
     * removed.
     */
    num = zc * (zc * yl * tanh_num + tanh_den)
	  - conj(z0) * (zc * yl * tanh_den + tanh_num);
    den = zc * (yl * z0 + 1.0) * tanh_den
	  + (zc * zc * yl + z0) * tanh_num;

    return num / den;
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
    double complex zl = I * 2.0 * M_PI * f * L;

    return add_tline_from_zl(vcdp, z0, f, zl);
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
    double complex yl = I * 2.0 * M_PI * f * C;

    return add_tline_from_yl(vcdp, z0, f, yl);
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
    return add_tline_from_zl(vcdp, z0, f, vcdp->vcd_zl);
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
 * eval_calkit_standard: evaluate the standard into result_matrix
 *   @stdp: vnacal_standard_t structure
 *   @function: name of user-called function
 *   @z0_vector: reference impedance vector
 *   @frequency: frequency (in Hz) to evaluate
 *   @result_matrix: caller-allocated matrix to hold result
 */
static int eval_calkit_standard(vnacal_standard_t *stdp, const char *function,
	const double complex *z0_vector, double frequency,
	double complex *result_matrix)
{
    vnacal_calkit_standard_t *cstdp;

    assert(stdp->std_ops->stdo_type == VNACAL_CALKIT);
    cstdp = (vnacal_calkit_standard_t *)stdp;
    switch (cstdp->cstd_calkit_data.vcd_type) {
    case VNACAL_CALKIT_SHORT:
	result_matrix[0] = eval_calkit_short(&cstdp->cstd_calkit_data,
		z0_vector[0], frequency);
	break;

    case VNACAL_CALKIT_OPEN:
	result_matrix[0] = eval_calkit_open(&cstdp->cstd_calkit_data,
		z0_vector[0], frequency);
	break;

    case VNACAL_CALKIT_LOAD:
	result_matrix[0] = eval_calkit_load(&cstdp->cstd_calkit_data,
		z0_vector[0], frequency);
	break;

    case VNACAL_CALKIT_THROUGH:
	eval_calkit_through(&cstdp->cstd_calkit_data, z0_vector,
		frequency, result_matrix);
	break;

    default:
	abort();
    }
    return 0;
}

/*
 * _calkit_ops: subclass operations on vnacal_calkit_standard_t
 */
static const vnacal_standard_ops_t _calkit_ops = {
    .stdo_type = VNACAL_CALKIT,
    .stdo_eval = eval_calkit_standard,
    .stdo_free = NULL
};

/*
 * vnacal_make_calkit_parameter_matrix: make parameter matrix for kit standard
 *   @function: name of user-called function
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
static int _vnacal_make_calkit_parameter_matrix(const char *function,
	vnacal_t *vcp, const vnacal_calkit_data_t *vcdp,
	int *parameter_matrix, size_t parameter_matrix_size)
{
    vnacal_calkit_standard_t *cstdp = NULL;
    vnacal_standard_t *stdp = NULL;
    const char *name;
    int ports = 0;

    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (vcdp == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: vcdp cannot be NULL", function);
	return -1;
    }
    if ((name = _vnacal_get_calkit_name(vcdp, &ports)) == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: vnacal_calkit_data_t structure is not valid", function);
	return -1;
    }
    if (ports * ports * sizeof(int) > parameter_matrix_size) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: insufficient result matrix allocation", function);
	return -1;
    }
    _vnacal_init_parameter_matrix(parameter_matrix, ports, ports);

    /*
     * Allocate and init the standard.
     */
    if ((cstdp = _vnacal_alloc_standard(function, vcp, &_calkit_ops,
		    ports, sizeof(vnacal_calkit_standard_t))) == NULL) {
	goto error;
    }
    stdp = &cstdp->cstd_base;
    if ((stdp->std_name = strdup(name)) == NULL) {
	goto error;
    }
    cstdp->cstd_calkit_data = *vcdp;

    /*
     * Add the parameter structures.
     */
    if (_vnacal_fill_standard_parameter_matrix(function, stdp,
		parameter_matrix) == -1) {
	goto error;
    }
    _vnacal_release_standard(&stdp);	/* release initial reference */
    assert(stdp != NULL);
    return ports;

error:
    if (stdp != NULL) {
	_vnacal_release_standard(&stdp); /* release initial reference */
	assert(stdp == NULL);
    }
    return -1;
}

/*
 * vnacal_make_calkit_parameter: make a parameter for a one-port kit standard
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @vcdp: data describing the calkit standard
 */
int vnacal_make_calkit_parameter(vnacal_t *vcp,
	const vnacal_calkit_data_t *vcdp)
{
    int parameter;

    if (_vnacal_make_calkit_parameter_matrix(__func__, vcp, vcdp,
	    &parameter, sizeof(parameter)) == -1) {
	return -1;
    }
    return parameter;
}

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
int vnacal_make_calkit_parameter_matrix(vnacal_t *vcp,
	const vnacal_calkit_data_t *vcdp, int *parameter_matrix,
	size_t parameter_matrix_size)
{
    return _vnacal_make_calkit_parameter_matrix(__func__, vcp, vcdp,
	    parameter_matrix, parameter_matrix_size);
}
