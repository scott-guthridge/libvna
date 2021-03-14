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

/*
 * Example of "through", "short", "delay" (TSD) calibration in 10-term
 * T and E parameters.  TSD calibration closely resembles "through",
 * "reflect", "line" (TRL) calibration except that the reflect and line
 * standards must be fully known.
 */

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vnacal.h>
#include <vnaconv.h>
#include <vnafile.h>

char *progname;

/*
 * PI, Z0: misc constants
 *   PI: used below to convert from Hz to angular frequency
 *   C:  speed of light in vacuum
 *   Z0: system impedance
 *   FC: center frequency of the line standard
 */
#define PI		3.14159265358979
#define Z0		50.0			/* ohms */
#define FC		18e+9			/* Hz */

/*
 * ACTUAL_FILE: file containing the actual DUT s-parameters
 */
#define ACTUAL_FILE	"MwT-1F.s2p"

/*
 * C_FMIN, C_FMAX, C_FREQUENCIES: cal frequency range and number of points
 *
 * Highest calibration frequency over lowest is restricted to a factor
 * of 8 so that the phase shift of the delay standard can remain within
 * a range of 20..160 degrees.
 */
#define C_FMIN		 4.0e+9		/* Hz */
#define C_FMAX		32.0e+9
#define C_FREQUENCIES	50

/*
 * row_type: a row of a 2x2 complex matrix
 *     This is used to simplify parameter declarations below.
 */
typedef double complex row_type[2];

/*
 * multiply: multiply matrices a and b, producing c
 *   @a: first matrix
 *   @b: second matrix
 *   @c: result matrix
 */
static void multiply(const row_type *a, const row_type *b, row_type *c)
{
    for (int i = 0; i < 2; ++i) {
	for (int k = 0; k < 2; ++k) {
	    double complex sum = 0.0;

	    for (int j = 0; j < 2; ++j) {
		sum += a[i][j] * b[j][k];
	    }
	    c[i][k] = sum;
	}
    }
}

/*
 * z0: system impedances for vnaconv_*
 */
static const double complex z0[2] = { Z0, Z0 };

/*
 * C1: shunt capacitance at VNA port 1
 * L2: series inductance at VNA port 2
 *
 * These values are used below to introduce errors to our VNA.
 */
#define C1		265.258e-15		/* 265 femtofarads */
#define L2		663.146e-12		/* 663 picohenries */

/*
 * measurement_t: selcts which measurement vna_measure returns
 */
typedef enum measurement {
    MEASURE_THROUGH,
    MEASURE_SHORT,
    MEASURE_DELAY,
    MEASURE_DUT
} measurement_t;

/*
 * vna_measure: simulate the VNA making the requested measurement
 *   @measurement: which measurement to simulate
 *   @frequencies: number of frequency points to return
 *   @frequency_vector: optional vector to receive frequencies
 *   @a_result: 2x2 matrix of vectors to receive forward voltages
 *   @b_result: 2x2 matrix of vectors to receive reflected voltages
 *   @delay_abcd: per-frequency vector of ABCD matrices describing delay
 *   @dut_actual: actual DUT s-parameters (needed only for MEASURE_DUT)
 *
 * The VNA measurements have the form 2x2 matrix of pointers to vectors
 * of values, one per frequency.  Conversely, delay_abcd and dut_actual
 * are vectors, one entry per frequency, of 2x2 matrices.  Be careful
 * of the distinction.
 */
static int vna_measure(measurement_t measurement,
	int frequencies, double *frequency_vector,
	double complex *(*a_result)[2], double complex *(*b_result)[2],
	double complex (*delay_abcd)[2][2], const vnadata_t *dut_actual)
{
    /*
     * Validate parameters.
     */
    switch (measurement) {
    case MEASURE_DELAY:
	assert(delay_abcd != NULL);
	break;

    case MEASURE_DUT:
	assert(dut_actual != NULL);
	if (frequencies > vnadata_get_frequencies(dut_actual)) {
	    (void)fprintf(stderr, "%s: vna_measure: %d frequencies "
		    "requested, but only %d known for DUT\n",
		    progname, frequencies,
		    vnadata_get_frequencies(dut_actual));
	    exit(2);
	}
	break;

    default:
	break;
    }

    /*
     * For each frequency...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f;
	double complex s;
	double complex a[2][2], b[2][2];
	double complex port1_abcd[2][2], port2_abcd[2][2];
	double complex temp1[2][2], temp2[2][2];

	/*
	 * Find the frequency.	For calibration measurements, space
	 * the frequencies linearly between C_FMIN and C_FMAX; for DUT
	 * measurements, take them from dut_actual.
	 */
	if (measurement != MEASURE_DUT) {
	    f = C_FMIN + (C_FMAX - C_FMIN) * (double)findex /
		(frequencies - 1);
	} else {
	    f = vnadata_get_frequency(dut_actual, findex);
	}
	s = 2.0 * PI * I * f;

	/*
	 * For all measurements, fill in the "a" matrix to simulate
	 * leakage in the VNA switch.  Send 2/3 of the signal to the
	 * intended port and 1/3 to the other.  Despite this simple
	 * example, the a matrix can be any arbitrary non-singular
	 * complex matrix.
	 */
	a[0][0] = 2.0 / 3.0;
	a[0][1] = 1.0 / 3.0;
	a[1][0] = 1.0 / 3.0;
	a[1][1] = 2.0 / 3.0;

	/*
	 * Fill in ABCD parameters representing the errors at VNA port 1.
	 * Here, we just add a shunt capacitance of C1.
	 */
	port1_abcd[0][0] = 1.0;
	port1_abcd[0][1] = 0.0;
	port1_abcd[1][0] = s * C1;
	port1_abcd[1][1] = 1.0;

	/*
	 * Fill in ABCD parametes representing the errors at VNA port 2.
	 * Here, we just add a series inductance of L2.
	 */
	port2_abcd[0][0] = 1.0;
	port2_abcd[0][1] = s * L2;
	port2_abcd[1][0] = 0.0;
	port2_abcd[1][1] = 1.0;

	/*
	 * Calculate b for the requested measurement.
	 */
	switch (measurement) {
	case MEASURE_THROUGH:
	    /*
	     * Multiply the ABCD parameters of the two error boxes,
	     * convert to s-parameters and find b = s a.
	     */
	    multiply(port1_abcd, port2_abcd, temp1);
	    vnaconv_atos(temp1, temp1, z0);
	    multiply(temp1, a, b);
	    break;

	case MEASURE_SHORT:
	    /*
	     * The reflection coefficient looking into an error box
	     * in ABCD parameters with the other port shorted has this
	     * simple form.
	     */
	    temp1[0][0] = (port1_abcd[0][1] - port1_abcd[1][1] * Z0) /
	                  (port1_abcd[0][1] + port1_abcd[1][1] * Z0);
	    temp1[0][1] = 0.0;
	    temp1[1][0] = 0.0;
	    temp1[1][1] = (port2_abcd[0][1] - port2_abcd[0][0] * Z0) /
	                  (port2_abcd[0][1] + port2_abcd[0][0] * Z0);
	    multiply(temp1, a, b);
	    break;

	case MEASURE_DELAY:
	    /*
	     * Multiply the ABCD parameters of the first error box, the
	     * delay and the second error box, convert to s-parameters
	     * and find b = s a.
	     */
	    multiply(port1_abcd, delay_abcd[findex], temp2);
	    multiply(temp2, port2_abcd, temp1);
	    vnaconv_atos(temp1, temp1, z0);
	    multiply(temp1, a, b);
	    break;

	case MEASURE_DUT:
	    /*
	     * Convert the actual s-parameters of the DUT to ABCD
	     * parameters, multiply the ABCD parameters of the first
	     * error, the DUT and the second error box, convert to
	     * s-parameters and find b = s a.
	     */
	    {
		const double complex *dut_s;
		double complex dut_abcd[2][2];

		dut_s = vnadata_get_matrix(dut_actual, findex);
		assert(dut_s != NULL);
		vnaconv_stoa((const row_type *)dut_s, dut_abcd, z0);
		multiply(port1_abcd, dut_abcd, temp2);
		multiply(temp2, port2_abcd, temp1);
		vnaconv_atos(temp1, temp1, z0);
		multiply(temp1, a, b);
	    }
	    break;

	default:
	    abort();
	}

	/*
	 * Copy the results to the caller's arrays.
	 */
	if (frequency_vector != NULL) {
	    frequency_vector[findex] = f;
	}
	for (int row = 0; row < 2; ++row) {
	    for (int column = 0; column < 2; ++column) {
		a_result[row][column][findex] = a[row][column];
	    }
	}
	for (int row = 0; row < 2; ++row) {
	    for (int column = 0; column < 2; ++column) {
		b_result[row][column][findex] = b[row][column];
	    }
	}
    }
    return 0;
}

/*
 * get_delay_parameters: find ABCD and S parameters of the delay standard
 *   @vcp: pointer returned from vnacal_create
 *   @frequencies: length of frequency_vector
 *   @frequency_vector: vector of calibration frequencies
 *   @delay_abcd: caller allocated vector of matrices to receive A parameters
 *   @delay_s: caller allocated matrix of int to receive s parameter indices
 *
 * Notes
 *   - The returned delay_abcd vector is used only by our simulated VNA;
 *     if using a real VNA, we wouldn't need to return this.
 *
 *   - The returned delay_s matrix doesn't actually contain s-parameters
 *     but rather integer indices into parameters we store into the
 *     vnacal_t structure.  The reason for this extra layer of indirection
 *     is to support unknown parameters that must be found by the library.
 */
static void get_delay_parameters(vnacal_t *vcp,
	int frequencies, const double *frequency_vector,
	double complex (*delay_abcd)[2][2], int (*delay_s)[2])
{

    double complex s[2][2][frequencies];	/* temp matrix of vectors */

    /*
     * Find the ABCD and S parameters of the delay standard at each
     * frequency.  For simplicity, we assume a lossless delay; however,
     * gl can be given a real component to represent a lossy delay.
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = frequency_vector[findex];
	double complex gl = I * PI * f / FC;
	double complex temp_s[2][2];

	delay_abcd[findex][0][0] = ccosh(gl);
	delay_abcd[findex][0][1] = csinh(gl) * Z0;
	delay_abcd[findex][1][0] = csinh(gl) / Z0;
	delay_abcd[findex][1][1] = ccosh(gl);

	/*
	 * Convert to S and store as matrix of vectors.
	 */
	vnaconv_atos(delay_abcd[findex], temp_s, z0);
	for (int row = 0; row < 2; ++row) {
	    for (int column = 0; column < 2; ++column) {
		s[row][column][findex] = temp_s[row][column];
	    }
	}
    }

    /*
     * Use vnacal_make_vector_parameter to store the s-parameters into
     * the vnacal_t structure and fill delay_s with the returned integer
     * indices that represent them.
     */
    for (int row = 0; row < 2; ++row) {
	for (int column = 0; column < 2; ++column) {
	    if ((delay_s[row][column] = vnacal_make_vector_parameter(vcp,
			    frequencies, frequency_vector,
			    s[row][column])) == -1) {
		exit(3);
	    }
	}
    }
}

/*
 * delete_delay_s: delete the delay_s parameters from the vnacal_t structure
 *   @vcp: pointer returned from vnacal_create
 *   @delay_s: matrix of parameter indices
 */
static void delete_delay_s(vnacal_t *vcp, int (*delay_s)[2])
{
    for (int row = 0; row < 2; ++row) {
	for (int column = 0; column < 2; ++column) {
	    vnacal_delete_parameter(vcp, delay_s[row][column]);
	    delay_s[row][column] = -1;
	}
    }
}

/*
 * error_fn: print errors from the vnacal library
 *   @category: category of error (ignored here)
 *   @message: single line error message without a newline
 *   @error_arg: passed through to the error function (unused here)
 */
static void error_fn(vnaerr_category_t category, const char *message,
	void *error_arg)
{
    (void)fprintf(stderr, "example: %s\n", message);
}

/*
 * make_calibration: make a calibration file for the simulated VNA
 */
static void make_calibration()
{
    vnacal_t *vcp;
    vnacal_new_t *vnp;
    double frequency_vector[C_FREQUENCIES];
    double complex a11_vector[C_FREQUENCIES];
    double complex a12_vector[C_FREQUENCIES];
    double complex a21_vector[C_FREQUENCIES];
    double complex a22_vector[C_FREQUENCIES];
    double complex b11_vector[C_FREQUENCIES];
    double complex b12_vector[C_FREQUENCIES];
    double complex b21_vector[C_FREQUENCIES];
    double complex b22_vector[C_FREQUENCIES];
    double complex *a[2][2] = {
	{ a11_vector, a12_vector },
	{ a21_vector, a22_vector }
    };
    double complex *b[2][2] = {
	{ b11_vector, b12_vector },
	{ b21_vector, b22_vector }
    };
    double complex delay_abcd[C_FREQUENCIES][2][2];
    int delay_s[2][2];

    /*
     * Create the calibration container structure.
     */
    if ((vcp = vnacal_create(error_fn, /*error_arg=*/NULL)) == NULL) {
	exit(4);
    }

    /*
     * Start a new calibration.
     */
    if ((vnp = vnacal_new_alloc(vcp, VNACAL_TE10, 2, 2,
		    C_FREQUENCIES)) == NULL) {
	exit(5);
    }

    /*
     * Make the calibration measurements for through, short and delay
     * standards.  Normally, we would interact with the user between
     * each of these steps to get the user to connect each standard
     * in sequence.  In our simulated environment, we skip this.
     * The frequency m_vector is filled from the first measurement only;
     * the frequencies for the other calibration steps have to be the
     * same as the first.
     */

    /*
     * Add through standard and set frequency vector.
     */
    vna_measure(MEASURE_THROUGH, C_FREQUENCIES, frequency_vector,
	    a, b, /*delay_abcd*/NULL, /*vdp_actual*/NULL);
    vnacal_new_set_frequency_vector(vnp, frequency_vector);
    vnacal_new_add_through(vnp, &a[0][0], 2, 2, &b[0][0], 2, 2,
	    /*port1*/1, /*port2*/2);

    /*
     * Add short standard.
     */
    vna_measure(MEASURE_SHORT, C_FREQUENCIES, /*frequency_vector*/NULL,
	    a, b, /*delay_abcd*/NULL, /*vdp_actual*/NULL);
    vnacal_new_add_double_reflect(vnp, &a[0][0], 2, 2, &b[0][0], 2, 2,
	    VNACAL_SHORT, VNACAL_SHORT, 1, 2);

    /*
     * Add delay standard.
     */
    get_delay_parameters(vcp, C_FREQUENCIES, frequency_vector,
	    delay_abcd, delay_s);
    vna_measure(MEASURE_DELAY, C_FREQUENCIES, /*frequency_vector*/NULL,
	    a, b, delay_abcd, /*vdp_actual*/NULL);
    vnacal_new_add_line(vnp, &a[0][0], 2, 2, &b[0][0], 2, 2, &delay_s[0][0],
	    /*port1*/1, /*port2*/2);
    delete_delay_s(vcp, delay_s);

    /*
     * Solve for the error terms.
     */
    if (vnacal_new_solve(vnp) == -1) {
	exit(6);
    }

    /*
     * Add the new calibration to the vancal_t structure and save.
     */
    if (vnacal_add_calibration(vcp, "cal-T8", vnp) == -1) {
	exit(7);
    }
    if (vnacal_save(vcp, "TSD.vnacal") == -1) {
	exit(8);
    }
    vnacal_new_free(vnp);
    vnacal_free(vcp);
    vcp = NULL;
}

/*
 * apply_calibration: apply the calibration to the DUT
 *     Normally, make_calibration and apply_calibration would be in
 *     separate programs, but to keep the example simple, we've just
 *     made them separate functions.
 */
static void apply_calibration()
{
    vnacal_t *vcp = NULL;
    vnafile_t *vfp = NULL;
    vnadata_t *vdp_actual = NULL;
    vnadata_t *vdp_corrected = NULL;
    int frequencies = 0;
    double frequency_vector[C_FREQUENCIES];
    double complex a11_vector[C_FREQUENCIES];
    double complex a12_vector[C_FREQUENCIES];
    double complex a21_vector[C_FREQUENCIES];
    double complex a22_vector[C_FREQUENCIES];
    double complex b11_vector[C_FREQUENCIES];
    double complex b12_vector[C_FREQUENCIES];
    double complex b21_vector[C_FREQUENCIES];
    double complex b22_vector[C_FREQUENCIES];
    double complex *a[2][2] = {
	{ a11_vector, a12_vector },
	{ a21_vector, a22_vector }
    };
    double complex *b[2][2] = {
	{ b11_vector, b12_vector },
	{ b21_vector, b22_vector }
    };

    /*
     * Load the calibration file.
     */
    if ((vcp = vnacal_load("TSD.vnacal", error_fn,
		    /*error_arg=*/NULL)) == NULL) {
	exit(9);
    }

    /*
     * Allocate a vnadata_t structure to hold the actual s-parameters
     * of the DUT.
     */
    if ((vdp_actual = vnadata_alloc()) == NULL) {
	(void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
		progname, strerror(errno));
	exit(10);
    }

    /*
     * Load the actual S-parameters.
     */
    if ((vfp = vnafile_load(ACTUAL_FILE, VNAFILE_AUTO,
		    error_fn, /*error_arg*/NULL, vdp_actual)) == NULL) {
	exit(11);
    }

    /*
     * Convert to S parameters if not already S.
     */
    if (vnadata_convert(vdp_actual, vdp_actual, VPT_S) == -1) {
	(void)fprintf(stderr, "%s: vnadata_convert: %s: %s\n",
		progname, ACTUAL_FILE, strerror(errno));
	exit(12);
    }
    vnafile_free(vfp);
    vfp = NULL;

    /*
     * Use our simulated VNA to measure the DUT with errors.
     */
    frequencies = vnadata_get_frequencies(vdp_actual);
    vna_measure(MEASURE_DUT, frequencies, frequency_vector,
	    a, b, /*delay_abcd*/NULL, vdp_actual);

    /*
     * Allocate a vnadata_t structure to receive the computed S parameters.
     */
    if ((vdp_corrected = vnadata_alloc()) == NULL) {
	(void)fprintf(stderr, "example: vnadata_alloc: %s\n",
		strerror(errno));
	exit(13);
    }

    /*
     * Apply the calibration and report the corrected values.
     */
    if (vnacal_apply(vcp, /*index*/0, frequency_vector, frequencies,
		&a[0][0], 2, 2, &b[0][0], 2, 2, vdp_corrected) == -1) {
	exit(14);
    }

    /*
     * Print the actual S-parameters from the device under test.
     */
    (void)printf("# actual\n");
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = vnadata_get_frequency(vdp_actual, findex);
	double complex s11 = vnadata_get_cell(vdp_actual, findex, 0, 0);
	double complex s12 = vnadata_get_cell(vdp_actual, findex, 0, 1);
	double complex s21 = vnadata_get_cell(vdp_actual, findex, 1, 0);
	double complex s22 = vnadata_get_cell(vdp_actual, findex, 1, 1);

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		f,
		creal(s11), cimag(s11),
		creal(s12), cimag(s12),
		creal(s21), cimag(s21),
		creal(s22), cimag(s22));
    }
    (void)printf("\n\n");

    /*
     * Print the "b" values as measured from the imperfect VNA.
     */
    (void)printf("# measured\n");
    for (int findex = 0; findex < frequencies; ++findex) {
	double complex b11 = b[0][0][findex];
	double complex b12 = b[0][1][findex];
	double complex b21 = b[1][0][findex];
	double complex b22 = b[1][1][findex];

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[findex],
		creal(b11), cimag(b11),
		creal(b12), cimag(b12),
		creal(b21), cimag(b21),
		creal(b22), cimag(b22));
    }
    (void)printf("\n\n");

    /*
     * Print the corrected values.
     */
    (void)printf("# corrected\n");
    for (int findex = 0; findex < frequencies; ++findex) {
	double complex s11 = vnadata_get_cell(vdp_corrected, findex, 0, 0);
	double complex s12 = vnadata_get_cell(vdp_corrected, findex, 0, 1);
	double complex s21 = vnadata_get_cell(vdp_corrected, findex, 1, 0);
	double complex s22 = vnadata_get_cell(vdp_corrected, findex, 1, 1);

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[findex],
		creal(s11), cimag(s11),
		creal(s12), cimag(s12),
		creal(s21), cimag(s21),
		creal(s22), cimag(s22));
    }

    /*
     * Free memory.
     */
    vnadata_free(vdp_corrected);
    vnadata_free(vdp_actual);
    vnacal_free(vcp);
}

/*
 * main
 */
int main(int argc, char **argv)
{
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    make_calibration();
    apply_calibration();
    exit(0);
}
