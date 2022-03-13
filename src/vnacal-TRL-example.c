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

/*
 * Example of "through", "reflect", "line" (TRL) calibration in 10-term T
 * and E parameters, where the reflection parameter and line parameter
 * are only partially known.
 */

#include <assert.h>
#include <complex.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vnacal.h>
#include <vnaconv.h>
#include <vnadata.h>

char *progname;

/*
 * I, C, Z0, ER_EFF: misc constants
 *   I: used below to convert from Hz to angular frequency
 *   C:  speed of light in vacuum
 *   Z0: system impedance
 *   ER_EFF: effective permittivity
 */
#define PI		3.14159265358979
#define C		2.9979246e+08		/* m/s */
#define Z0		50.0			/* ohms */
#define NP_PER_DB	0.11512925		/* neper per dB */
#define MM_PER_M	1000.0			/* mm per meter */
#define ER_EFF		8.25

/*
 * C_FMIN, C_FMAX, C_FREQUENCIES: cal frequency range and number of points
 */
#define C_FMIN		 1.0e+9			/* Hz */
#define C_FMAX		 8.0e+9
#define C_FREQUENCIES	50

/*
 * M_FREQUENCIES: number of actual DUT frequencies
 */
#define M_FREQUENCIES	339

/*
 * ACTUAL_FILE: file containing the actual DUT s-parameters
 */
#define ACTUAL_FILE	"BFCV-4085+_Plus25DegC.s2p"

/*
 * VNA Port 1 parastic elements
 *   From the directional coupler, L1 and R1 are in series
 *   and C1 is shunted across the port.
 */
#define R1		10.0			/* ohms */
#define L1		3.979e-9		/* henries */
#define C1		1.592e-12		/* farads */

/*
 * VNA Port 2 parasitic elements
 *   From the directional coupler, L2 and C2 are in series
 *   and R2 is shunted across the port.
 */
#define R2		100.0			/* ohms */
#define L2		1.326e-9		/* henries */
#define C2		530.5e-15		/* farads */

/*
 * Errors in the Reflect Standard
 *    Resistor RR in series with inductor, RL
 */
#define RR		5.0			/* ohms */
#define RL		707.4e-12		/* henries */

/*
 * Errors in the Line Standard
 */
#define LINE_LOSS	0.5			/* dB/mm */
#define PHASE_ERROR	10.			/* degrees */

/*
 * Calculate length of the line standard in meters.
 */
#define FC		((C_FMIN + C_FMAX) / 2.0)
#define KAPPA		(1.0 / sqrt(ER_EFF))
#define LINE_LENGTH	(0.25 * C / FC * KAPPA)

/*
 * Calculate ideal and actual line gamma in meters^-1.
 */
#define IDEAL_GAMMA(f) \
	(I * 2.0 * PI * (f) / (C * KAPPA))
#define ACTUAL_GAMMA(f) \
	(IDEAL_GAMMA(f) * cexp(I * PI * PHASE_ERROR / 180.0) + \
	 NP_PER_DB * MM_PER_M * LINE_LOSS)

/*
 * row_type: a row of a 2x2 complex matrix
 */
typedef double complex row_type[2];

/*
 * z0: system impedances for vnaconv_*
 */
static const double complex z0[2] = { Z0, Z0 };

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
 * error_fn: print errors from the vnacal library
 *   @message: single line error message without a newline
 *   @error_arg: passed through to the error function (unused here)
 *   @category: category of error (ignored here)
 */
static void error_fn(const char *message, void *error_arg,
	vnaerr_category_t category)
{
    (void)fprintf(stderr, "example: %s\n", message);
}

/*
 * measurement_type: selcts which measurement vna_measure returns
 */
typedef enum measurement {
    MEASURE_THROUGH,
    MEASURE_REFLECT,
    MEASURE_LINE,
    MEASURE_DUT
} measurement_type;

/*
 * dut_info_type: information on the actual DUT used by vna_measure
 */
typedef struct dut_info {
    vnadata_t *dut_vdp;		/* actual S parameters of DUT */
    int        dut_offset;	/* first offset in frequency range */
    int        dut_frequencies;	/* number of frequencies in range */
} dut_info_type;

/*
 * dut_setup: set up the simulated VNA to measure the device under test
 */
static void dut_setup(dut_info_type *dutp)
{
    vnadata_t *vdp_actual = NULL;
    int frequencies;

    /*
     * Load the actual S-parameters.  Note that this is just for the
     * simulated VNA.  Normally, we wouldn't know these.
     */
    if ((vdp_actual = vnadata_alloc(error_fn, /*error_arg*/NULL)) == NULL) {
	exit(20);
    }
    if (vnadata_load(vdp_actual, ACTUAL_FILE) == -1) {
	exit(21);
    }
    if (vnadata_convert(vdp_actual, vdp_actual, VPT_S) == -1) {
	(void)fprintf(stderr, "%s: vnadata_convert: %s: %s\n",
		progname, ACTUAL_FILE, strerror(errno));
	exit(22);
    }

    /*
     * Find the subrange of frequencies within C_FMIN..C_FMAX and fill
     * in the dut_info structure.
     */
    (void)memset((void *)dutp, 0, sizeof(*dutp));
    dutp->dut_vdp = vdp_actual;
    frequencies = vnadata_get_frequencies(vdp_actual);
    {
	int i = 0;
	int offset;

	while (i < frequencies &&
		vnadata_get_frequency(vdp_actual, i) < C_FMIN)
	    ++i;
	dutp->dut_offset = offset = i;
	while (i < frequencies &&
		vnadata_get_frequency(vdp_actual, i) <= C_FMAX)
	    ++i;
	frequencies = i - dutp->dut_offset;
	if (frequencies > M_FREQUENCIES) {
	    frequencies = M_FREQUENCIES;
	}
	dutp->dut_frequencies = frequencies;
    }
}

/*
 * Make short aliases for the matrix elements below.
 */
#define A11	a[0][0]
#define A12	a[0][1]
#define A21	a[1][0]
#define A22	a[1][1]
#define U11	u[0][0]
#define U12	u[0][1]
#define U21	u[1][0]
#define U22	u[1][1]
#define V11	v[0][0]
#define V12	v[0][1]
#define V21	v[1][0]
#define V22	v[1][1]
#define W11	w[0][0]
#define W12	w[0][1]
#define W21	w[1][0]
#define W22	w[1][1]

/*
 * vna_measure: simulate the VNA making the requested measurement
 *   @dutp: info on the actual DUT (for MEASURE_DUT only)
 *   @measurement: which measurement to simulate
 *   @frequency_vector: optional vector to receive frequencies
 *   @a_result: 2x2 matrix of vectors to receive forward voltages
 *   @b_result: 2x2 matrix of vectors to receive reflected voltages
 */
static int vna_measure(dut_info_type *dutp,
	measurement_type measurement, double *frequency_vector,
	double complex *(*a_result)[2], double complex *(*b_result)[2])
{
    int frequencies;

    /*
     * Determine the number of frequencies.
     */
    if (measurement == MEASURE_DUT) {
	frequencies = dutp->dut_frequencies;
    } else {
	frequencies = C_FREQUENCIES;
    }

    /*
     * For each frequency...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f;
	double complex s;
	double complex port1_abcd[2][2], port2_abcd[2][2];
	double complex u[2][2], v[2][2], w[2][2];	/* temporary */
	double complex a[2][2], b[2][2];		/* result */

	/*
	 * Find the frequency.	For calibration measurements, space
	 * the frequencies linearly between C_FMIN and C_FMAX; for DUT
	 * measurements, take them from the known values.
	 */
	if (measurement == MEASURE_DUT) {
	    f = vnadata_get_frequency(dutp->dut_vdp,
		    dutp->dut_offset + findex);
	} else {
	    f = C_FMIN + (C_FMAX - C_FMIN) * (double)findex /
		(frequencies - 1);
	}
	s = 2.0 * PI * I * f;

	/*
	 * For all measurements, fill in the "a" matrix to simulate
	 * leakage in the VNA switch.  Send 2/3 of the signal to the
	 * intended port and 1/3 to the other.
	 */
	A11 = 2.0 / 3.0;
	A12 = 1.0 / 3.0;
	A21 = 1.0 / 3.0;
	A22 = 2.0 / 3.0;

	/*
	 * Find ABCD parameters for the errors at VNA port 1.
	 *   Detector is on the left; DUT is on the right.
	 */
	/* series inductor L1 */
	U11 = 1.0;
	U12 = L1 * s;
	U21 = 0.0;
	U22 = 1.0;

	/* series resistor R1 */
	V11 = 1.0;
	V12 = R1;
	V21 = 0.0;
	V22 = 1.0;
	multiply(u, v, w);

	/* shunt capacitor C1 */
	U11 = 1.0;
	U12 = 0.0;
	U21 = C1 * s;
	U22 = 1.0;
	multiply(w, u, port1_abcd);

	/*
	 * Find ABCD parameters for the errors at VNA port 2.
	 *   DUT is on the left; detector is on the right.
	 */
	/* shunt resistor R2 */
	U11 = 1.0;
	U12 = 0.0;
	U21 = 1.0 / R2;
	U22 = 1.0;

	/* series inductor L2 */
	V11 = 1.0;
	V12 = L2 * s;
	V21 = 0.0;
	V22 = 1.0;
	multiply(u, v, w);

	/* series capacitor C2 */
	U11 = 1.0;
	U12 = 1.0 / (C2 * s);
	U21 = 0.0;
	U22 = 1.0;
	multiply(w, u, port2_abcd);

	/*
	 * Calculate the b matrix for the requested measurement.
	 */
	switch (measurement) {
	case MEASURE_THROUGH:
	    /*
	     * Multiply the ABCD parameters of the two error boxes,
	     * convert to s-parameters and find b = s a.
	     */
	    multiply(port1_abcd, port2_abcd, u);
	    vnaconv_atos(u, u, z0);
	    multiply(u, a, b);
	    break;

	case MEASURE_REFLECT:
	    /*
	     * Calculate the measured reflect values.
	     */
	    {
		double complex zr = RR + RL * s;
		double complex gamma;

		vnaconv_ztosn(&zr, &gamma, z0, 1);
		vnaconv_atos(port1_abcd, u, z0);
		V11 = U11 + U12 * U21 * gamma / (1.0 - U22 * gamma);
		V12 = 0.0;
		V21 = 0.0;
		vnaconv_atos(port2_abcd, u, z0);
		V22 = U22 + U12 * U21 * gamma / (1.0 - U11 * gamma);
		multiply(v, a, b);
	    }
	    break;

	case MEASURE_LINE:
	    /*
	     * Multiply the ABCD parameters of the first error box,
	     * the line and the second error box. Next, convert to
	     * s-parameters and find b = s a.
	     */
	    {
		double complex gl;

		gl = LINE_LENGTH * ACTUAL_GAMMA(f);
		U11 = ccosh(gl);
		U12 = csinh(gl) * Z0;
		U21 = csinh(gl) / Z0;
		U22 = ccosh(gl);

		multiply(port1_abcd, u, v);
		multiply(v, port2_abcd, u);
		vnaconv_atos(u, u, z0);
		multiply(u, a, b);
	    }
	    break;

	case MEASURE_DUT:
	    /*
	     * Convert the actual s-parameters of the DUT to ABCD
	     * parameters.  Multiply the ABCD parameters of the first
	     * error box, the DUT and the second error box.  Finally,
	     * convert to s-parameters and find b = s a.
	     */
	    {
		int dut_index = dutp->dut_offset + findex;
		const double complex *dut_s;

		dut_s = vnadata_get_matrix(dutp->dut_vdp, dut_index);
		vnaconv_stoa((const row_type *)dut_s, u, z0);
		multiply(port1_abcd, u, v);
		multiply(v, port2_abcd, u);
		vnaconv_atos(u, u, z0);
		multiply(u, a, b);
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
    double complex s21_vector[C_FREQUENCIES];	/* vector of line s12=s21 */
    int line_guess;
    int unknown_reflect, unknown_line;
    int line_s[2][2];

    /*
     * Create the calibration container.
     */
    if ((vcp = vnacal_create(error_fn, /*error_arg*/NULL)) == NULL) {
	exit(2);
    }

    /*
     * Start a new calibration.
     */
    if ((vnp = vnacal_new_alloc(vcp, VNACAL_TE10, 2, 2,
		    C_FREQUENCIES)) == NULL) {
	exit(3);
    }

    /*
     * Print the length of the line standard.
     */
    (void)printf("# line length: %f mm\n", 1000.0 * LINE_LENGTH);

    /*
     * Make the calibration measurements for through, reflect and line
     * standards.  Normally, we would interact with the user between
     * each of these steps to get the user to connect each standard
     * in sequence.  In our simulated environment, we can skip this.
     * The frequency_vector is filled from the first measurement only;
     * the frequencies for the other calibration steps have to be the
     * same as the first.
     */

    /*
     * Add the through standard and set frequency vector.
     */
    vna_measure(/*dutp*/NULL, MEASURE_THROUGH, frequency_vector, a, b);
    vnacal_new_set_frequency_vector(vnp, frequency_vector);
    vnacal_new_add_through(vnp, &a[0][0], 2, 2, &b[0][0], 2, 2,
	    /*port1*/1, /*port2*/2);

    /*
     * Add the reflect standard.  We know it's symmetrical and
     * approximately a short, but we don't know it exactly.
     */
    unknown_reflect = vnacal_make_unknown_parameter(vcp, VNACAL_SHORT);
    vna_measure(/*dutp*/NULL, MEASURE_REFLECT, /*frequency_vector*/NULL, a, b);
    vnacal_new_add_double_reflect(vnp, &a[0][0], 2, 2, &b[0][0], 2, 2,
	    unknown_reflect, unknown_reflect, 1, 2);

    /*
     * Find the ideal S12 == S21 parameters of the line.  From them,
     * form a vector parameter we'll use as the initial guess of the
     * actual line parameter, make the unknown line parameter from the
     * initial guess, and delete the guess.
     */
    for (int findex = 0; findex < C_FREQUENCIES; ++findex) {
	double f = frequency_vector[findex];
	double complex gl = LINE_LENGTH * IDEAL_GAMMA(f);
	double complex abcd[2][2];
	double complex s[2][2];

	abcd[0][0] = ccosh(gl);
	abcd[0][1] = csinh(gl) * Z0;
	abcd[1][0] = csinh(gl) / Z0;
	abcd[1][1] = ccosh(gl);
	vnaconv_atos(abcd, s, z0);
	s21_vector[findex] = s[1][0];
    }
    if ((line_guess = vnacal_make_vector_parameter(vcp,
		    frequency_vector, C_FREQUENCIES, s21_vector)) == -1) {
	exit(4);
    }
    if ((unknown_line = vnacal_make_unknown_parameter(vcp, line_guess)) == -1) {
	exit(5);
    }

    /*
     * Add the line standard.  We know it's matched.  We know the length.
     * But we don't know the exact propagation constant.
     */
    line_s[0][0] = VNACAL_MATCH;
    line_s[0][1] = unknown_line;
    line_s[1][0] = unknown_line;
    line_s[1][1] = VNACAL_MATCH;
    vna_measure(/*dutp*/NULL, MEASURE_LINE, /*frequency_vector*/NULL, a, b);
    vnacal_new_add_line(vnp, &a[0][0], 2, 2, &b[0][0], 2, 2, &line_s[0][0],
	    /*port1*/1, /*port2*/2);

    /*
     * Solve for the error terms.
     */
    if (vnacal_new_solve(vnp) == -1) {
	exit(7);
    }

    /*
     * Print the initial guesses for the transmission and reflection
     * coefficients.
     */
    (void)printf("# initial guess values for T, R\n");
    for (int findex = 0; findex < C_FREQUENCIES; ++findex) {
	double f = C_FMIN + (C_FMAX - C_FMIN) * (double)findex /
	    (C_FREQUENCIES - 1);
	double complex t;

	t = vnacal_get_parameter_value(vcp, line_guess, f);
	(void)printf("%e %+e %+e %+e %+e\n",
		f, creal(t), cimag(t), -1.0, 0.0);
    }
    (void)printf("\n\n");

    /*
     * Print the solved transmission and reflection coefficients.
     */
    (void)printf("# solved T, R\n");
    for (int findex = 0; findex < C_FREQUENCIES; ++findex) {
	double f = C_FMIN + (C_FMAX - C_FMIN) * (double)findex /
	    (C_FREQUENCIES - 1);
	double complex t, r;

	t = vnacal_get_parameter_value(vcp, unknown_line, f);
	r = vnacal_get_parameter_value(vcp, unknown_reflect, f);
	(void)printf("%e %+e %+e %+e %+e\n",
		f, creal(t), cimag(t), creal(r), cimag(r));
    }
    (void)printf("\n\n");

    /*
     * Print the actual transmission and reflection coefficients.
     */
    (void)printf("# actual T, R\n");
    for (int findex = 0; findex < C_FREQUENCIES; ++findex) {
	double f = C_FMIN + (C_FMAX - C_FMIN) * (double)findex /
	    (C_FREQUENCIES - 1);
	double complex s = 2.0 * PI * I * f;
	double complex gl = LINE_LENGTH * ACTUAL_GAMMA(f);
	double complex zr = RR + RL * s;
	double complex u[2][2], v[2][2];
	double complex gamma;

	U11 = ccosh(gl);
	U12 = csinh(gl) * Z0;
	U21 = csinh(gl) / Z0;
	U22 = ccosh(gl);
	vnaconv_atos(u, v, z0);
	vnaconv_ztosn(&zr, &gamma, z0, 1);
	(void)printf("%e %+e %+e %+e %+e\n",
		f, creal(V21), cimag(V21), creal(gamma), cimag(gamma));
    }
    (void)printf("\n\n");

    /*
     * Add the new calibration to the vnacal_t structure and save.
     */
    if (vnacal_add_calibration(vcp, "cal-TE10", vnp) == -1) {
	exit(8);
    }
    if (vnacal_save(vcp, "TRL.vnacal") == -1) {
	exit(9);
    }
    vnacal_delete_parameter(vcp, unknown_line);
    vnacal_delete_parameter(vcp, line_guess);
    vnacal_delete_parameter(vcp, unknown_reflect);
    vnacal_new_free(vnp);
    vnacal_free(vcp);
    vcp = NULL;
}

/*
 * apply_calibration: apply the calibration to the DUT
 *     Normally, make_calibration and apply_calibration would be in
 *     separate programs, but to keep the example simple, we've only
 *     made them separate functions.
 */
static void apply_calibration()
{
    vnacal_t *vcp = NULL;
    vnadata_t *vdp_corrected = NULL;
    vnadata_t *vdp_actual = NULL;
    int frequencies;
    int offset;
    double frequency_vector[M_FREQUENCIES];
    double complex a11_vector[M_FREQUENCIES];
    double complex a12_vector[M_FREQUENCIES];
    double complex a21_vector[M_FREQUENCIES];
    double complex a22_vector[M_FREQUENCIES];
    double complex b11_vector[M_FREQUENCIES];
    double complex b12_vector[M_FREQUENCIES];
    double complex b21_vector[M_FREQUENCIES];
    double complex b22_vector[M_FREQUENCIES];
    double complex *a[2][2] = {
	{ a11_vector, a12_vector },
	{ a21_vector, a22_vector }
    };
    double complex *b[2][2] = {
	{ b11_vector, b12_vector },
	{ b21_vector, b22_vector }
    };
    dut_info_type dut_info;

    /*
     * Load the calibration file.
     */
    if ((vcp = vnacal_load("TRL.vnacal", error_fn,
		    /*error_arg*/NULL)) == NULL) {
	exit(10);
    }

    /*
     * Set up the simulated VNA.
     */
    dut_setup(&dut_info);
    vdp_actual = dut_info.dut_vdp;
    frequencies = dut_info.dut_frequencies;
    offset = dut_info.dut_offset;

    /*
     * Measure the DUT with errors using the simulated VNA.
     */
    vna_measure(&dut_info, MEASURE_DUT, frequency_vector, a, b);

    /*
     * Allocate a vnadata_t structure to receive the corrected S parameters.
     */
    if ((vdp_corrected = vnadata_alloc(error_fn, /*error_arg*/NULL)) == NULL) {
	exit(11);
    }

    /*
     * Apply the correction.
     */
    if (vnacal_apply(vcp, /*index*/0, frequency_vector, frequencies,
		&a[0][0], 2, 2, &b[0][0], 2, 2, vdp_corrected) == -1) {
	exit(12);
    }

    /*
     * Print the actual S-parameters from the device under test.
     */
    (void)printf("# actual\n");
    for (int findex = 0; findex < frequencies; ++findex) {
	int dindex = offset + findex;
	double f = vnadata_get_frequency(vdp_actual, dindex);
	double complex s11 = vnadata_get_cell(vdp_actual, dindex, 0, 0);
	double complex s12 = vnadata_get_cell(vdp_actual, dindex, 0, 1);
	double complex s21 = vnadata_get_cell(vdp_actual, dindex, 1, 0);
	double complex s22 = vnadata_get_cell(vdp_actual, dindex, 1, 1);

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
    vnadata_free(dut_info.dut_vdp);
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
