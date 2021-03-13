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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <vnacal.h>

/*
 * MAX: find maximum of two numbers
 */
#ifndef MAX
#define MAX(x, y)	((x) >= (y) ? (x) : (y))
#endif /* MAX */

/*
 * FMIN, FMAX: frequency range of the VNA in Hz
 */
#define FMIN		10e+3
#define FMAX		100e+6

/*
 * C_ROWS, C_COLUMNS, C_FREQUENCIES: calibration dimensions
 *   Calibration matrix is 2x1, i.e. the VNA drives signal and
 *   measures reflected power on the first port only, and measures
 *   forward power on the second port only.  C_FREQUENCIES is the
 *   number of frequency points used for the calibration.
 */
#define C_ROWS		  2
#define C_COLUMNS	  1
#define C_FREQUENCIES	 79

/*
 * M_ROWS, M_COLUMNS, C_FREQUENCIES: measurement dimensions
 *   We measure full 2x2 S-parameters from the device under test.
 *   The number of frequency points used in the measurement doesn't
 *   have to match the calibration -- the library interpolates
 *   between error parameters when necessary.
 */
#define M_ROWS		  2
#define M_COLUMNS	  2
#define M_FREQUENCIES	100

/*
 * PI, Z0, W1, W2: misc constants
 *   PI is used below to convert from Hz to angular frequency
 *   Z0 is the system impedance
 *   W1 is the undamped natural frequency of the errors in our VNA
 *   W2 is the undamped natural frequency of our simulated DUT
 */
#define PI	3.14159265358979323846264338327950288419716939937508
#define Z0	50.0
#define W1	(2 * PI * 10e+6)
#define W2	(2 * PI * 1e+6)

/*
 * measurement_t: which simulated measurement should vna_measure return
 */
typedef enum measurement {
    SHORT_CALIBRATION,
    OPEN_CALIBRATION,
    LOAD_CALIBRATION,
    THROUGH_CALIBRATION,
    FORWARD_MEASUREMENT,
    REVERSE_MEASUREMENT
} measurement_t;

/*
 * vna_measure: simulate the requested VNA measurement
 *   @measurement: which measurement to simulate
 *   @frequencies: number of frequency points
 *   @m_frequency_vector: returned vector of frequencies
 *   @detector1_vector: returned voltages from detector 1
 *   @detector2_vector: returned voltages from detector 2
 *
 *   Our simulated VNA has two flaws: first, there is a
 *   stray capacitance of 1 / (Z0 * W1) [318pF] between
 *   port 1 and ground; second, there is an inductance of
 *   Z0 / W1 [796 nH] in series with port 2.
 *
 *   The simulated device under test (DUT) is a second
 *   order LC divider low pass filter with L = Z0 / W2
 *   [7.96 μH] and C = 1 / (Z0 * W2) [3.18nF].
 *
 */
static int vna_measure(measurement_t measurement,
	int frequencies, double *m_frequency_vector,
	double complex *detector1_vector,
	double complex *detector2_vector)
{
    double c = log(FMAX / FMIN);

    /*
     * For each frequency FMIN to FMAX spaced uniformly on a log
     * scale...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = FMIN * exp((double)findex / (frequencies - 1) * c);
	double complex s = I * 2 * PI * f;
	double complex d, detector1, detector2;

	switch (measurement) {
	case SHORT_CALIBRATION:
	    /*
	     * The shorted calibration standard on port 1 shunts
	     * out the stray capacitance, giving a perfect gamma
	     * value of -1.  Port 1 is connected to a terminator
	     * and receives no signal, but the detector picks up
	     * a bit of internal noise.
	     */
	    detector1 = -1.0;
	    detector2 =  0.1;
	    break;

	case OPEN_CALIBRATION:
	    /*
	     * The open calibration standard exposes the stray
	     * capacitance on port 1.  Port 1 continues to pick up
	     * internal noise.
	     */
	    detector1 = (1.0 - s/W1) / (1.0 + s/W1);
	    detector2 = -0.3;
	    break;

	case LOAD_CALIBRATION:
	    /*
	     * The load calibration is in parallel with the stray
             * capacitance on port 1.  Port 1 picks up yet more
	     * internal noise.
	     */
	    detector1 = -s / (s + 2*W1);
	    detector2 = 0.2;
	    break;

	case THROUGH_CALIBRATION:
	    /*
	     * In the through configuration, the stray capacitance
	     * on port 1 and stray inductance on port 2 form a
	     * resonant circuit with a high-pass reflected signal
	     * and low-pass transmitted signal.
	     */
	    d = s*s + 2*W1*s + 2*W1*W1;
	    detector1 = -s*s / d;
	    detector2 = 2*W1*W1 / d;
	    break;

	case FORWARD_MEASUREMENT:
	    /*
	     * In the forward configuration, the DUT forms a fourth
	     * order resonant circuit with the stray impedances of
	     * the VNA.
	     */
	    d = s*s*s*s + 2*W1*s*s*s + (W1+W2)*(W1+W2)*s*s
		+ 2*W1*W2*(W1+W2)*s + 2*W1*W1*W2*W2;
	    detector1 = -(s*s*s*s - (W1*W1 - 2*W1*W2 - W2*W2)*s*s) / d;
	    detector2 = 2*W1*W1*W2*W2 / d;
	    break;

	case REVERSE_MEASUREMENT:
	    /*
	     * In the reverse configuration, the stray capacitance on
	     * port 1 is in parallel with the DUT capacitor and the
	     * stray inductance on port 2 is in series with the DUT
	     * inductor forming only a second order resonant circuit.
	     */
	    d = s*s + 2*W1*W2/(W1+W2)*s + 2*W1*W1*W2*W2/((W1+W2)*(W1+W2));
	    detector1 = -s*s / d;
	    detector2 = 2*W1*W1*W2*W2/((W1+W2)*(W1+W2)) / d;
	    break;

	default:
	    abort();
	}

	/*
	 * Return the requested vectors.
	 */
	if (m_frequency_vector != NULL)
	    m_frequency_vector[findex] = f;
	if (detector1_vector != NULL)
	    detector1_vector[findex] = detector1;
	if (detector2_vector != NULL)
	    detector2_vector[findex] = detector2;
    }
    return 0;
}

/*
 * error_fn: error printing function for the library
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
    vnacal_new_t *vnp;
    double frequency_vector[C_FREQUENCIES];
    double complex m_vector1[C_FREQUENCIES];
    double complex m_vector2[C_FREQUENCIES];
    double complex *const m[C_ROWS][C_COLUMNS] = {
	{ m_vector1 },
	{ m_vector2 }
    };
    vnacal_t *vcp;

    /*
     * Create the calibration container structure.
     */
    if ((vcp = vnacal_create(error_fn, /*error_arg=*/NULL)) == NULL) {
	exit(1);
    }

    /*
     * Create a new calibration.
     */
    if ((vnp = vnacal_new_alloc(vcp, VNACAL_E12,
		    C_ROWS, C_COLUMNS, C_FREQUENCIES)) == NULL) {
	exit(2);
    }

    /*
     * Make the calibration measurements for short, open, load and
     * through standards.  Normally, we would interact with the
     * user between each of these steps to get the user to connect
     * each standard in sequence.  In our simulated environment,
     * we can skip that part.  The frequency m_vector is filled
     * from the first measurement only -- the frequencies for the
     * other calibration steps have to be the same as the first.
     * The three leakage measurements are averaged.
     */

    /*
     * Short calibration
     */
    vna_measure(SHORT_CALIBRATION, C_FREQUENCIES,
	    frequency_vector, m_vector1, m_vector2);
    vnacal_new_set_frequency_vector(vnp, frequency_vector);
    vnacal_new_add_single_reflect_m(vnp, &m[0][0], C_ROWS, C_COLUMNS,
	    VNACAL_SHORT, 1);

    /*
     * Open calibration
     */
    vna_measure(OPEN_CALIBRATION, C_FREQUENCIES, NULL, m_vector1, m_vector2);
    vnacal_new_add_single_reflect_m(vnp, &m[0][0], C_ROWS, C_COLUMNS,
	    VNACAL_OPEN, 1);

    /*
     * Load calibration
     */
    vna_measure(LOAD_CALIBRATION, C_FREQUENCIES, NULL, m_vector1, m_vector2);
    vnacal_new_add_single_reflect_m(vnp, &m[0][0], C_ROWS, C_COLUMNS,
	    VNACAL_MATCH, 1);

    /*
     * Through calibration.
     */
    vna_measure(THROUGH_CALIBRATION, C_FREQUENCIES,
        NULL, m_vector1, m_vector2);
    vnacal_new_add_through_m(vnp, &m[0][0], C_ROWS, C_COLUMNS, 1, 2);

    /*
     * Solve for the error terms.
     */
    if (vnacal_new_solve(vnp) == -1) {
	exit(3);
    }

    /*
     * Add the calibration and save to a file.
     */
    if (vnacal_add_calibration(vcp, "cal_2x1", vnp) == -1) {
	exit(4);
    }
    if (vnacal_save(vcp, "SOLT.vnacal") == -1) {
	exit(5);
    }
    vnacal_new_free(vnp);
    vnacal_free(vcp);
    vcp = NULL;
}

/*
 * apply_calibration: apply the calibration to the simulated device
 *
 *   Normally, make_calibration and apply_calibration would be in
 *   separate programs, but to keep the example simple, we've just made
 *   them separate functions.
 */
static void apply_calibration()
{
    vnacal_t *vcp;
    double frequency_vector[M_FREQUENCIES];
    double complex m_vector11[M_FREQUENCIES];
    double complex m_vector12[M_FREQUENCIES];
    double complex m_vector21[M_FREQUENCIES];
    double complex m_vector22[M_FREQUENCIES];
    double complex *m[M_ROWS][M_COLUMNS] = {
	{ m_vector11, m_vector12 },
	{ m_vector21, m_vector22 }
    };
    vnadata_t *s_matrix;

    /*
     * Load the calibration file.
     */
    if ((vcp = vnacal_load("SOLT.vnacal", error_fn,
		    /*error_arg=*/NULL)) == NULL) {
	exit(6);
    }

    /*
     * Make the forward and reverse measurements of the device under
     * test.  We would normally have to interact with the user between
     * these steps in order to get the user to swap the connections.
     * Alternatively, if the VNA has a relay to swap ports automatically,
     * we would send different relay codes for these two measurements.
     * Note though, that if the VNA has a relay to swap ports, we'd
     * want to make a 2x2 calibration matrix above instead of 2x1 so
     * that the calibration also covers the relay.
     */

    /*
     * Make the forward measurement.
     */
    vna_measure(FORWARD_MEASUREMENT, M_FREQUENCIES, frequency_vector,
	    m_vector11, m_vector21);

    /*
     * Make the reverse measurement.
     */
    vna_measure(REVERSE_MEASUREMENT, M_FREQUENCIES, NULL,
	    m_vector22, m_vector12);

    /*
     * First, calculate and print the S-parameters we would expect
     * from the device under test if we measured them with a
     * perfect VNA.
     */
    (void)printf("# expected\n");
    for (int i = 0; i < M_FREQUENCIES; ++i) {
	double complex s = 2 * PI * I * frequency_vector[i];
	double complex d = s*s + 2*W2*s + 2*W2*W2;
	double complex s11 =  s*s     / d;
	double complex s12 =  2*W2*W2 / d;
	double complex s21 =  2*W2*W2 / d;
	double complex s22 = -s*s     / d;

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[i],
		creal(s11), cimag(s11),
		creal(s12), cimag(s12),
		creal(s21), cimag(s21),
		creal(s22), cimag(s22));
    }
    (void)printf("\n\n");

    /*
     * Next, print the values as measured from the imperfect VNA.
     */
    (void)printf("# measured\n");
    for (int i = 0; i < M_FREQUENCIES; ++i) {
	double complex m11 = m_vector11[i];
	double complex m12 = m_vector12[i];
	double complex m21 = m_vector21[i];
	double complex m22 = m_vector22[i];

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[i],
		creal(m11), cimag(m11),
		creal(m12), cimag(m12),
		creal(m21), cimag(m21),
		creal(m22), cimag(m22));
    }
    (void)printf("\n\n");

    /*
     * Allocate a vnadata_t structure to receive the computed S parameters.
     */
    if ((s_matrix = vnadata_alloc()) == NULL) {
	(void)fprintf(stderr, "example: vnadata_alloc: %s\n",
		strerror(errno));
	exit(7);
    }

    /*
     * Apply the calibration and report the corrected values.
     */
    if (vnacal_apply_m(vcp, /*index*/0, frequency_vector, M_FREQUENCIES,
		&m[0][0], M_ROWS, M_COLUMNS, s_matrix) == -1) {
	exit(8);
    }
    (void)printf("# corrected\n");
    for (int i = 0; i < M_FREQUENCIES; ++i) {
	double complex s11, s12, s21, s22;

	s11 = vnadata_get_cell(s_matrix, i, 0, 0);
	s12 = vnadata_get_cell(s_matrix, i, 0, 1);
	s21 = vnadata_get_cell(s_matrix, i, 1, 0);
	s22 = vnadata_get_cell(s_matrix, i, 1, 1);
	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[i],
		creal(s11), cimag(s11),
		creal(s12), cimag(s12),
		creal(s21), cimag(s21),
		creal(s22), cimag(s22));
    }
    vnadata_free(s_matrix);
    vnacal_free(vcp);
}

/*
 * main
 */
int main(int argc, char **argv)
{
    make_calibration();
    apply_calibration();

    exit(0);
}
