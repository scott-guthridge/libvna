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
 *   measures reflected power on the first port only.  It measures
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
 *   between error parameters if necessary.
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
 *   @detector0_vector: returned voltages from detector 0
 *   @detector1_vector: returned voltages from detector 1
 *
 *   To avoid confusion, we refer to the two ports of the VNA as
 *   port 0 and port 1 (as opposed to 1 and 2) to match C array
 *   indices, which start with zero.
 *
 *   Our simulated VNA has two flaws: first, there is a stray
 *   capacitance of 1 / (Z0 * W1) [318pF] between port 0 and ground;
 *   second, there is an inductance of Z0 / W1 [796 nH] in series
 *   with port 1.
 *
 *   The simulated device under test (DUT) is a second order
 *   LC divider low pass filter with L = Z0 / W2 [7.96 μH] and
 *   C = 1 / (Z0 * W2) [3.18nF].
 *
 */
static int vna_measure(measurement_t measurement,
	int frequencies, double *m_frequency_vector,
	double complex *detector0_vector,
	double complex *detector1_vector)
{
    double c = log(FMAX / FMIN);

    /*
     * For each frequency FMIN to FMAX spaced uniformly on a log
     * scale...
     */
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = FMIN * exp((double)findex / (frequencies - 1) * c);
	double complex s = I * 2 * PI * f;
	double complex d, detector0, detector1;

	switch (measurement) {
	case SHORT_CALIBRATION:
	    /*
	     * The shorted calibration standard on port 0 shunts
	     * out the stray capacitance, giving a perfect gamma
	     * value of -1.  Port 1 is connected to a terminator
	     * and receives no signal, but the detector picks up
	     * a bit of internal noise.
	     */
	    detector0 = -1.0;
	    detector1 =  0.1;
	    break;

	case OPEN_CALIBRATION:
	    /*
	     * The open calibration standard exposes the stray
	     * capacitance on port 0.  Port 1 continues to pick up
	     * internal noise.
	     */
	    detector0 = (1.0 - s/W1) / (1.0 + s/W1);
	    detector1 = -0.3;
	    break;

	case LOAD_CALIBRATION:
	    /*
	     * The load calibration is in parallel with the stray
             * capacitance on port 0.  Port 1 picks up yet more
	     * internal noise.
	     */
	    detector0 = -s / (s + 2*W1);
	    detector1 = 0.2;
	    break;

	case THROUGH_CALIBRATION:
	    /*
	     * In the through configuration, the stray capacitance
	     * on port 0 and stray inductance on port 1 form a
	     * resonant circuit with a high-pass reflected signal
	     * and low-pass transmitted signal.
	     */
	    d = s*s + 2*W1*s + 2*W1*W1;
	    detector0 = -s*s / d;
	    detector1 = 2*W1*W1 / d;
	    break;

	case FORWARD_MEASUREMENT:
	    /*
	     * In the forward configuration, the DUT forms a fourth
	     * order resonant circuit with the stray impedances of
	     * the VNA.
	     */
	    d = s*s*s*s + 2*W1*s*s*s + (W1+W2)*(W1+W2)*s*s
		+ 2*W1*W2*(W1+W2)*s + 2*W1*W1*W2*W2;
	    detector0 = -(s*s*s*s - (W1*W1 - 2*W1*W2 - W2*W2)*s*s) / d;
	    detector1 = 2*W1*W1*W2*W2 / d;
	    break;

	case REVERSE_MEASUREMENT:
	    /*
	     * In the reverse configuration, the stray capacitance on
	     * port 0 is in parallel with the DUT capacitor and the
	     * stray inductance on port 1 is in series with the DUT
	     * inductor forming only a second order resonant circuit.
	     */
	    d = s*s + 2*W1*W2/(W1+W2)*s + 2*W1*W1*W2*W2/((W1+W2)*(W1+W2));
	    detector0 = -s*s / d;
	    detector1 = 2*W1*W1*W2*W2/((W1+W2)*(W1+W2)) / d;
	    break;

	default:
	    abort();
	}

	/*
	 * Return the requested vectors.
	 */
	if (m_frequency_vector != NULL)
	    m_frequency_vector[findex] = f;
	if (detector0_vector != NULL)
	    detector0_vector[findex] = detector0;
	if (detector1_vector != NULL)
	    detector1_vector[findex] = detector1;
    }
    return 0;
}

/*
 * error_fn: error printing function for the library
 *   @message: single line error message without a newline
 *   @error_arg: passed through to the error function (unused here)
 */
static void error_fn(const char *message, void *error_arg)
{
    (void)fprintf(stderr, "example: %s\n", message);
}

/*
 * main
 */
int main(int argc, char **argv)
{
    vnacal_calset_t *vcsp;
    double frequency_vector[MAX(C_FREQUENCIES, M_FREQUENCIES)];
    double complex vector0[MAX(C_FREQUENCIES, M_FREQUENCIES)];
    double complex vector1[MAX(C_FREQUENCIES, M_FREQUENCIES)];
    vnacal_input_t *vip;
    vnadata_t *s_matrix;
    vnacal_t *vcp;

    /*
     * Allocate the structure to hold the calibration measurements.
     */
    if ((vcsp = vnacal_calset_alloc(/*setname=*/"default",
		    C_ROWS, C_COLUMNS, C_FREQUENCIES,
		    error_fn, /*error_arg=*/NULL)) == NULL) {
	(void)fprintf(stderr, "vnacal_calset_alloc: %s\n",
		strerror(errno));
	exit(2);
    }

    /*
     * Make the calibration measurements for short, open, load and
     * through standards.  Normally, we would interact with the
     * user between each of these steps to get the user to connect
     * each standard in sequence.  In our simulated environment,
     * we can skip that part.  The frequency vector is filled
     * from the first measurement only -- the frequencies for the
     * other calibration steps have to be the same as the first.
     * The three leakage measurements are averaged.
     */

    /*
     * Short calibration
     */
    vna_measure(SHORT_CALIBRATION, C_FREQUENCIES,
	    frequency_vector, vector0, vector1);
    vnacal_calset_set_frequency_vector(vcsp, frequency_vector);
    vnacal_calset_add_vector(vcsp, 0, 0, VNACAL_Sii_REF0, vector0);
    vnacal_calset_add_vector(vcsp, 1, 0, VNACAL_Sij_LEAKAGE, vector1);

    /*
     * Open calibration
     */
    vna_measure(OPEN_CALIBRATION, C_FREQUENCIES, NULL, vector0, vector1);
    vnacal_calset_add_vector(vcsp, 0, 0, VNACAL_Sii_REF1, vector0);
    vnacal_calset_add_vector(vcsp, 1, 0, VNACAL_Sij_LEAKAGE, vector1);

    /*
     * Load calibration
     */
    vna_measure(LOAD_CALIBRATION, C_FREQUENCIES, NULL, vector0, vector1);
    vnacal_calset_add_vector(vcsp, 0, 0, VNACAL_Sii_REF2, vector0);
    vnacal_calset_add_vector(vcsp, 1, 0, VNACAL_Sij_LEAKAGE, vector1);

    /*
     * Through calibration.
     */
    vna_measure(THROUGH_CALIBRATION, C_FREQUENCIES, NULL, vector0, vector1);
    vnacal_calset_add_vector(vcsp, 1, 0, VNACAL_Sjj_THROUGH, vector0);
    vnacal_calset_add_vector(vcsp, 1, 0, VNACAL_Sij_THROUGH, vector1);

    /*
     * Create the calibration from the measurements and save it to
     * a file.
     */
    if ((vcp = vnacal_create(/*sets=*/1, &vcsp, error_fn,
		    /*error_arg=*/NULL)) == NULL) {
	(void)fprintf(stderr, "vnacal_create: %s\n",
		strerror(errno));
	exit(3);
    }
    if (vnacal_save(vcp, "example.vnacal", ".excal") == -1) {
	(void)fprintf(stderr, "vnacal_save: %s\n",
		strerror(errno));
	exit(4);
    }
    vnacal_calset_free(vcsp);
    vnacal_free(vcp);
    vcp = NULL;

    /*
     * Now, use the calibration we made above to correct imperfect
     * measurements of the device under test.  Starting here, we
     * would normally be in a different program, but to keep the
     * example shorter we've combined them.
     *
     * Begin by loading the saved calibration.
     */
    if ((vcp = vnacal_load("example.vnacal", ".excal", error_fn,
		    /*error_arg=*/NULL)) == NULL) {
	(void)fprintf(stderr, "vnacal_load: %s\n",
		strerror(errno));
	exit(5);
    }

    /*
     * Allocate a vnacal_input_t object to apply the calibration
     * to measured values.
     */
    if ((vip = vnacal_input_alloc(vcp, /*set=*/0,
		    M_ROWS, M_COLUMNS, M_FREQUENCIES)) == NULL) {
	(void)fprintf(stderr, "example: vnadata_alloc_and_init: %s\n",
		strerror(errno));
	exit(6);
    }

    /*
     * Allocate a vnadata_t object to hold the S parameters.
     */
    if ((s_matrix = vnadata_alloc()) == NULL) {
	(void)fprintf(stderr, "example: vnadata_alloc: %s\n",
		strerror(errno));
	exit(7);
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
     * Forward measurement
     */
    vna_measure(FORWARD_MEASUREMENT, M_FREQUENCIES, frequency_vector,
	    vector0, vector1);
    vnacal_input_set_frequency_vector(vip, frequency_vector);
    vnacal_input_add_vector(vip, 0, 0, vector0);
    vnacal_input_add_vector(vip, 1, 0, vector1);

    /*
     * Reverse measurement
     */
    vna_measure(REVERSE_MEASUREMENT, M_FREQUENCIES, NULL,
	    vector0, vector1);
    vnacal_input_add_vector(vip, 1, 1, vector0);
    vnacal_input_add_vector(vip, 0, 1, vector1);

    /*
     * First, calculate and print the S-parameters we would expect
     * from the device under test if we measured them with a
     * perfect VNA.
     */
    (void)printf("# expected\n");
    for (int i = 0; i < M_FREQUENCIES; ++i) {
	double complex s = 2 * PI * I * frequency_vector[i];
	double complex d = s*s + 2*W2*s + 2*W2*W2;
	double complex s00 =  s*s     / d;
	double complex s01 =  2*W2*W2 / d;
	double complex s10 =  2*W2*W2 / d;
	double complex s11 = -s*s     / d;

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[i],
		creal(s00), cimag(s00),
		creal(s01), cimag(s01),
		creal(s10), cimag(s10),
		creal(s11), cimag(s11));
    }
    (void)printf("\n\n");

    /*
     * Now print the values as measured from the imperfect VNA.
     */
    (void)printf("# measured\n");
    for (int i = 0; i < M_FREQUENCIES; ++i) {
	double complex m00 = vnacal_input_get_value(vip, 0, 0, i);
	double complex m01 = vnacal_input_get_value(vip, 0, 1, i);
	double complex m10 = vnacal_input_get_value(vip, 1, 0, i);
	double complex m11 = vnacal_input_get_value(vip, 1, 1, i);

	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[i],
		creal(m00), cimag(m00),
		creal(m01), cimag(m01),
		creal(m10), cimag(m10),
		creal(m11), cimag(m11));
    }
    (void)printf("\n\n");

    /*
     * Apply the calibration to the measured data and print the
     * corrected s_matrix values.
     */
    if (vnacal_input_apply(vip, s_matrix) == -1) {
	(void)fprintf(stderr, "vnacal_input_apply: %s\n",
		strerror(errno));
	exit(8);
    }
    (void)printf("# calibrated\n");
    for (int i = 0; i < M_FREQUENCIES; ++i) {
	double complex s00, s01, s10, s11;

	s00 = vnadata_get_cell(s_matrix, i, 0, 0);
	s01 = vnadata_get_cell(s_matrix, i, 0, 1);
	s10 = vnadata_get_cell(s_matrix, i, 1, 0);
	s11 = vnadata_get_cell(s_matrix, i, 1, 1);
	(void)printf("%e %+e %+e %+e %+e %+e %+e %+e %+e\n",
		frequency_vector[i],
		creal(s00), cimag(s00),
		creal(s01), cimag(s01),
		creal(s10), cimag(s10),
		creal(s11), cimag(s11));
    }
    vnacal_input_free(vip);
    vnacal_free(vcp);

    exit(0);
    /*NOTREACHED*/
}
