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

#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#include <vnacal.h>
#include <vnadata.h>

/*
 * Touchstone parameters of a 3.5mm male-to-male RF adapter
 */
#define MALE_TO_MALE_ADAPTER	"MtoM-3.5mm-adapter.s2p"

/*
 * Frequency range and reference impedance.
 */
#define F_MIN	100.0e+3
#define F_MAX	8.5e+9
#define F_STEPS	50
#define Z0	50.0

/*
 * VNALab EF35-101 Female 3.5mm Short Standard
 */
static const vnacal_calkit_data_t short_data = {
    .vcd_type = VNACAL_CALKIT_SHORT,
    .vcd_offset_delay = 33.340790e-12,
    .vcd_offset_loss = 5.460953e+9,
    .vcd_offset_z0 = 51.259595,
    .vcd_fmin = 0.0,
    .vcd_fmax = 8.5e+9,
    .vcd_l_coefficients = {
	-119.006943e-12,
	-1.310397249e-21,
	1.511773982-30,
	-91.480400e-42
    }
};

/*
 * error_fn: error printing function for the library
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
 * embed_example: apply male-to-male adapter to 3.5mm short standard
 *    Example of applying a male-to-make coaxial adapter to a female
 *    calibration standard to make it match a VNA expecting a male device.
 */
static void embed_example()
{
    vnacal_t *vcp = NULL;
    int short_parameter = -1;
    int fixture_parameter_matrix[2][2];
    int embedded_short_parameter = -1;

    /*
     * Create the calibration container structure.
     */
    if ((vcp = vnacal_create(error_fn, /*error_arg=*/NULL)) == NULL) {
	exit(1);
    }

    /*
     * Create a calibration parameter for the calkit short standard.
     */
    short_parameter = vnacal_make_calkit_parameter(vcp, &short_data);
    if (short_parameter == -1) {
	exit(2);
    }

    /*
     * Load the network parameter data for the RF adapter and
     * build a parameter matrix.
     */
    if (vnacal_load_data_parameter_matrix(vcp, MALE_TO_MALE_ADAPTER,
		&fixture_parameter_matrix[0][0],
		sizeof(fixture_parameter_matrix)) == -1) {
	exit(3);
    }

    /*
     * Embed the short standard in the fixture.
     */
    if ((embedded_short_parameter = vnacal_embed_parameter(vcp,
		    short_parameter, fixture_parameter_matrix)) == -1) {
	exit(4);
    }

    /*
     * Now, we can use embedded_short_parameter for calibration by
     * passing it to vnacal_new_add_single_reflect().  To keep the example
     * simple, we'll just print the S-parameters of the standard before
     * and after embedding.
     */
    (void)printf("# embed:\n");
    for (int findex = 0; findex < F_STEPS; ++findex) {
	double f = F_MIN + (F_MAX - F_MIN) / (double)(F_STEPS - 1) * findex;
	double complex original_value, embedded_value;

	if ((original_value = vnacal_eval_parameter(vcp,
			short_parameter, f, Z0)) == HUGE_VAL) {
	    exit(5);
	}
	if ((embedded_value = vnacal_eval_parameter(vcp,
			embedded_short_parameter, f, Z0)) == HUGE_VAL) {
	    exit(6);
	}
	(void)printf("%e %+f %+f %+f %+f\n", f,
		creal(original_value), cimag(original_value),
		creal(embedded_value), cimag(embedded_value));
    }
    (void)printf("\n\n");

    /*
     * Clean up.
     */
    vnacal_free(vcp);	/* implicitly deletes the parameters */
}

/*
 * deembed_example: remove adapter from measured data
 *    Example of recovering network parameter data from a device
 *    that was measured through the RF adapter above.
 */
static void deembed_example()
{
    int frequencies;
    vnacal_t *vcp = NULL;
    vnadata_t *vdp_measured = NULL;
    vnadata_t *vdp_dut = NULL;
    int fixture_parameter_matrix[2][2];

    /*
     * Create the calibration container structure.
     */
    if ((vcp = vnacal_create(error_fn, /*error_arg=*/NULL)) == NULL) {
	exit(10);
    }

    /*
     * Load the embedded DUT.  We could have obtained this data from
     * vnacal_apply() rather than reading from a file.
     */
    if ((vdp_measured = vnadata_alloc(error_fn, /*error_arg=*/NULL)) == NULL) {
	exit(11);
    }
    if (vnadata_load(vdp_measured, "embedded-DUT.npd") == -1) {
	exit(12);
    }
    frequencies = vnadata_get_frequencies(vdp_measured);

    /*
     * Build a parameter matrix for the male-to-make RF adapter.
     */
    if (vnacal_load_data_parameter_matrix(vcp, MALE_TO_MALE_ADAPTER,
		&fixture_parameter_matrix[0][0],
		sizeof(fixture_parameter_matrix)) == -1) {
	exit(13);
    }

    /*
     * Allocate a vnadata_t struct to hold the result.  Set
     * dimensions, reference impedances and frequency vector.
     * If we didn't still need the original data, we could skip
     * this step and do an in-place de-embedding on vdp_measured.
     */
    if ((vdp_dut = vnadata_alloc_and_init(error_fn, /*error_arg=*/NULL,
		    VPT_S, 1, 1, frequencies)) == NULL) {
	exit(14);
    }
    if (vnadata_set_z0_vector(vdp_dut,
		vnadata_get_z0_vector(vdp_measured)) == -1) {
	exit(15);
    }
    if (vnadata_set_frequency_vector(vdp_dut,
		vnadata_get_frequency_vector(vdp_measured)) == -1) {
	exit(16);
    }

    /*
     * De-embed the DUT from the adapter.
     */
    if (vnacal_deembed(vcp, vdp_measured, vdp_dut,
		&fixture_parameter_matrix[0][0], 2) == -1) {
	exit(17);
    }

    /*
     * Print the DUT S-parameters before and after de-embedding.
     */
    (void)printf("# de-embed:\n");
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = vnadata_get_frequency(vdp_measured, findex);
	double complex measured_value, deembedded_value;

	if ((measured_value = vnadata_get_cell(vdp_measured,
			findex, 0, 0)) == HUGE_VAL) {
	    exit(18);
	}
	if ((deembedded_value = vnadata_get_cell(vdp_dut,
			findex, 0, 0)) == HUGE_VAL) {
	    exit(19);
	}
	(void)printf("%e %+f %+f %+f %+f\n", f,
		creal(measured_value), cimag(measured_value),
		creal(deembedded_value), cimag(deembedded_value));
    }
    (void)printf("\n\n");

    /*
     * Clean up.
     */
    vnadata_free(vdp_dut);
    vnadata_free(vdp_measured);
    vnacal_free(vcp);
}

/*
 * main
 */
int main(int argc, char **argv)
{
    embed_example();
    deembed_example();

    exit(0);
}
