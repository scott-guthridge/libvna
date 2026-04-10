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
 * eval_embed: evaluate standard embedded in a fixture
 *   @stdp: vnacal_standard_t structure
 *   @function: name of user-called function
 *   @z0_vector: reference impedance vector
 *   @frequency: frequency (in Hz) to evaluate
 *   @result_matrix: caller-allocated matrix to hold result
 */
static int eval_embed(struct vnacal_standard *stdp, const char *function,
	const double complex *z0_vector, double frequency,
	double complex *result_matrix)
{
    vnacal_embed_standard_t *estdp = (vnacal_embed_standard_t *)stdp;
    const int target_ports = stdp->std_ports;
    const int target_cells = target_ports * target_ports;
    const int fixture_ports = 2 * target_ports;
    const int fixture_half_ports = fixture_ports / 2;
    const int fixture_cells = fixture_ports * fixture_ports;
    double complex target[target_cells];
    double complex fixture_z0_vector[fixture_ports];
    double complex fixture[fixture_cells];
    double complex el[target_ports][target_ports];
    double complex er[target_ports][target_ports];
    double complex et[target_ports][target_ports];
    double complex em[target_ports][target_ports];
    double complex u[target_ports][target_ports];
    double complex v[target_ports][target_ports];
    double complex w[target_ports][target_ports];

    /*
     * Evaluate the standard.
     */
    if (_vnacal_eval_parameter_matrix_i(function, estdp->estd_target_map,
		frequency, z0_vector, target) == -1) {
	return -1;
    }

    /*
     * Form the z0 vector for the fixture with standard-facing target_ports
     * in conjugate match with the standard, and evaluate the fixture.
     */
    assert(stdp->std_ops->stdo_type == VNACAL_EMBED);
    for (int i = 0; i < fixture_ports / 2; ++i) {
	double complex z0 = z0_vector[i];

	fixture_z0_vector[i] = z0;
	fixture_z0_vector[i + fixture_ports / 2] = conj(z0);
    }
    if (_vnacal_eval_parameter_matrix_i(function, estdp->estd_fixture_map,
		frequency, fixture_z0_vector, fixture) == -1) {
	return -1;
    }

    /*
     * Form the fixture el, er, et and em matrices.
     */
#define FIXTURE(r, c)	(fixture[fixture_ports * (r) + (c)])
    for (int i = 0; i < fixture_half_ports; ++i) {
	for (int j = 0; j < fixture_half_ports; ++j) {
	    el[i][j] = FIXTURE(i, j);
	    er[i][j] = FIXTURE(i, j + fixture_half_ports);
	    et[i][j] = FIXTURE(i + fixture_half_ports, j);
	    em[i][j] = FIXTURE(i + fixture_half_ports,
			       j + fixture_half_ports);
	}
    }
#undef FIXTURE

    /*
     * Compute: er * target * (I - em * target)^-1 * et + el
     *     u = er * target
     *     v = em * target
     *     v = I - v
     *     w = u / v
     *     result_matrix = w * et
     *     result_matrix += el
     */
    _vnacommon_mmultiply(&u[0][0], &er[0][0], target,
	    target_ports, target_ports, target_ports);
    _vnacommon_mmultiply(&v[0][0], &em[0][0], target,
	    target_ports, target_ports, target_ports);
    for (int row = 0; row < target_ports; ++row) {
	for (int column = 0; column < target_ports; ++column) {
	    v[row][column] = ((row == column) ? 1.0 : 0.0) - v[row][column];
	}
    }
    _vnacommon_mrdivide(&w[0][0], &u[0][0], &v[0][0],
	    target_ports, target_ports);
    _vnacommon_mmultiply(result_matrix, &w[0][0], &et[0][0],
	    target_ports, target_ports, target_ports);
    for (int row = 0; row < target_ports; ++row) {
	for (int column = 0; column < target_ports; ++column) {
	    result_matrix[target_ports * row + column] += el[row][column];
	}
    }
    return 0;
}

/*
 * eval_deembed: evaluate standard de-embedded from a fixture
 *   @stdp: vnacal_standard_t structure
 *   @function: name of user-called function
 *   @z0_vector: reference impedance vector
 *   @frequency: frequency (in Hz) to evaluate
 *   @result_matrix: caller-allocated matrix to hold result
 */
static int eval_deembed(struct vnacal_standard *stdp, const char *function,
	const double complex *z0_vector, double frequency,
	double complex *result_matrix)
{
    vnacal_embed_standard_t *estdp = (vnacal_embed_standard_t *)stdp;
    const int target_ports = stdp->std_ports;
    const int target_cells = target_ports * target_ports;
    const int fixture_ports = 2 * target_ports;
    const int fixture_half_ports = fixture_ports / 2;
    const int fixture_cells = fixture_ports * fixture_ports;
    double complex target[target_cells];
    double complex fixture_z0_vector[fixture_ports];
    double complex fixture[fixture_cells];
    double complex el[target_ports][target_ports];
    double complex er[target_ports][target_ports];
    double complex et[target_ports][target_ports];
    double complex em[target_ports][target_ports];
    double complex u[target_ports][target_ports];
    double complex v[target_ports][target_ports];
    double complex w[target_ports][target_ports];

    /*
     * Evaluate the embedded standard.
     */
    if (_vnacal_eval_parameter_matrix_i(function, estdp->estd_target_map,
		frequency, z0_vector, target) == -1) {
	return -1;
    }


    /*
     * Form the z0 vector for the fixture with standard-facing target_ports
     * in conjugate match with the standard, and evaluate the fixture.
     */
    assert(stdp->std_ops->stdo_type == VNACAL_DEEMBED);
    for (int i = 0; i < fixture_ports / 2; ++i) {
	double complex z0 = z0_vector[i];

	fixture_z0_vector[i] = z0;
	fixture_z0_vector[i + fixture_ports / 2] = conj(z0);
    }
    if (_vnacal_eval_parameter_matrix_i(function, estdp->estd_fixture_map,
		frequency, fixture_z0_vector, fixture) == -1) {
	return -1;
    }

    /*
     * Form the fixture el, er, et and em matrices.
     */
#define FIXTURE(r, c)	(fixture[fixture_ports * (r) + (c)])
    for (int i = 0; i < fixture_half_ports; ++i) {
	for (int j = 0; j < fixture_half_ports; ++j) {
	    el[i][j] = FIXTURE(i, j);
	    er[i][j] = FIXTURE(i, j + fixture_half_ports);
	    et[i][j] = FIXTURE(i + fixture_half_ports, j);
	    em[i][j] = FIXTURE(i + fixture_half_ports,
			       j + fixture_half_ports);
	}
    }
#undef FIXTURE

    /*
     * Compute: (er + v em) \ v
     *          where v = (target - el) / et
     *
     *   u = target - el
     *   v = u / et
     *   w = v * em
     *   w += er
     *   result_matrix = w \ v
     */
    for (int row = 0; row < target_ports; ++row) {
	for (int column = 0; column < target_ports; ++column) {
	    u[row][column] = target[target_ports * row + column] -
			     el[row][column];
	}
    }
    _vnacommon_mrdivide(&v[0][0], &u[0][0], &et[0][0],
	    target_ports, target_ports);
    _vnacommon_mmultiply(&w[0][0], &v[0][0], &em[0][0],
	    target_ports, target_ports, target_ports);
    for (int row = 0; row < target_ports; ++row) {
	for (int column = 0; column < target_ports; ++column) {
	    w[row][column] += er[row][column];
	}
    }
    _vnacommon_mldivide(result_matrix, &w[0][0], &v[0][0],
	    target_ports, target_ports);
    return 0;
}

/*
 * free_embed_standard: destruct the derived portion of vnacal_embed_standard_t
 *   @stdp: vnacal_standard_t structure
 */
static void free_embed_standard(vnacal_standard_t *stdp)
{
    vnacal_embed_standard_t *estdp;
    vnacal_parameter_type_t type = stdp->std_ops->stdo_type;
    const int ports = stdp->std_ports;
    const int target_cells = ports * ports;
    const int fixture_cells = 4 * target_cells;

    assert(type == VNACAL_EMBED || type == VNACAL_DEEMBED);
    estdp = (vnacal_embed_standard_t *)stdp;
    if (estdp != NULL) {
        _vnacal_free_parameter_matrix_map(estdp->estd_target_map);
	estdp->estd_target_map = NULL;
	if (estdp->estd_target_matrix != NULL) {
	    for (int cell = 0; cell < target_cells; ++cell) {
		vnacal_parameter_t **vpmrpp = &estdp->estd_target_matrix[cell];

		if (*vpmrpp != NULL) {
		    _vnacal_release_parameter(*vpmrpp);
		    *vpmrpp = NULL;
		}
	    }
	}
        _vnacal_free_parameter_matrix_map(estdp->estd_fixture_map);
	estdp->estd_fixture_map = NULL;
	if (estdp->estd_fixture_matrix != NULL) {
	    for (int cell = 0; cell < fixture_cells; ++cell) {
		vnacal_parameter_t **vpmrpp = &estdp->estd_fixture_matrix[cell];

		if (*vpmrpp != NULL) {
		    _vnacal_release_parameter(*vpmrpp);
		    *vpmrpp = NULL;
		}
	    }
	}
    }
}

/*
 * _vnacal_embed_ops: subclass operations for embed
 */
const vnacal_standard_ops_t _vnacal_embed_ops = {
    .stdo_type = VNACAL_EMBED,
    .stdo_eval = eval_embed,
    .stdo_free = free_embed_standard
};

/*
 * _vnacal_deembed_ops: subclass operations for deembed
 */
const vnacal_standard_ops_t _vnacal_deembed_ops = {
    .stdo_type = VNACAL_DEEMBED,
    .stdo_eval = eval_deembed,
    .stdo_free = free_embed_standard
};

/*
 * _vnacal_embed_parameter_matrix: embed/de-embed a DUT into/from a fixture
 *   @function: called function name
 *   @stdop: embed or de-embed subclass operations
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @target_matrix: parameter matrix of DUT/standard
 *   @target_ports: number of DUT/standard ports
 *   @fixture_matrix: parameter matrix of fixture
 *   @result_matrix: caller-allocated matrix to receive result
 *   @result_matrix_size: allocation in bytes of the result matrix
 *
 * Takes target_matrix describing a device under test or calibration
 * standard, embeds it using fixture_matrix, and fills result_matrix with
 * a new parameter matrix representing the embedded device, with the same
 * dimensions as target_matrix.  The fixture_matrix has dimensions
 * (2 * target_ports) x (2 * target_ports).  The first half of the fixture
 * ports face the VNA; the second half face the DUT.  Caller can delete
 * the returned parameters by a call to vnacal_delete_parameter_matrix.
 */
int _vnacal_embed_parameter_matrix(const char *function,
	const vnacal_standard_ops_t *stdop, vnacal_t *vcp,
	const int *target_matrix, int target_ports,
	const int *fixture_matrix,
	int *result_matrix, size_t result_matrix_size)
{

    vnacal_parameter_type_t type = stdop->stdo_type;
    vnacal_standard_t *stdp = NULL;
    vnacal_embed_standard_t *estdp = NULL;
    const int fixture_ports = 2 * target_ports;
    const int target_cells = target_ports * target_ports;
    const int fixture_cells = fixture_ports * fixture_ports;
    double fmin = 0.0;
    double fmax = INFINITY;

    assert(type == VNACAL_EMBED || type == VNACAL_DEEMBED);
    if (vcp == NULL || vcp->vc_magic != VC_MAGIC) {
	errno = EINVAL;
	return -1;
    }
    if (target_matrix == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: fixture_matrix cannot be NULL", function);
	return -1;
    }
    if (fixture_matrix == NULL) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: fixture_matrix cannot be NULL", function);
	return -1;
    }
    if (target_ports < 0) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: target_ports and fixture_ports must be non-negative",
		function);
	return -1;
    }
    if (target_ports * target_ports * sizeof(int) > result_matrix_size) {
	_vnacal_error(vcp, VNAERR_USAGE,
		"%s: insufficient allocation for result_matrix", function);
	return -1;
    }
    _vnacal_init_parameter_matrix(result_matrix, target_ports, target_ports);

    /*
     * Allocate and init the vnacal_standard_t structure.
     */
    if ((estdp = _vnacal_alloc_standard(function, vcp, stdop,
		    target_ports, sizeof(vnacal_embed_standard_t))) == NULL) {
	goto error;
    }
    stdp = &estdp->estd_base;
    if ((estdp->estd_target_matrix = calloc(target_cells,
		    sizeof(vnacal_parameter_t *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto error;
    }
    if ((estdp->estd_fixture_matrix = calloc(fixture_cells,
		    sizeof(vnacal_parameter_t *))) == NULL) {
	_vnacal_error(vcp, VNAERR_SYSTEM, "calloc: %s", strerror(errno));
	goto error;
    }
    for (int cell = 0; cell < target_cells; ++cell) {
	vnacal_parameter_t *vpmrp;
	double temp_fmin, temp_fmax;

	vpmrp = _vnacal_get_parameter(vcp, target_matrix[cell]);
	if (vpmrp == NULL) {
	    goto error;
	}
	_vnacal_get_parameter_frange(vpmrp, &temp_fmin, &temp_fmax);
	if (temp_fmin > fmin) {
	    fmin = temp_fmin;
	}
	if (temp_fmax < fmax) {
	    fmax = temp_fmax;
	}
	_vnacal_hold_parameter(vpmrp);
	estdp->estd_target_matrix[cell] = vpmrp;
    }
    for (int cell = 0; cell < fixture_cells; ++cell) {
	vnacal_parameter_t *vpmrp;
	double temp_fmin, temp_fmax;

	vpmrp = _vnacal_get_parameter(vcp, fixture_matrix[cell]);
	if (vpmrp == NULL) {
	    goto error;
	}
	_vnacal_get_parameter_frange(vpmrp, &temp_fmin, &temp_fmax);
	if (temp_fmin > fmin) {
	    fmin = temp_fmin;
	}
	if (temp_fmax < fmax) {
	    fmax = temp_fmax;
	}
	_vnacal_hold_parameter(vpmrp);
	estdp->estd_fixture_matrix[cell] = vpmrp;
    }
    estdp->estd_fmin = fmin;
    estdp->estd_fmax = fmax;
    if ((estdp->estd_target_map = _vnacal_analyze_parameter_matrix(__func__,
		    vcp, estdp->estd_target_matrix, target_ports, target_ports,
		    /*initial=*/false)) == NULL) {
	goto error;
    }
    if ((estdp->estd_fixture_map = _vnacal_analyze_parameter_matrix(__func__,
		    vcp, estdp->estd_fixture_matrix,
		    fixture_ports, fixture_ports,
		    /*initial=*/false)) == NULL) {
	goto error;
    }

    /*
     * Fill in parameter name.
     */
    {
	bool embed = stdp->std_ops->stdo_type == VNACAL_EMBED;
	char target_name[PARAMETER_BUFFER_ALLOC];
	char fixture_name[PARAMETER_BUFFER_ALLOC];

	_vnacal_get_parameter_name(estdp->estd_target_matrix[0],
		/*with_sxx=*/false, target_name);

	_vnacal_get_parameter_name(estdp->estd_fixture_matrix[0],
		/*with_sxx=*/false, fixture_name);
	if (asprintf(&stdp->std_name, "%s[%s, %s]",
		    embed ? "embed" : "deembed",
		    target_name, fixture_name) == -1) {
	    stdp->std_name = NULL;
	    _vnacal_error(vcp, VNAERR_SYSTEM, "asprintf: %s", strerror(errno));
	    goto error;
	}
    }

    /*
     * Fill result_matrix.
     */
    if (_vnacal_fill_standard_parameter_matrix(function, stdp,
		result_matrix) == -1) {
	goto error;
    }
    _vnacal_release_standard(&stdp);	/* release initial reference */
    assert(stdp != NULL);
    return target_ports;

error:
    if (stdp != NULL) {
	_vnacal_release_standard(&stdp); /* release initial reference */
	assert(stdp == NULL);
    }
    return -1;
}

/*
 * vnacal_embed_parameter: embed a single port DUT into a two-port fixture
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @dut: parameter describing the DUT or standard
 *   @fixture_matrix: 2x2 matrix describing the test fixture
 */
int vnacal_embed_parameter(vnacal_t *vcp, int dut,
	const int (*fixture_matrix)[2])
{
    int parameter;

    if (_vnacal_embed_parameter_matrix(__func__, &_vnacal_embed_ops,
		vcp, &dut, 1, *fixture_matrix,
		&parameter, sizeof(parameter)) == -1) {
	return -1;
    }
    return parameter;
}

/*
 * vnacal_embed_parameter_matrix: embed a DUT into a test fixture
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @dut_matrix: parameter matrix of DUT/standard
 *   @dut_ports: number of DUT/standard ports
 *   @fixture_matrix: parameter matrix of fixture
 *   @result_matrix: caller-allocated matrix to receive result
 *   @result_matrix_size: allocation in bytes of the result matrix
 *
 * Fill result_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_* functions.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
int vnacal_embed_parameter_matrix(vnacal_t *vcp,
	const int *dut_matrix, int dut_ports,
	const int *fixture_matrix,
	int *result_matrix, size_t result_matrix_size)
{
    return _vnacal_embed_parameter_matrix(__func__, &_vnacal_embed_ops,
	    vcp, dut_matrix, dut_ports, fixture_matrix,
	    result_matrix, result_matrix_size);
}

/*
 * vnacal_deembed_parameter: de-embed a single port DUT from two-port fixture
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @embedded: parameter describing the DUT or standard
 *   @fixture_matrix: 2x2 matrix describing the test fixture
 */
int vnacal_deembed_parameter(vnacal_t *vcp, int embedded,
	const int (*fixture_matrix)[2])
{
    int parameter;

    if (_vnacal_embed_parameter_matrix(__func__, &_vnacal_deembed_ops,
		vcp, &embedded, 1, *fixture_matrix,
		&parameter, sizeof(parameter)) == -1) {
	return -1;
    }
    return parameter;
}

/*
 * vnacal_deembed_parameter_matrix: de-embed a DUT from a test fixture
 *   @vcp: pointer returned from vnacal_create or vnacal_load
 *   @embedded_matrix: parameter matrix of DUT/standard
 *   @embedded_ports: number of DUT/standard ports
 *   @fixture_matrix: parameter matrix of fixture
 *   @result_matrix: caller-allocated matrix to receive result
 *   @result_matrix_size: allocation in bytes of the result matrix
 *
 * Fill result_matrix with parameter indices suitable for passing to
 * the vnacal_new_add_* functions.
 *
 * Returns the number of ports (rows and columns) of the standard.
 * Caller can delete the returned parameters by a call to
 * vnacal_delete_parameter_matrix.
 */
int vnacal_deembed_parameter_matrix(vnacal_t *vcp,
	const int *embedded_matrix, int embedded_ports,
	const int *fixture_matrix,
	int *result_matrix, size_t result_matrix_size)
{
    return _vnacal_embed_parameter_matrix(__func__, &_vnacal_deembed_ops,
	    vcp, embedded_matrix, embedded_ports, fixture_matrix, 
	    result_matrix, result_matrix_size);
}
