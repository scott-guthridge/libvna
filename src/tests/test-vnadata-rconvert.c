/*
 * Vector Network Analyzer Library
 * Copyright © 2020-2024 D Scott Guthridge <scott_guthridge@rompromity.net>
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
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "vnacommon_internal.h"
#include "vnadata_internal.h"
#include "libt.h"
#include "libt_crand.h"

#define FREQUENICES	10
#define TRIALS		50

/*
 * TEST_EQUAL: fail the test if x and y are not equal
 *   Assumes local variable "result" and label "out" are defined.
 */
#define TEST_EQUAL(x, y, label) \
    if (opt_a) { \
	assert(libt_isequal_label((x), (y), (label))); \
    } else { \
	if (!libt_isequal_label((x), (y), (label))) { \
	    result = T_FAIL; \
	    goto out; \
	} \
    }

/*
 * CACP_ALIGN: find the extra alignment (if any) of a double
 * 	complex following a double complex pointer
 */
struct cacp_ {
    double complex *p;
    double complex v;
};
#define CACP_ALIGN	\
    (sizeof(struct cacp_) - \
     sizeof(double complex *) - \
     sizeof(double complex))

/*
 * common_t: common parameters
 */
typedef struct common {
    int frequencies;		/* number of frequencies */
    int ports;			/* number of ports */
    bool use_fz01;		/* use frequency-dependent z01 */
    bool use_fz02;		/* use frequency-dependent z02 */
    double *frequency_vector;	/* vector of frequencies */
    double complex **z01;	/* reference impedances 1 */
    double complex **z02;	/* reference impedances 2 */
    double complex **a1;	/* incident root power in reference 1 */
    double complex **b1;	/* reflected root power in reference 1 */
    double complex **v;		/* voltage matrix */
    double complex **i;		/* current matrix */
    double complex **a2;	/* incident root power in reference 2 */
    double complex **b2;	/* reflected root power in reference 2 */
    vnadata_t *vdp1;		/* data before conversion */
    vnadata_t *vdp2;		/* data after conversion */
} common_t;

/*
 * Options
 */
char *progname;
static const char options[] = "av";
static const char *const usage[] = {
    "[-av]",
    NULL
};
static const char *const help[] = {
    "-a	 abort on data miscompare",
    "-v	 show verbose output",
    NULL
};
bool opt_a = false;
int opt_v = 0;

/*
 * check_cmatrix: check that two complexes matrices are equal
 *   @lable: test context
 *   @a: serialized first matrix
 *   @b: serialized second matrix
 *   @ports: dimension (assuming square)
 */
static libt_result_t check_cmatrix(const char *label,
	const double complex *a, const double complex *b, int ports)
{
    libt_result_t result = T_SKIPPED;

    for (int row = 0; row < ports; ++row) {
	for (int column = 0; column < ports; ++column) {
	    int cell = ports * row + column;

	    TEST_EQUAL(a[cell], b[cell], label);
	}
    }
    result = T_PASS;

out:
    return result;
}

/*
 * alloc_matrix_vector: allocate a vector of pointers to complex matrix
 *   @frequencies: number of frequencies
 *   @elements: number of elements in each matrix
 *
 * The two objects are placed in a single allocation; thus the data
 * structure can be freed by a single call to free.
 */
static double complex **alloc_matrix_vector(size_t frequencies,
					    size_t elements)
{
    const size_t size1 = frequencies * sizeof(double complex *) + CACP_ALIGN;
    const size_t size2 = frequencies * elements * sizeof(double complex);
    double complex **result;
    double complex *matrix;

    if ((result = malloc(size1 + size2)) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n",
		progname, strerror(errno));
	exit(3);
    }
    matrix = (double complex *)((char *)result + size1);
    for (int findex = 0; findex < frequencies; ++findex) {
	result[findex] = matrix;
	matrix += elements;
    }
    return result;
}

/*
 * import_data2: copy the data matrix from cp->vdp1 to data
 *   @cp: common argument structure
 *   @data: vector of pointers to matrix to receive data
 */
static void import_data2(const common_t *cp, double complex **data)
{
    const int ports = cp->ports;
    const vnadata_t *vdp = cp->vdp2;

    for (int findex = 0; findex < cp->frequencies; ++findex) {
	for (int row = 0; row < ports; ++row) {
	    for (int column = 0; column < ports; ++column) {
		const int cell = ports * row + column;

		data[findex][cell] = vnadata_get_cell(vdp, findex, row, column);
	    }
	}
    }
}

/*
 * fill_data1: fill cp->vdp1
 *   @cp: common argument structure
 *   @type: type of the data
 *   @data: vector of pointers (one per frequency) to data matrices
 */
static libt_result_t fill_data1(const common_t *cp,
	vnadata_parameter_type_t type, double complex **data)
{
    vnadata_t *vdp = cp->vdp1;
    const int ports = cp->ports;
    const int frequencies = cp->frequencies;
    libt_result_t result = T_SKIPPED;

    if (vnadata_init(vdp, type, ports, ports, frequencies) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (vnadata_set_frequency_vector(vdp, cp->frequency_vector) == -1) {
	result = T_FAIL;
	goto out;
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	if (vnadata_set_matrix(vdp, findex, data[findex]) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    }
    if (cp->use_fz01) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *z0_vector = cp->z01[findex];

	    if (vnadata_set_fz0_vector(vdp, findex, z0_vector) == -1) {
		result = T_FAIL;
		goto out;
	    }
	}
    } else {
	if (vnadata_set_z0_vector(vdp, cp->z01[0]) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    }
    result = T_PASS;

out:
    return result;
}

/*
 * check_data2: check that cp->vd2 has correct information
 *   @cp: common argument structure
 *   @label: context information for error messages
 *   @type: expected type of the data
 *   @reference: which references to use: 1 for z01, 2 for z02
 */
static libt_result_t check_data2(const common_t *cp, const char *label,
	vnadata_parameter_type_t type, int reference)
{
    const int frequencies = cp->frequencies;
    const int ports = cp->ports;
    const vnadata_t *vdp = cp->vdp2;
    bool use_fz0;
    double complex **z0, **a, **b;
    libt_result_t result = T_SKIPPED;

    switch (reference) {
    case 1:
	use_fz0 = cp->use_fz01;
	z0 = cp->z01;
	a = cp->a1;
	b = cp->b1;
	break;
    case 2:
	use_fz0 = cp->use_fz02;
	z0 = cp->z02;
	a = cp->a2;
	b = cp->b2;
	break;
    default:
	abort();
    }
    if (vnadata_get_type(vdp) != type) {
	(void)printf("%s: type: %d != %d\n",
		label, (int)vnadata_get_type(vdp), (int)type);
	result = T_FAIL;
	goto out;
    }
    if (vnadata_get_rows(vdp) != ports) {
	(void)printf("%s: rows: %d != %d\n",
		label, vnadata_get_rows(vdp), ports);
	result = T_FAIL;
	goto out;
    }
    if (vnadata_get_columns(vdp) != ports) {
	(void)printf("%s: columns: %d != %d\n",
		label, vnadata_get_columns(vdp), ports);
	result = T_FAIL;
	goto out;
    }
    if (vnadata_get_frequencies(vdp) != frequencies) {
	(void)printf("%s: number of frequencies: %d != %d\n",
	    label, vnadata_get_frequencies(vdp), frequencies);
	result = T_FAIL;
	goto out;
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	const double f1 = vnadata_get_frequency(vdp, findex);
	const double f2 = cp->frequency_vector[findex];

	if (f1 != f2) {
	    (void)printf("%s: frequency[%d]: %e != %e\n",
		    label, findex, f1, f2);
	    result = T_FAIL;
	    goto out;
	}
    }
    if (frequencies > 1 && vnadata_has_fz0(vdp) != use_fz0) {
	(void)printf("%s: has_fz0: %s != %s\n", label,
		use_fz0 ? "true" : "false",
		vnadata_has_fz0(vdp) ? "true" : "false");
	result = T_FAIL;
	goto out;
    }
    if (use_fz0) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int port = 0; port < ports; ++port) {
		double complex v1 = vnadata_get_fz0(vdp, findex, port);
		double complex v2 = z0[findex][port];

		if (v1 != v2) {
		    (void)printf("%s: fz0[%d][%d]: %-f%+fj != %-f%+fj\n",
			    label, findex, port,
			    creal(v1), cimag(v1),
			    creal(v2), cimag(v2));
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
    } else {
	for (int port = 0; port < ports; ++port) {
	    double complex v1 = vnadata_get_z0(vdp, port);
	    double complex v2 = z0[0][port];

	    if (v1 != v2) {
		(void)printf("%s: z0[%d]: %-f%+fj != %-f%+fj\n",
			label, port,
			creal(v1), cimag(v1),
			creal(v2), cimag(v2));
		result = T_FAIL;
		goto out;
	    }
	}
    }
    switch (type) {
    case VPT_S:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *s = vnadata_get_matrix(vdp, findex);
	    double complex q[ports][ports];

	    _vnacommon_mmultiply(*q, s, a[findex], ports, ports, ports);
	    result = check_cmatrix(label, *q, b[findex], ports);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_Z:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *z = vnadata_get_matrix(vdp, findex);
	    double complex q[ports][ports];

	    _vnacommon_mmultiply(*q, z, cp->i[findex], ports, ports, ports);
	    result = check_cmatrix(label, *q, cp->v[findex], ports);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_Y:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *y = vnadata_get_matrix(vdp, findex);
	    double complex q[ports][ports];

	    _vnacommon_mmultiply(*q, y, cp->v[findex], ports, ports, ports);
	    result = check_cmatrix(label, *q, cp->i[findex], ports);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_T:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *t = vnadata_get_matrix(vdp, findex);
	    double complex x[2][2], y[2][2], q[2][2];

	    for (int k = 0; k < 2; ++k) {
		y[0][k] = b[findex][0 + k];
		y[1][k] = a[findex][0 + k];
		x[0][k] = a[findex][2 + k];
		x[1][k] = b[findex][2 + k];
	    }
	    _vnacommon_mmultiply(*q, t, *x, 2, 2, 2);
	    result = check_cmatrix(label, *q, *y, 2);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_U:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *u = vnadata_get_matrix(vdp, findex);
	    double complex x[2][2], y[2][2], q[2][2];

	    for (int k = 0; k < 2; ++k) {
		y[0][k] = a[findex][2 + k];
		y[1][k] = b[findex][2 + k];
		x[0][k] = b[findex][0 + k];
		x[1][k] = a[findex][0 + k];
	    }
	    _vnacommon_mmultiply(*q, u, *x, 2, 2, 2);
	    result = check_cmatrix(label, *q, *y, 2);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_H:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *h = vnadata_get_matrix(vdp, findex);
	    double complex x[2][2], y[2][2], q[2][2];

	    for (int k = 0; k < 2; ++k) {
		y[0][k] = cp->v[findex][0 + k];
		y[1][k] = cp->i[findex][2 + k];
		x[0][k] = cp->i[findex][0 + k];
		x[1][k] = cp->v[findex][2 + k];
	    }
	    _vnacommon_mmultiply(*q, h, *x, 2, 2, 2);
	    result = check_cmatrix(label, *q, *y, 2);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_G:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *g = vnadata_get_matrix(vdp, findex);
	    double complex x[2][2], y[2][2], q[2][2];

	    for (int k = 0; k < 2; ++k) {
		y[0][k] = cp->i[findex][0 + k];
		y[1][k] = cp->v[findex][2 + k];
		x[0][k] = cp->v[findex][0 + k];
		x[1][k] = cp->i[findex][2 + k];
	    }
	    _vnacommon_mmultiply(*q, g, *x, 2, 2, 2);
	    result = check_cmatrix(label, *q, *y, 2);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_A:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *a = vnadata_get_matrix(vdp, findex);
	    double complex x[2][2], y[2][2], q[2][2];

	    for (int k = 0; k < 2; ++k) {
		y[0][k] =  cp->v[findex][0 + k];
		y[1][k] =  cp->i[findex][0 + k];
		x[0][k] =  cp->v[findex][2 + k];
		x[1][k] = -cp->i[findex][2 + k];
	    }
	    _vnacommon_mmultiply(*q, a, *x, 2, 2, 2);
	    result = check_cmatrix(label, *q, *y, 2);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    case VPT_B:
	for (int findex = 0; findex < frequencies; ++findex) {
	    const double complex *b = vnadata_get_matrix(vdp, findex);
	    double complex x[2][2], y[2][2], q[2][2];

	    for (int k = 0; k < 2; ++k) {
		y[0][k] =  cp->v[findex][2 + k];
		y[1][k] = -cp->i[findex][2 + k];
		x[0][k] =  cp->v[findex][0 + k];
		x[1][k] =  cp->i[findex][0 + k];
	    }
	    _vnacommon_mmultiply(*q, b, *x, 2, 2, 2);
	    result = check_cmatrix(label, *q, *y, 2);
	    if (result != T_PASS) {
		goto out;
	    }
	}
	break;

    default:
	assert(!"case not handled");
    }
    result = T_PASS;

out:
    return result;
}

/*
 * type_list: list of types to convert from/to
 */
static vnadata_parameter_type_t type_list[] = {
    VPT_S,
    VPT_Z,
    VPT_Y,
    /* start 2x2 */
    VPT_T,
    VPT_U,
    VPT_H,
    VPT_G,
    VPT_A,
    VPT_B,
    VPT_UNDEF
};

/*
 * test_conversion: test conversion from type1 to each other type
 *   @cp: common argument structure
 *   @type1: source type
 *   @data: source data
 */
static libt_result_t test_conversion(common_t *cp,
	vnadata_parameter_type_t type1, double complex **data)
{
    const int frequencies = cp->frequencies;
    const int ports = cp->ports;
    char label[100];
    libt_result_t result = T_SKIPPED;

    for (int i = 0; type_list[i] != VPT_UNDEF; ++i) {
	vnadata_parameter_type_t type2 = type_list[i];
	const int z0_length = (cp->use_fz02 ? frequencies : 1) * ports;

	if (type2 == VPT_T && ports != 2) {
	    break;
	}
	if ((result = fill_data1(cp, type1, data)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(cp->vdp1, cp->vdp2, type2,
		    cp->z02[0], z0_length) == -1) {
	    goto out;
	}
	(void)sprintf(label, "%s -> %s",
		vnadata_get_type_name(type1),
		vnadata_get_type_name(type2));
	if ((result = check_data2(cp, label, type2, 2)) != T_PASS) {
	    goto out;
	}
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		const double complex *mp = vnadata_get_matrix(cp->vdp2, findex);

		(void)sprintf(label, "%s[%d]",
			vnadata_get_type_name(type2), findex);
		libt_print_cmatrix(label, mp, ports, ports);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}
    }
    result = T_PASS;

out:
    return result;
}

/*
 * error_fn: error reporting function
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)printf("error: %s: %s\n", progname, message);
}

/*
 * run_trial: test one trial of vnadata_rconvert()
 */
static libt_result_t run_test(int ports, bool use_fz01, bool use_fz02,
	bool inplace)
{
    const int frequencies = FREQUENICES;
    common_t c;
    double complex **s1 = NULL;
    double complex **t1 = NULL;
    double complex **u1 = NULL;
    double complex **z = NULL;
    double complex **y = NULL;
    double complex **h = NULL;
    double complex **g = NULL;
    double complex **a = NULL;
    double complex **b = NULL;
    char label[100];
    libt_result_t result = T_SKIPPED;

    /*
     * Init common arguments.
     */
    (void)memset((void *)&c, 0, sizeof(c));
    c.frequencies = frequencies;
    c.ports = ports;
    c.use_fz01 = use_fz01;
    c.use_fz02 = use_fz02;

    /*
     * Init the frequency vector.
     */
    if ((c.frequency_vector = calloc(frequencies, sizeof(double))) == NULL) {
	exit(3);
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	c.frequency_vector[findex] = findex * 1.0e+9;
    }

    /*
     * Generate z01 and z02.
     */
    if (use_fz01) {
	c.z01 = alloc_matrix_vector(frequencies, ports);
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int port = 0; port < ports; ++port) {
		c.z01[findex][port] = libt_crandn();
	    }
	    if (opt_v > 1) {
		(void)sprintf(label, "fz01[%d]", findex);
		libt_print_cmatrix(label, c.z01[findex], ports, 1);
	    }
	}
    } else {
	c.z01 = alloc_matrix_vector(1, ports);
	for (int port = 0; port < ports; ++port) {
	    c.z01[0][port] = libt_crandn();
	}
	if (opt_v > 1) {
	    libt_print_cmatrix("z01", c.z01[0], ports, 1);
	}
    }
    if (use_fz02) {
	c.z02 = alloc_matrix_vector(frequencies, ports);
	for (int findex = 0; findex < frequencies; ++findex) {
	    for (int port = 0; port < ports; ++port) {
		c.z02[findex][port] = libt_crandn();
	    }
	    if (opt_v > 1) {
		(void)sprintf(label, "fz02[%d]", findex);
		libt_print_cmatrix(label, c.z02[findex], ports, 1);
	    }
	}
    } else {
	c.z02 = alloc_matrix_vector(1, ports);
	for (int port = 0; port < ports; ++port) {
	    c.z02[0][port] = libt_crandn();
	}
	if (opt_v > 1) {
	    libt_print_cmatrix("z02", c.z02[0], ports, 1);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * Generate s1.
     */
    s1 = alloc_matrix_vector(frequencies, ports * ports);
    for (int findex = 0; findex < frequencies; ++findex) {
	for (int row = 0; row < ports; ++row) {
	    for (int column = 0; column < ports; ++column) {
		int cell = ports * row + column;

		s1[findex][cell] = libt_crandn();
	    }
	}
	if (opt_v > 1) {
	    (void)sprintf(label, "s1[%d]", findex);
	    libt_print_cmatrix(label, s1[findex], ports, ports);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * Generate a1 and b1.
     */
    c.a1 = alloc_matrix_vector(frequencies, ports * ports);
    c.b1 = alloc_matrix_vector(frequencies, ports * ports);
    for (int findex = 0; findex < frequencies; ++findex) {
	for (int row = 0; row < ports; ++row) {
	    for (int column = 0; column < ports; ++column) {
		int cell = ports * row + column;

		c.a1[findex][cell] = libt_crandn();
	    }
	}
	_vnacommon_mmultiply(c.b1[findex], s1[findex], c.a1[findex],
		ports, ports, ports); // b1 = S * a1
	if (opt_v > 1) {
	    (void)sprintf(label, "a1[%d]", findex);
	    libt_print_cmatrix(label, c.a1[findex], ports, ports);
	    (void)sprintf(label, "b1[%d]", findex);
	    libt_print_cmatrix(label, c.b1[findex], ports, ports);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * Find v and i.
     */
    c.v = alloc_matrix_vector(frequencies, ports * ports);
    c.i = alloc_matrix_vector(frequencies, ports * ports);
    for (int findex = 0; findex < frequencies; ++findex) {
	const int z0_findex = use_fz01 ? findex : 0;

	for (int row = 0; row < ports; ++row) {
	    double complex z0 = c.z01[z0_findex][row];
	    double complex z0c = conj(z0);
	    double r = creal(z0);
	    double k = sqrt(fabs(r)) / r;
	    for (int column = 0; column < ports; ++column) {
		int cell = ports * row + column;
		const double complex a1 = c.a1[findex][cell];
		const double complex b1 = c.b1[findex][cell];

		c.v[findex][cell] = k * (z0c * a1 + z0 * b1);
		c.i[findex][cell] = k * (a1 - b1);
	    }
	}
	if (opt_v > 1) {
	    (void)sprintf(label, "v[%d]", findex);
	    libt_print_cmatrix(label, c.v[findex], ports, ports);
	    (void)sprintf(label, "i[%d]", findex);
	    libt_print_cmatrix(label, c.i[findex], ports, ports);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * Find a1 and a2.
     */
    c.a2 = alloc_matrix_vector(frequencies, ports * ports);
    c.b2 = alloc_matrix_vector(frequencies, ports * ports);
    for (int findex = 0; findex < frequencies; ++findex) {
	const int z0_findex = use_fz02 ? findex : 0;

	for (int row = 0; row < ports; ++row) {
	    double complex z0 = c.z02[z0_findex][row];
	    double complex z0c = conj(z0);
	    double r = creal(z0);
	    double k = 1.0 / (2.0 * sqrt(fabs(r)));
	    for (int column = 0; column < ports; ++column) {
		int cell = ports * row + column;
		double complex v = c.v[findex][cell];
		double complex i = c.i[findex][cell];

		c.a2[findex][cell] = k * (v + z0  * i);
		c.b2[findex][cell] = k * (v - z0c * i);
	    }
	}
	if (opt_v > 1) {
	    (void)sprintf(label, "a2[%d]", findex);
	    libt_print_cmatrix(label, c.a2[findex], ports, ports);
	    (void)sprintf(label, "b2[%d]", findex);
	    libt_print_cmatrix(label, c.b2[findex], ports, ports);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * Allocate vdp1 and vdp2.
     */
    if ((c.vdp1 = vnadata_alloc(error_fn, NULL)) == NULL) {
	exit(3);
    }
    if (inplace) {
	c.vdp2 = c.vdp1;
    } else {
	if ((c.vdp2 = vnadata_alloc(error_fn, NULL)) == NULL) {
	    exit(3);
	}
    }

    /*
     * Find z
     */
    z = alloc_matrix_vector(frequencies, ports * ports);
    if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	goto out;
    }
    if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_Z, NULL, 0) == -1) {
	goto out;
    }
    if ((result = check_data2(&c, "s1 -> z1", VPT_Z, 1)) != T_PASS) {
	goto out;
    }
    import_data2(&c, z);
    if (opt_v > 1) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    (void)sprintf(label, "z[%d]", findex);
	    libt_print_cmatrix(label, z[findex], ports, ports);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * Find y
     */
    y = alloc_matrix_vector(frequencies, ports * ports);
    if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	goto out;
    }
    if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_Y, NULL, 0) == -1) {
	goto out;
    }
    if ((result = check_data2(&c, "s1 -> y1", VPT_Y, 1)) != T_PASS) {
	goto out;
    }
    import_data2(&c, y);
    if (opt_v > 1) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    (void)sprintf(label, "y[%d]", findex);
	    libt_print_cmatrix(label, y[findex], ports, ports);
	}
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }

    /*
     * If 2-port, find t1, u1, g, h, a and b.
     */
    if (ports == 2) {
	/*
	 * Find t1
	 */
	t1 = alloc_matrix_vector(frequencies, ports * ports);
	if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_T, NULL, 0) == -1) {
	    goto out;
	}
	if ((result = check_data2(&c, "s1 -> t1", VPT_T, 1)) != T_PASS) {
	    goto out;
	}
	import_data2(&c, t1);
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		(void)sprintf(label, "s1 -> t1[%d]", findex);
		libt_print_cmatrix(label, t1[findex], 2, 2);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}

	/*
	 * Find u1
	 */
	u1 = alloc_matrix_vector(frequencies, ports * ports);
	if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_U, NULL, 0) == -1) {
	    goto out;
	}
	if ((result = check_data2(&c, "s1 -> u1", VPT_U, 1)) != T_PASS) {
	    goto out;
	}
	import_data2(&c, u1);
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		(void)sprintf(label, "s1 -> u1[%d]", findex);
		libt_print_cmatrix(label, u1[findex], 2, 2);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}

	/*
	 * Find g
	 */
	g = alloc_matrix_vector(frequencies, ports * ports);
	if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_G, NULL, 0) == -1) {
	    goto out;
	}
	if ((result = check_data2(&c, "s1 -> g", VPT_G, 1)) != T_PASS) {
	    goto out;
	}
	import_data2(&c, g);
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		(void)sprintf(label, "s1 -> g[%d]", findex);
		libt_print_cmatrix(label, g[findex], 2, 2);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}

	/*
	 * Find h
	 */
	h = alloc_matrix_vector(frequencies, ports * ports);
	if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_H, NULL, 0) == -1) {
	    goto out;
	}
	if ((result = check_data2(&c, "s1 -> h", VPT_H, 1)) != T_PASS) {
	    goto out;
	}
	import_data2(&c, h);
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		(void)sprintf(label, "s1 -> h[%d]", findex);
		libt_print_cmatrix(label, h[findex], 2, 2);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}

	/*
	 * Find a
	 */
	a = alloc_matrix_vector(frequencies, ports * ports);
	if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_A, NULL, 0) == -1) {
	    goto out;
	}
	if ((result = check_data2(&c, "s1 -> a", VPT_A, 1)) != T_PASS) {
	    goto out;
	}
	import_data2(&c, a);
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		(void)sprintf(label, "s1 -> a[%d]", findex);
		libt_print_cmatrix(label, a[findex], 2, 2);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}

	/*
	 * Find b
	 */
	b = alloc_matrix_vector(frequencies, ports * ports);
	if ((result = fill_data1(&c, VPT_S, s1)) != T_PASS) {
	    goto out;
	}
	if (vnadata_rconvert(c.vdp1, c.vdp2, VPT_B, NULL, 0) == -1) {
	    goto out;
	}
	if ((result = check_data2(&c, "s1 -> b", VPT_B, 1)) != T_PASS) {
	    goto out;
	}
	import_data2(&c, b);
	if (opt_v > 1) {
	    for (int findex = 0; findex < frequencies; ++findex) {
		(void)sprintf(label, "s1 -> b[%d]", findex);
		libt_print_cmatrix(label, b[findex], 2, 2);
	    }
	}
	if (opt_v > 1) {
	    (void)printf("\n");
	}
    }

    /*
     * Convert from each type.
     */
    if ((result = test_conversion(&c, VPT_S, s1)) != T_PASS) {
	goto out;
    }
    if ((result = test_conversion(&c, VPT_Z, z)) != T_PASS) {
	goto out;
    }
    if ((result = test_conversion(&c, VPT_Y, y)) != T_PASS) {
	goto out;
    }
    if (ports == 2) {
	if ((result = test_conversion(&c, VPT_T, t1)) != T_PASS) {
	    goto out;
	}
	if ((result = test_conversion(&c, VPT_U, u1)) != T_PASS) {
	    goto out;
	}
	if ((result = test_conversion(&c, VPT_H, h)) != T_PASS) {
	    goto out;
	}
	if ((result = test_conversion(&c, VPT_G, g)) != T_PASS) {
	    goto out;
	}
	if ((result = test_conversion(&c, VPT_A, a)) != T_PASS) {
	    goto out;
	}
	if ((result = test_conversion(&c, VPT_B, b)) != T_PASS) {
	    goto out;
	}
    }
    result = T_PASS;

out:
    free((void *)b);
    free((void *)a);
    free((void *)g);
    free((void *)h);
    free((void *)y);
    free((void *)z);
    free((void *)u1);
    free((void *)t1);
    free((void *)s1);
    if (c.vdp2 != c.vdp1) {
	vnadata_free(c.vdp2);
    }
    vnadata_free(c.vdp1);
    free((void *)c.b2);
    free((void *)c.a2);
    free((void *)c.i);
    free((void *)c.v);
    free((void *)c.b1);
    free((void *)c.a1);
    free((void *)c.z02);
    free((void *)c.z01);

    return result;
}

/*
 * run_trials: run all trials of the test
 */
static libt_result_t run_trials()
{
    libt_result_t result = T_SKIPPED;

    for (int trial = 1; trial <= TRIALS; ++trial) {
	for (int ports = 1; ports <= 5; ++ports) {
	    for (bool use_fz01 = false, stop = false; !stop;
		    stop = use_fz01, use_fz01 = true) {
		for (bool use_fz02 = false, stop = false; !stop;
			stop = use_fz02, use_fz02 = true) {
		    for (bool inplace = false, stop = false; !stop;
			    stop = inplace, inplace = true) {
			if (opt_v) {
			    (void)printf("Test rconvert: trial %3d ports %d "
				    "use_fz01 %d use_fz02 %d inplace %d\n",
				    trial, ports, (int)use_fz01, (int)use_fz02,
				    (int)inplace);
			}
			result = run_test(ports, use_fz01, use_fz02, inplace);
			if (result != T_PASS) {
			    goto out;
			}
		    }
		}
	    }
	}
	if (opt_v) {
	    printf("-------------\n");
	}
    }
    result = T_PASS;

out:
    libt_report(result);
    return result;
}

/*
 * print_usage: print a usage message and exit
 */
static void print_usage()
{
    const char *const *cpp;

    for (cpp = usage; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s: usage %s\n", progname, *cpp);
    }
    for (cpp = help; *cpp != NULL; ++cpp) {
	(void)fprintf(stderr, "%s\n", *cpp);
    }
    exit(2);
}

/*
 * main: test program
 */
int
main(int argc, char **argv)
{
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case 'a':
	    opt_a = 1;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	case -1:
	    break;

	default:
	    print_usage();
	}
	break;
    }
    argc -= optind;
    argv += optind;
    if (argc != 0) {
	print_usage();
    }
    libt_isequal_init();
    exit(run_trials());
}
