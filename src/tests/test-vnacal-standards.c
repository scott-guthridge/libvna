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
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <stdio.h>
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif
#include "libt_crand.h"
#include "libt.h"
#include "vnacal_internal.h"
#include "vnaconv.h"


#define NTRIALS	 	30
#define NRECTANGULAR	2

#define MIN_FPOINTS	30
#define MAX_FPOINTS    600

#define FMIN	 1.0e+9
#define FMAX	18.0e+9

/*
 * Command Line Options
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
int  opt_v = 0;

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
 * error_fn: error reporting function
 *   @message: error message
 *   @arg: (unused)
 *   @category: error category (unused)
 */
static void error_fn(const char *message, void *arg, vnaerr_category_t category)
{
    (void)printf("%s: %s\n", progname, message);
}

/*
 * tf2_t: coefficients of a 2nd order transfer function
 */
typedef struct tf2 {
    double complex tf2_n0, tf2_n1, tf2_n2;
    double complex tf2_d1, tf2_d2;
} tf2_t;

/*
 * tf2_init: create random 2nd order transfer function
 */
static void tf2_init(tf2_t *p, double fmin, double fmax)
{
    double fc = (fmin + fmax) / 2.0;
    double complex z1 = fc * libt_rand_nsmm(0.832557, 0.5, 0.01, 100.0);
    double complex z2 = fc * libt_rand_nsmm(0.832557, 0.5, 0.01, 100.0);
    double complex p1 = fc * libt_rand_nsmm(0.832557, 0.5, 0.01, 100.0);
    double complex p2 = fc * libt_rand_nsmm(0.832557, 0.5, 0.01, 100.0);
    double complex d = p1 * p2;

    /*
     * Convert random poles and zeros to ratio of polynomial form.
     * Note that we didn't force the complex poles and zeros above
     * into conjugate pairs because all we want is a mathematically
     * consistent ratio of smooth functions of f.
     */
    p->tf2_n0 =  z1 * z2   / d;
    p->tf2_n1 = -(z1 + z2) / d;
    p->tf2_n2 =  1.0       / d;
    p->tf2_d1 = -(p1 + p2) / d;
    p->tf2_d2 =  1.0       / d;
}

/*
 * tf2_eval: evaluate the transfer function at f
 */
static double complex tf2_eval(const tf2_t *p, double f)
{
    return (p->tf2_n0 + f * (p->tf2_n1 + f * p->tf2_n2)) /
	     (1.0 + f * (p->tf2_d1 + f * p->tf2_d2));
}

/*
 * tf2_print: print the transfer function coefficients
 */
static void tf2_print(const tf2_t *p)
{
    (void)printf("    n0: %+e %+ej\n", creal(p->tf2_n0), cimag(p->tf2_n0));
    (void)printf("    n1: %+e %+ej\n", creal(p->tf2_n1), cimag(p->tf2_n1));
    (void)printf("    n2: %+e %+ej\n", creal(p->tf2_n2), cimag(p->tf2_n2));
    (void)printf("    d1: %+e %+ej\n", creal(p->tf2_d1), cimag(p->tf2_d1));
    (void)printf("    d2: %+e %+ej\n", creal(p->tf2_d2), cimag(p->tf2_d2));
}

/*
 * eval_fz0: return frequency-dependent z0 value at frequency f
 *   @zc: impedance returned at center frequency (ohms)
 *   @fc: center frequency (Hz)
 *   @f: frequency
 *
 * Models z0 as a resistor in parallel with a capacitor with given
 * impedance at center frequency.
 */
static double complex eval_fz0(double complex zc, double fc, double f)
{
    double rc = creal(zc);
    double xc = cimag(zc);

    return (rc*rc + xc*xc) / (rc - I * xc * f / fc);
}

/*
 * keysight_calc_tline_parameters0: calc Zc, gl (classic version)
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
static double complex keysight_calc_tline_parameters0(
	const vnacal_calkit_data_t *vcdp, double f, double complex *Zc)
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
    *Zc = offset_z0 + (1.0 - I) * offset_loss * fGrt / (2.0 * w);

    return gamma_l;
}

/*
 * keysight_calc_tline_parameters: calc Z, gl (revised version)
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zc: address of complex to receive the characteristic impedance
 *
 * Returns the transmission coefficient times electrical length
 * (gamma el).  This is the revised version described here:
 * https://www.keysight.com/us/en/assets/7018-01375/application-notes/
 * 5989-4840.pdf
 */
static double complex keysight_calc_tline_parameters(
	const vnacal_calkit_data_t *vcdp, double f, double complex *Zc)
{
    double complex temp;
    double offset_delay = vcdp->vcd_offset_delay;	/* s */
    double offset_loss = vcdp->vcd_offset_loss;		/* Ω/s */
    double offset_z0 = vcdp->vcd_offset_z0;		/* Ω */

    temp = csqrt(1.0 +
		 (1.0 - I) * offset_loss /
		 (2.0 * M_PI * sqrt(1.0e+9 * f) * offset_z0));
    *Zc = offset_z0 * temp;
    return I * 2.0 * M_PI * f * offset_delay * temp;
}

/*
 * keysight_add_tline: find reflection coefficient of a load at end of a tline
 *   @Gt: reflection coefficient of the termination (load)
 *   @Zr: reference impedance (only real numbers)
 *   @Zc: characteristic impedance of the transmission line
 *   @gl: transmission coefficient times electrical length
 */
static double complex keysight_add_tline(double complex Gt, double Zr,
	double complex Zc, double complex gl)
{
    double complex G1 = (Zc - Zr) / (Zc + Zr);
    double complex em2gl = cexp(-2.0 * gl);

    return (G1 * (1.0 - em2gl - G1 * Gt) + em2gl * Gt) /
	   (1.0 - G1 * (em2gl * G1 + Gt * (1.0 - em2gl)));
}

/*
 * keysight_short: evaluate a calkit short standard at f, Zr
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zr: reference impedance (real, positive)
 */
static double complex keysight_short(const vnacal_calkit_data_t *vcdp,
	double f, double Zr)
{
    double Ls = vcdp->vcd_l_coefficients[0] +
		vcdp->vcd_l_coefficients[1] * f +
		vcdp->vcd_l_coefficients[2] * f*f +
		vcdp->vcd_l_coefficients[3] * f*f*f;
    double complex Zt = I * 2.0 * M_PI * f * Ls;
    double complex Gt = (Zt - Zr) / (Zt + Zr);
    double complex Zc;
    double complex gl;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = keysight_calc_tline_parameters0(vcdp, f, &Zc);
    } else {
	gl = keysight_calc_tline_parameters(vcdp, f, &Zc);
    }
    return keysight_add_tline(Gt, Zr, Zc, gl);
}

/*
 * keysight_open: evaluate a calkit open standard at f, Zr
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zr: reference impedance (real, positive)
 */
static double complex keysight_open(const vnacal_calkit_data_t *vcdp,
	double f, double Zr)
{
    double Co = vcdp->vcd_c_coefficients[0] +
		vcdp->vcd_c_coefficients[1] * f +
		vcdp->vcd_c_coefficients[2] * f*f +
		vcdp->vcd_c_coefficients[3] * f*f*f;
    double complex Zt = 1.0 / (I * 2.0 * M_PI * f * Co);
    double complex Gt = (Zt - Zr) / (Zt + Zr);
    double complex Zc;
    double complex gl;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = keysight_calc_tline_parameters0(vcdp, f, &Zc);
    } else {
	gl = keysight_calc_tline_parameters(vcdp, f, &Zc);
    }
    return keysight_add_tline(Gt, Zr, Zc, gl);
}

/*
 * keysight_load: evaluate a calkit load standard at f, Zr
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zr: reference impedance (real, positive)
 */
static double complex keysight_load(const vnacal_calkit_data_t *vcdp,
	double f, double Zr)
{
    double complex Zt = vcdp->vcd_zl;
    double complex Gt = (Zt - Zr) / (Zt + Zr);
    double complex Zc;
    double complex gl;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = keysight_calc_tline_parameters0(vcdp, f, &Zc);
    } else {
	gl = keysight_calc_tline_parameters(vcdp, f, &Zc);
    }
    return keysight_add_tline(Gt, Zr, Zc, gl);
}

/*
 * keysight_through: evaluate a calkit through standard at f, Zr
 *   @vcdp: calibration kit data
 *   @f: frequency in Hz
 *   @Zr: reference impedance (real, positive)
 *   @result: caller-allocated vector to receive result
 */
static void keysight_through(const vnacal_calkit_data_t *vcdp,
	double f, double Zr, double complex (*result)[2])
{
    double complex Zc;
    double complex gl;
    double complex G;
    double complex p;
    double complex d;
    double complex s11;
    double complex s12;

    if (vcdp->vcd_flags & VNACAL_CKF_TRADITIONAL) {
	gl = keysight_calc_tline_parameters0(vcdp, f, &Zc);
    } else {
	gl = keysight_calc_tline_parameters(vcdp, f, &Zc);
    }
    G = (Zc - Zr) / (Zc + Zr);
    p = cexp(-gl);
    d = 1.0 - p * p * G * G;
    s11 = G * (1.0 - p * p) / d;
    s12 = p * (1.0 - G * G) / d;
    result[0][0] = s11;
    result[0][1] = s12;
    result[1][0] = s12;
    result[1][1] = s11;
}

/*
 * compute_l_coefficients: fit an inductor in parallel with small capacitor
 *     Find a cubic polynomial fit for the inductance of an inductor (l)
 *     in parallel with a small capacitor (c) over a frequency range:
 *     f1..f2.  The capacitor must be small enough that its resonant
 *     frequency with the inductor lies beyond f2.
 *
 *     If l and c are reversed, the same function finds the dual polynomial
 *     of capacitance for a capacitor in series with a small inductor over
 *     the frequency range.
 *
 * == Theory ==
 *
 * Impedance of an inductor in parallel with a capacitor is:
 *
 *        j ω L           j ω L
 *   --------------- = -----------
 *   1 + (j ω)^2 L C   1 - ω^2 L C
 *
 * If we set this equal to the effective impedance of an inductor, L_eff(ω):
 *
 *      j ω L
 *   ----------- = j ω L_eff(ω)
 *   1 - ω^2 L C
 *
 * and solve for L_eff(ω), we have:
 *
 *                  L
 *   L_eff(ω) = -----------
 *              1 - ω^2 L C
 *
 * substituting: ω = 2 pi f:
 *
 *                      L
 *   L_eff(f) = ------------------
 *              1 - (2 pi f)^2 L C
 *
 * we want polynomial coefficients L0, L1, L2 and L3 that minimize
 * the squared error:
 *
 *   Integral f1..f2 of (L0 + L1 f + L2 f^2 + L3 f^3 - L_eff(f))^2 df
 *
 * We can take that integral then take the derivatives with respect to
 * each of L0, L1, L2 and L3.  Setting those derivatives to zero forms a
 * linear system of four equations and four unknowns.  But this can be
 * expressed more compactly as a least-squares minimization problem in
 * the Hilbert space L^2([f1, f2]) with normal equations:
 *
 *   Sum_{j=0}^3 <f^i, f^j> L[j] = <f^i, L_eff>, i = 0..3
 *
 * with coefficient matrix:
 *
 *   A[i,j] = Integral f1..f2 of f^(i+j) df, i = 0..3, j = 0..3
 *
 *          = 1 / k (f2^k - f1^k), k = i + j + 1, i = 0..3, j = 0..3
 *
 * and right-hand side vector:
 *
 *   b[i] = Integral f1..f2 of f^i Leff(f) df, i = 0..3
 *
 * Taking the integrals above, we have:
 *
 * A:
 *  (f2 - f1),         (f2^2 - f1^2) / 2, (f2^3 - f1^3) / 3, (f2^4 - f1^4) / 4
 *  (f2^2 - f1^2) / 2, (f2^3 - f1^3) / 3, (f2^4 - f1^4) / 4, (f2^5 - f1^5) / 5
 *  (f2^3 - f1^3) / 3, (f2^4 - f1^4) / 4, (f2^5 - f1^5) / 5, (f2^6 - f1^6) / 6
 *  (f2^4 - f1^4) / 4, (f2^5 - f1^5) / 5, (f2^6 - f1^6) / 6, (f2^7 - f1^7) / 7
 *
 * b:
 *   L atanh(2 pi f sqrt(L C))
 *   -------------------------
 *         2 pi sqrt(L C)
 *
 *   -log(1 - 4 pi^2 L C f^2)
 *   -------------------------
 *           9 C pi^2
 *
 *   atanh(2 pi sqrt(L C) f) - 2 pi f sqrt(L C)
 *   ------------------------------------------
 *                8 pi^3 C sqrt(L C)
 *
 *  -log(1 - 4 pi^2 C f^2) - 4 pi^2 L C f^2
 *  ---------------------------------------
 *                32 pi^4 L C^2
 *
 * over the interval f = f1,f2.
 *
 * The A matrix is positive-definite so can be decomposed using Cholesky
 * decomposition, making the system easily solvable using forward and
 * backward substitution.
 */
static void compute_l_coefficients(double f1, double f2, double l, double c,
	double result[4])
{
    assert(f1 > 0.0);
    assert(f2 > f1);
    assert(l > 0.0);
    assert(c >= 0.0);

    /*
     * Special-case c of zero to avoid divide by zero below.
     */
    if (c == 0.0) {
	result[0] = l;
	for (int i = 1; i < 4; ++i) {
	    result[i] = 0.0;
	}
	return;
    }
    {
	const double f1_2 = f1 * f1;
	const double f1_3 = f1_2 * f1;
	const double f1_4 = f1_3 * f1;
	const double f2_2 = f2 * f2;
	const double f2_3 = f2_2 * f2;
	const double f2_4 = f2_3 * f2;
	const double df = f2 - f1;
	const double df_12 = sqrt(df);
	const double df_32 = pow(df, 3.0 / 2.0);
	const double df_52 = pow(df, 5.0 / 2.0);
	const double df_72 = pow(df, 7.0 / 2.0);
	const double sqrt3 = sqrt(3.0);
	const double sqrt5 = sqrt(5.0);
	const double sqrt7 = sqrt(7.0);
	const double sqrt_lc = sqrt(l * c);
	const double atanh_temp1 = atanh(2.0 * M_PI * sqrt_lc * f1);
	const double atanh_temp2 = atanh(2.0 * M_PI * sqrt_lc * f2);
	const double log_temp1 = log(1.0 - 4.0 * M_PI * M_PI * l * c * f1_2);
	const double log_temp2 = log(1.0 - 4.0 * M_PI * M_PI * l * c * f2_2);
	double u[4][4] = {	/* Cholesky decomposition of A */
	    {
		df_12,
		(f2_2 - f1_2) / (2.0 * df_12),
		(f2_3 - f1_3) / (3.0 * df_12),
		(f2_4 - f1_4) / (4.0 * df_12)
	    },
	    {
		0.0,
		df_32 / (2.0 * sqrt3),
		df_32 * (f1 + f2) / (2.0 * sqrt3),
		sqrt3 * df_32 * (3.0 * f1_2 + 4.0 * f1 * f2 + 3.0 * f2_2) / 20.0
	    },
	    {
		0.0,
		0.0,
		df_52 / (6.0 * sqrt5),
		df_52 * (f1 + f2) / (4.0 * sqrt5)
	    },
	    {
		0.0,
		0.0,
		0.0,
		df_72 / (20.0 * sqrt7)
	    }
	};
	double b[4] = {
	    l * (atanh_temp2 - atanh_temp1) / (2.0 * M_PI * sqrt_lc),
	    -(log_temp2 - log_temp1) / (8.0 * M_PI * M_PI * c),
	    (atanh_temp2 - atanh_temp1 - 2.0 * M_PI * sqrt_lc * (f2 - f1)) /
		(8.0 * M_PI * M_PI * M_PI * c * sqrt_lc),
	    -(log_temp2 - log_temp1 +
		    4.0 * M_PI * M_PI * l * c * (f2_2 - f1_2)) /
		(32.0 * M_PI * M_PI * M_PI * M_PI * l * c * c)
	};

	/*
	 * Normalize the rows of the lower-triangular matrix and b so
	 * that the major diagonal (not stored) is all 1's.  Store the
	 * L entries in the lower triangle of u.
	 */
	for (int i = 0; i < 4; ++i) {
	    double d = u[i][i];

	    for (int j = 0; j < i; ++j) {
		u[i][j] = u[j][i] / d;
	    }
	    b[i] /= d;
	}

	/*
	 * Use forward substitution to find the intermediate X' such that
	 * L X' = B.
	 */
	for (int i = 0; i < 4; ++i) {
	    double s = b[i];

	    for (int k = 0; k < i; ++k) {
		s -= u[i][k] * result[k];
		assert(i > k);
	    }
	    result[i] = s;
	}

	/*
	 * Use back substitution to find the result X, such that U X = X'.
	 */
	for (int i = 3; i >= 0; --i) {
	    double s = result[i];

	    for (int k = i + 1; k < 4; ++k) {
		s -= u[i][k] * result[k];
	    }
	    result[i] = s / u[i][i];
	}
    }
}

/*
 * standard_type_t: type of standard
 */
typedef enum standard_type {
    T_TRADITIONAL_SCALAR,
    T_TRADITIONAL_VECTOR,
    T_CALKIT,
    T_DATA,
} standard_type_t;

/*
 * test_standard_t: describes any standard
 */
typedef struct test_standard {
    struct test   *ts_tp;	/* associated test test */
    standard_type_t	  ts_type;	/* type of standard */
    int			  ts_ports;	/* number of ports in standard */
    union {
	int		 *port_map;	/* std port -> parameter matrix port */
	int               cell;		/* parameter matrix cell (trad) */
    } u1;
    struct test_standard *ts_next;	/* next standard in list */
    union {
	double complex scalar;		/* constant value */
	tf2_t vector;			/* transfer function to compute vals */
	vnacal_calkit_data_t calkit;	/* calkit parameters */
	struct {
	    bool has_fz0;		/* frequency-dependent z0's */
	    double complex *z0_vector;	/* z0 of port or z0 at fc of port */
	    tf2_t *matrix;		/* matrix of transfer functions */
	} data;
    } u2;
} test_standard_t;
#define ts_port_map	u1.port_map
#define ts_cell		u1.cell
#define ts_scalar	u2.scalar
#define ts_vector	u2.vector
#define ts_calkit	u2.calkit
#define ts_has_fz0	u2.data.has_fz0
#define ts_z0_vector	u2.data.z0_vector
#define ts_matrix	u2.data.matrix

/*
 * test_t: common state for each test test
 */
typedef struct test {
    vnacal_t           *t_vcp;		/* vnacal_t structure */
    int                 t_rows;		/* rows in parameter matrix */
    int			t_columns;	/* columns in parameter matrix */
    int			t_remaining;	/* count of remaining ports */
    int		       *t_port_vector;	/* vector of unused/used ports */
    uint32_t		t_traditional;	/* bitmask of traditional ports */
    int                *t_parameter_matrix; /* matrix of parameter handles */
    double              t_fmin;		/* minimum frequency */
    double		t_fmax;		/* maximum frequency */
    test_standard_t    *t_head;		/* list of standards */
    test_standard_t   **t_anchor;	/* where to insert next standard */
} test_t;

/*
 * test_standard_alloc: allocate and init a test_standard_t structure
 *   @type: type of standard
 */
test_standard_t *test_standard_alloc(test_t *tp, standard_type_t type)
{
    test_standard_t *tsp = NULL;

    if ((tsp = malloc(sizeof(test_standard_t))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	exit(T_ERROR);
    }
    (void)memset((void *)tsp, 0, sizeof(*tsp));
    tsp->ts_tp = tp;
    tsp->ts_type = type;
    return tsp;
}

/*
 * test_standard_free: free a test standard
 *   @tsp: standard struct to free
 */
void test_standard_free(test_standard_t *tsp)
{
    if (tsp != NULL) {
	switch (tsp->ts_type) {
	case T_TRADITIONAL_SCALAR:
	case T_TRADITIONAL_VECTOR:
	    break;
	case T_CALKIT:
	    free((void *)tsp->ts_port_map);
	    break;
	case T_DATA:
	    free((void *)tsp->ts_port_map);
	    free((void *)tsp->ts_z0_vector);
	    free((void *)tsp->ts_matrix);
	    break;
	}
	free((void *)tsp);
    }
}

/*
 * test_standard_print: print a test standard
 *   @tsp: standard struct to print
 */
void test_standard_print(test_standard_t *tsp)
{
    test_t *tp = tsp->ts_tp;
    int columns = tp->t_columns;

    switch (tsp->ts_type) {
    case T_TRADITIONAL_SCALAR:
	{
	    int row = tsp->ts_cell / columns;
	    int column = tsp->ts_cell % columns;
	    (void)printf("standard: traditional scalar\n");
	    (void)printf("    s%d%d\n", row + 1, column + 1);
	    (void)printf("    %f%+fj\n",
		    creal(tsp->ts_scalar), cimag(tsp->ts_scalar));
	    (void)printf("\n");
	}
	break;

    case T_TRADITIONAL_VECTOR:
	{
	    int row = tsp->ts_cell / columns;
	    int column = tsp->ts_cell % columns;
	    (void)printf("standard: traditional vector\n");
	    (void)printf("    s%d%d\n", row + 1, column + 1);
	    tf2_print(&tsp->ts_vector);
	    (void)printf("\n");
	}
	break;

    case T_CALKIT:
	{
	    vnacal_calkit_data_t *vcdp = &tsp->ts_calkit;
	    const char *subtype;

	    switch (vcdp->vcd_type) {
	    case VNACAL_CALKIT_SHORT:
		subtype = "short";
		break;
	    case VNACAL_CALKIT_OPEN:
		subtype = "open";
		break;
	    case VNACAL_CALKIT_LOAD:
		subtype = "load";
		break;
	    case VNACAL_CALKIT_THROUGH:
		subtype = "through";
		break;
	    default:
		abort();
	    }
	    (void)printf("standard: calkit %s\n", subtype);
	    for (int r = 0; r < tsp->ts_ports; ++r) {
		int row = tsp->ts_port_map[r];

		(void)printf("   ");
		for (int c = 0; c < tsp->ts_ports; ++c) {
		    int column = tsp->ts_port_map[c];

		    (void)printf(" s%d%d", row + 1, column + 1);
		}
		(void)printf("\n");
	    }
	    (void)printf("\n");
	    (void)printf("    offset delay: %e\n", vcdp->vcd_offset_delay);
	    (void)printf("    offset loss:  %e\n", vcdp->vcd_offset_loss);
	    (void)printf("    offset z0:    %f\n", vcdp->vcd_offset_z0);
	    switch (vcdp->vcd_type) {
	    case VNACAL_CALKIT_SHORT:
		(void)printf("    l0: %e\n", vcdp->vcd_l_coefficients[0]);
		(void)printf("    l1: %e\n", vcdp->vcd_l_coefficients[1]);
		(void)printf("    l2: %e\n", vcdp->vcd_l_coefficients[2]);
		(void)printf("    l3: %e\n", vcdp->vcd_l_coefficients[3]);
		break;

	    case VNACAL_CALKIT_OPEN:
		(void)printf("    c0: %e\n", vcdp->vcd_c_coefficients[0]);
		(void)printf("    c1: %e\n", vcdp->vcd_c_coefficients[1]);
		(void)printf("    c2: %e\n", vcdp->vcd_c_coefficients[2]);
		(void)printf("    c3: %e\n", vcdp->vcd_c_coefficients[3]);
		break;

	    case VNACAL_CALKIT_LOAD:
		(void)printf("    zl: %f%+fj\n",
			creal(vcdp->vcd_zl), cimag(vcdp->vcd_zl));
		break;

	    default:
	        break;
	    }
	    (void)printf("\n");
	}
	break;

    case T_DATA:
	(void)printf("standard: data %d x %d\n", tsp->ts_ports, tsp->ts_ports);
	for (int r = 0; r < tsp->ts_ports; ++r) {
	    int row = tsp->ts_port_map[r];

	    (void)printf("   ");
	    for (int c = 0; c < tsp->ts_ports; ++c) {
		int column = tsp->ts_port_map[c];

		(void)printf(" s%d%d", row + 1, column + 1);
	    }
	    (void)printf("\n");
	}
	(void)printf("    fz0: %s\n", tsp->ts_has_fz0 ? "true" : "false");
	for (int i = 0; i < tsp->ts_ports; ++i) {
	    (void)printf("    z0[%d]: %f%+fj\n",
		i, creal(tsp->ts_z0_vector[i]), cimag(tsp->ts_z0_vector[i]));
	}
	(void)printf("\n");
	for (int r = 0; r < tsp->ts_ports; ++r) {
	    for (int c = 0; c < tsp->ts_ports; ++c) {
		int cell = r * tsp->ts_ports + c;

		(void)printf("  s%d%d:\n", r + 1, c + 1);
		tf2_print(&tsp->ts_matrix[cell]);
		(void)printf("\n");
	    }
	}
	break;

    default:
        abort();
    }
}

/*
 * test_init: init a caller-allocated test_t struct
 *   @tp: test common state
 *   @rows: rows in parameter matrix
 *   @columns: columns in parameter matrix
 *   @fmin: minimum frequency
 *   @fmax: maximum frequency
 */
static libt_result_t test_init(test_t *tp, int rows, int columns,
	double fmin, double fmax)
{
    int ports = MAX(rows, columns);
    vnacal_t *vcp = NULL;
    libt_result_t result = T_FAIL;

    /*
     * Create the calibration structure.
     */
    if ((vcp = vnacal_create(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnacal_create: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }

    /*
     * Init the test structure.
     */
    (void)memset((void *)tp, 0, sizeof(*tp));
    tp->t_vcp = vcp;
    tp->t_rows = rows;
    tp->t_columns = columns;
    tp->t_remaining = ports;
    if ((tp->t_port_vector = malloc(ports * sizeof(int))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	exit(T_ERROR);
    }
    for (int i = 0; i < ports; ++i) {
	tp->t_port_vector[i] = i;
    }
    if ((tp->t_parameter_matrix = calloc(rows * columns,
		    sizeof(int))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	exit(T_ERROR);
    }
    tp->t_fmin = fmin;
    tp->t_fmax = fmax;
    tp->t_head = NULL;
    tp->t_anchor = &tp->t_head;
    result = T_PASS;

out:
    return result;
}

/*
 * test_free: free dynamic memory in the test_t struct
 *   @tp: test common state
 */
static void test_free(test_t *tp)
{
    while (tp->t_head != NULL) {
	test_standard_t *tsp = tp->t_head;

	tp->t_head = tsp->ts_next;
	test_standard_free(tsp);
    }
    free((void *)tp->t_parameter_matrix);
    free((void *)tp->t_port_vector);
    vnacal_free(tp->t_vcp);
}

/*
 * test_get_port: randomly choose an unsued port of the parameter matrix
 *   @tp: test common state
 */
static int test_get_port(test_t *tp)
{
    int src_index, dest_index, port;

    assert(tp->t_remaining > 0);
    if (tp->t_remaining > 1) {
	src_index = random() % tp->t_remaining;
    } else {
	src_index = 0;
    }
    dest_index = --tp->t_remaining;	/* index of last unused entry */

    /*
     * Get chosen port number.  Swap with last unused entry.
     */
    port = tp->t_port_vector[src_index];
    if (src_index != dest_index) {
	tp->t_port_vector[src_index] = tp->t_port_vector[dest_index];
	tp->t_port_vector[dest_index] = port;
    }
    return port;
}

/*
 * test_add_traditional_standard: add a traditional parameter standard
 *   @tp: test common state
 */
static libt_result_t test_add_traditional_standard(test_t *tp,
	int row, int column)
{
    int cell = row * tp->t_columns + column;
    test_standard_t *tsp = test_standard_alloc(tp, T_TRADITIONAL_SCALAR);
    int parameter, n;
    double *frequency_vector = NULL;
    double complex *gamma_vector = NULL;
    libt_result_t result = T_FAIL;

    tsp->ts_ports = 1;
    tsp->ts_cell = cell;
    switch (random() % 5) {
    case 0: /* short */
	tsp->ts_scalar = -1.0;
	tp->t_parameter_matrix[cell] = VNACAL_SHORT;
	break;
    case 1: /* open */
	tsp->ts_scalar =  1.0;
	tp->t_parameter_matrix[cell] = VNACAL_OPEN;
	break;
    case 2: /* match */
	tsp->ts_scalar =  0.0;
	tp->t_parameter_matrix[cell] = VNACAL_MATCH;
	break;
    case 3: /* arbitrary scalar */
	tsp->ts_scalar = libt_crandn();
	parameter = vnacal_make_scalar_parameter(tp->t_vcp, tsp->ts_scalar);
	if (parameter == -1) {
	    (void)fprintf(stderr, "%s: vnacal_make_scalar_parameter: %s\n",
		    progname, strerror(errno));
	    return T_FAIL;
	}
	tp->t_parameter_matrix[cell] = parameter;
	break;
    case 4: /* vector */
	tsp->ts_type = T_TRADITIONAL_VECTOR;
	tf2_init(&tsp->ts_vector, tp->t_fmin, tp->t_fmax);
	n = MIN_FPOINTS + random() % (MAX_FPOINTS - MIN_FPOINTS);
	if ((frequency_vector = calloc(n, sizeof(double))) == NULL) {
	    (void)fprintf(stderr, "%s: malloc: %s\n",
		    progname, strerror(errno));
	    exit(T_ERROR);
	}
	if ((gamma_vector = calloc(n, sizeof(double complex))) == NULL) {
	    (void)fprintf(stderr, "%s: malloc: %s\n",
		    progname, strerror(errno));
	    exit(T_ERROR);
	}
	for (int i = 0; i < n; ++i) {
	    double f = tp->t_fmin + (tp->t_fmax - tp->t_fmin) *
		(i / (double)(n - 1));

	    frequency_vector[i] = f;
	    gamma_vector[i] = tf2_eval(&tsp->ts_vector, f);
	}
	parameter = vnacal_make_vector_parameter(tp->t_vcp,
		frequency_vector, n, gamma_vector);
	if (parameter == -1) {
	    (void)fprintf(stderr, "%s: vnacal_make_vector_parameter: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	tp->t_parameter_matrix[cell] = parameter;
	break;
    }
    if (opt_v > 1) {
	test_standard_print(tsp);
    }
    *tp->t_anchor = tsp;
    tp->t_anchor = &tsp->ts_next;
    result = T_PASS;
out:
    free((void *)gamma_vector);
    free((void *)frequency_vector);
    if (result != T_PASS) {
	free((void *)tsp);
    }
    return result;
}

#define MIN_DELAY_CYCLES	(10.0 / 360.0)
#define MAX_DELAY_CYCLES	(5)
#define MIN_RESONANT_FACTOR	(1.5)
#define MAX_RESONANT_FACTOR	(10.0)

/*
 * test_add_calkit_standard: add a calkit standard
 *   @tp: test common state
 */
static libt_result_t test_add_calkit_standard(test_t *tp)
{
    double fmin = tp->t_fmin;
    double fmax = tp->t_fmax;
    double glr, gli, zcr, rt;
    double l, c;
    double k;
    test_standard_t *tsp;
    vnacal_calkit_data_t *vcdp;
    int matrix[2][2];
    int ports;
    libt_result_t result = T_FAIL;

    /*
     * Allocate standard.  Set to 1 port for now -- we will change
     * it if we use through below.
     */
    tsp = test_standard_alloc(tp, T_CALKIT);
    vcdp = &tsp->ts_calkit;
    tsp->ts_ports = ports = 1;

    /*
     * Work backwards choosing values for the real and imaginary parts
     * of gl (glr and gli, respectively) and the real part of zc (zcr)
     * at fmax.  From those, calculate the three offset parameters.
     *
     * First, choose a random imaginary part of gl that delays by at
     * most MAX_DELAY_CYCLES at fmax.  Then choose a random real part
     * that keeps the quantity in the square below (rt) root positive.
     */
    gli = 2.0 * M_PI * libt_randu(MIN_DELAY_CYCLES, MAX_DELAY_CYCLES);
    glr = libt_randu(0.0, (M_SQRT2 - 1.0) * 0.99 * gli);

    /*
     * Choose a random real part of zc between 5 and 500 ohms with mean
     * close to 50 ohms.
     */
    zcr = 50.0 * libt_rand_nsmm(0.832557, 0.5, 0.1, 10.0);

    /*
     * Convert to the offset parameters.
     */
    rt = sqrt(gli*gli - 2.0*glr*gli - glr*glr);
    vcdp->vcd_offset_delay = rt / (2.0 * M_PI * fmax);
    vcdp->vcd_offset_loss = 4.0 * M_PI * glr * zcr * sqrt(1.0e+9 * fmax) / rt;
    vcdp->vcd_offset_z0 = rt * zcr / gli;

    /*
     * Fill in the frequency range.
     */
    vcdp->vcd_fmin = FMIN;
    vcdp->vcd_fmax = FMAX;

    /*
     * Decide between the traditional or revised model.
     */
    if (random() & 1) {
	vcdp->vcd_flags |= VNACAL_CKF_TRADITIONAL;
    }

    /*
     * Find the subtype.
     */
    switch (random() % (tp->t_remaining > 1 ? 4 : 3)) {
    case 0: /* short */
	/*
	 * Calkit short.  Model the inductance of the short as an inductor
	 * in parallel with a small capacitor.	Make the inductor
	 * have 50 ohm reactance in the vicinity of center frequency.
	 * Select the capacitance so the resonant frequency is beyond
	 * fmax by a random amount.
	 */
	vcdp->vcd_type = VNACAL_CALKIT_SHORT;
	l = 50.0 / (2.0 * M_PI * libt_randu(fmin, 2.0 * fmax));
	k = libt_randu(MIN_RESONANT_FACTOR, MAX_RESONANT_FACTOR);
	c = 1.0 / (4.0 * M_PI * M_PI * fmax * fmax * k * k * l);
	compute_l_coefficients(fmin, fmax, l, c, vcdp->vcd_l_coefficients);
	break;

    case 1: /* open */
	/*
	 * Calkit open.  Model the capacitance as a capacitor in series
	 * with a small inductor.  Make the capacitor have 50 ohm reactace
	 * in the vicinity of center frequency.  Select the inductor so
	 * that the resonant frequency is beyond fmax by a random amount.
	 * Because of duality, we can reverse L and C and parallel and
	 * series and use the same function to find the coefficients as
	 * in the short case.
	 */
	vcdp->vcd_type = VNACAL_CALKIT_OPEN;
	c = 1.0 / (50.0 * 2.0 * M_PI * libt_randu(fmin, 2.0 * fmax));
	k = libt_randu(MIN_RESONANT_FACTOR, MAX_RESONANT_FACTOR);
	l = 1.0 / (4.0 * M_PI * M_PI * fmax * fmax * k * k * c);
	compute_l_coefficients(fmin, fmax, c, l, vcdp->vcd_c_coefficients);
	break;

    case 2: /* load */
	/*
	 * Calkit load.  Choose a random complex impedance
	 * anywhere in the complex plane with scale 50 ohms.
	 */
	vcdp->vcd_type = VNACAL_CALKIT_LOAD;
	vcdp->vcd_zl = 50.0 * libt_crandn();
	break;

    case 3: /* through */
	vcdp->vcd_type = VNACAL_CALKIT_THROUGH;
	tsp->ts_ports = ports = 2;
	break;
    }

    /*
     * Make the port map.
     */
    if ((tsp->ts_port_map = malloc(ports * sizeof(int))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n",
		progname, strerror(errno));
	exit(T_ERROR);
    }
    for (int i = 0; i < ports; ++i) {
	tsp->ts_port_map[i] = test_get_port(tp);
    }

    /*
     * Create the standard's parameter matrix and copy it into
     * the main parameter matrix.
     */
    if (ports == 1) {
	int parameter;

	parameter = vnacal_make_calkit_parameter(tp->t_vcp, vcdp);
	if (parameter == -1) {
	    (void)fprintf(stderr, "%s: vnacal_make_calkit_parameter: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	matrix[0][0] = parameter;
    } else {
	if (vnacal_make_calkit_parameter_matrix(tp->t_vcp, vcdp,
		    &matrix[0][0], sizeof(matrix)) == -1) {
	    (void)fprintf(stderr,
		    "%s: vnacal_make_calkit_parameter_matrix: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
    }
    for (int r = 0; r < ports; ++r) {
	int row = tsp->ts_port_map[r];

	for (int c = 0; c < ports; ++c) {
	    int column = tsp->ts_port_map[c];

	    if (row < tp->t_rows && column < tp->t_columns) {
		int cell = row * tp->t_columns + column;
		tp->t_parameter_matrix[cell] = matrix[r][c];
	    } else {
		vnacal_delete_parameter(tp->t_vcp, matrix[r][c]);
	    }
	}
    }
    if (opt_v > 1) {
	test_standard_print(tsp);
    }
    *tp->t_anchor = tsp;
    tp->t_anchor = &tsp->ts_next;
    result = T_PASS;

out:
    return result;
}

/*
 * test_add_data_standard: add a data standard
 *   @tp: test common state
 */
static libt_result_t test_add_data_standard(test_t *tp)
{
    int ports;
    int frequencies;;
    bool use_fz0;
    const double fmin = tp->t_fmin;
    const double fmax = tp->t_fmax;
    const double fc = (fmin + fmax) / 2.0;
    test_standard_t *tsp = NULL;
    vnadata_t *vdp = NULL;
    int *parameter_matrix = NULL;
    libt_result_t result = T_FAIL;

    /*
     * Allocate standard.
     */
    ports = random() % MIN(4, tp->t_remaining) + 1;
    tsp = test_standard_alloc(tp, T_DATA);
    tsp->ts_ports = ports;

    /*
     * Make the port map.
     */
    if ((tsp->ts_port_map = malloc(ports * sizeof(int))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n",
		progname, strerror(errno));
	exit(T_ERROR);
    }
    for (int i = 0; i < tsp->ts_ports; ++i) {
	tsp->ts_port_map[i] = test_get_port(tp);
    }

    /*
     * Create the z0 vector.  When using fz0, this is the impedance
     * at center frequency.  Our values have an average magnitude of
     * 50 ohms with a real part that has a magnitude of at least 5
     * ohms.
     */
    use_fz0 = random() & 1;	  /* use frequency-dependent z0? */
    tsp->ts_has_fz0 = use_fz0;
    if ((tsp->ts_z0_vector = calloc(ports, sizeof(double complex))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n",
		progname, strerror(errno));
	exit(T_ERROR);
    }
    for (int i = 0; i < ports; ++i) {
	tsp->ts_z0_vector[i] = 50.0 * libt_crand_nsmmra(0.832557, 0.5,
		0.3, 1000.0, 0.0, -140.0);
    }

    /*
     * Create transfer functions for the data elements.
     */
    if ((tsp->ts_matrix = calloc(ports * ports, sizeof(tf2_t))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n",
		progname, strerror(errno));
	exit(T_ERROR);
    }
    for (int i = 0; i < ports * ports; ++i) {
	tf2_init(&tsp->ts_matrix[i], fmin, fmax);
    }

    /*
     * Allocate and fill the data structure.
     */
    if ((vdp = vnadata_alloc(error_fn, NULL)) == NULL) {
	(void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
		progname, strerror(errno));
	result = T_FAIL;
	goto out;
    }
    frequencies = MIN_FPOINTS + random() % (MAX_FPOINTS - MIN_FPOINTS);
    if (vnadata_init(vdp, VPT_S, ports, ports, frequencies) == -1) {
	result = T_FAIL;
	goto out;
    }
    if (!use_fz0) {
	if (vnadata_set_z0_vector(vdp, tsp->ts_z0_vector) == -1) {
	    result = T_FAIL;
	    goto out;
	}
    }
    for (int findex = 0; findex < frequencies; ++findex) {
	double f = fmin + (fmax - fmin) * findex / (double)(frequencies - 1);
	double complex *data_matrix = vnadata_get_matrix(vdp, findex);

	if (use_fz0) {
	    for (int port = 0; port < ports; ++port) {
		double complex z0 = eval_fz0(tsp->ts_z0_vector[port], fc, f);

		if (vnadata_set_fz0(vdp, findex, port, z0) == -1) {
		    result = T_FAIL;
		    goto out;
		}
	    }
	}
	if (vnadata_set_frequency(vdp, findex, f) == -1) {
	    result = T_FAIL;
	    goto out;
	}
	for (int i = 0; i < ports; ++i) {
	    for (int j = 0; j < ports; ++j) {
		int cell = ports * i + j;

		data_matrix[cell] = tf2_eval(&tsp->ts_matrix[cell], f);
	    }
	}
    }

    /*
     * Create the local parameter matrix and copy parameters into
     * the common parameter matrix.
     */
    if ((parameter_matrix = malloc(ports * ports * sizeof(int))) == NULL) {
	(void)fprintf(stderr, "%s: malloc: %s\n", progname, strerror(errno));
	exit(T_ERROR);
    }
    for (int cell = 0; cell < ports * ports; ++cell) {
	parameter_matrix[cell] = -1;
    }
    if (vnacal_make_data_parameter_matrix(tp->t_vcp, vdp,
		parameter_matrix, ports * ports * sizeof(int)) == -1) {
	result = T_FAIL;
	goto out;
    }
    for (int r = 0; r < ports; ++r) {
	int row = tsp->ts_port_map[r];

	for (int c = 0; c < ports; ++c) {
	    int standard_cell = r * ports + c;
	    int column = tsp->ts_port_map[c];

	    if (row < tp->t_rows && column < tp->t_columns) {
		int cell = row * tp->t_columns + column;
		tp->t_parameter_matrix[cell] =
		    parameter_matrix[standard_cell];
	    } else {
		vnacal_delete_parameter(tp->t_vcp,
			parameter_matrix[standard_cell]);
	    }
	}
    }
    free((void *)parameter_matrix);
    if (opt_v > 1) {
	test_standard_print(tsp);
    }
    *tp->t_anchor = tsp;
    tp->t_anchor = &tsp->ts_next;
    result = T_PASS;

out:
    vnadata_free(vdp);
    if (result != T_PASS) {
	test_standard_free(tsp);
    }
    return result;
}

/*
 * test_add_standards: add random standards
 *   @tp: test common state
 */
static libt_result_t test_add_standards(test_t *tp)
{
    libt_result_t result = T_FAIL;

    while (tp->t_remaining > 0) {
	/*
	 * Choose between traditional parameters, a calkit standard
	 * or a data standard.
	 */
	switch (random() % 3) {
	case 0: /* traditional parameter */
	    /*
	     * Just reserve the port now.  We'll add standards for
	     * them below.
	     */
	    {
		int port = test_get_port(tp);

		tp->t_traditional |= 1U << port;
	    }
	    break;

	case 1: /* calkit standard */
	    result = test_add_calkit_standard(tp);
	    if (result != T_PASS)
		return result;
	    break;

	case 2: /* data standard */
	    result = test_add_data_standard(tp);
	    if (result != T_PASS)
		return result;
	    break;
	}
    }
    /*
     * At every intersection of ports reserved as traditional, add
     * a traditional standard.
     */
    for (int row = 0; row < tp->t_rows; ++row) {
	if (!(tp->t_traditional & 1U << row)) {
	    continue;
	}
	for (int column = 0; column < tp->t_columns; ++column) {
	    if (!(tp->t_traditional & 1U << column)) {
		continue;
	    }
	    result = test_add_traditional_standard(tp, row, column);
	    if (result != T_PASS)
		return result;
	}
    }
    return T_PASS;
}

/*
 * eval_standards: evaluate the standards at frequency f
 *   @tp: test common state
 *   @result_matrix: caller-allocated matrix to receive values
 *   @f: frequency at which to evaluate
 */
static void eval_standards(test_t *tp, double f,
	const double complex *z0_vector, double complex *result_matrix)
{
    int rows = tp->t_rows;
    int columns = tp->t_columns;

#if BINARY_ZERO_IS_DOUBLE_ZERO
    (void)memset((void *)result_matrix, 0,
	    rows * columns * sizeof(double complex));
#else
    for (int i = 0; i < rows * columns; ++i) {
	result[i] = 0.0;
    }
#endif
    for (test_standard_t *tsp = tp->t_head; tsp != NULL; tsp = tsp->ts_next) {
	vnacal_calkit_data_t *vcdp = NULL;
	double complex temp_z1[tsp->ts_ports]; /* standard's z0 */
	double complex temp_z2[tsp->ts_ports]; /* requested z0 */
	double complex temp_matrix[tsp->ts_ports][tsp->ts_ports];
	switch (tsp->ts_type) {
	case T_TRADITIONAL_SCALAR:
	    result_matrix[tsp->ts_cell] = tsp->ts_scalar;
	    continue;  /* skip renormalization below */

	case T_TRADITIONAL_VECTOR:
	    result_matrix[tsp->ts_cell] = tf2_eval(&tsp->ts_vector, f);
	    continue;  /* skip renormalization below */

	case T_CALKIT:
	    /*
	     * The keysight functions above work with only positive real
	     * z0, and the through requires the same z0 on both ports.
	     * To support general z0 values, evaluate all with a constant
	     * 50 ohm reference, then renormalize.
	     */
	    for (int i = 0; i < tsp->ts_ports; ++i) {
		temp_z1[i] = 50.0;
	    }
	    (void)memset((void *)temp_matrix, 0, sizeof(temp_matrix));
	    vcdp = &tsp->ts_calkit;
	    switch (vcdp->vcd_type) {
	    case VNACAL_CALKIT_SHORT:
		assert(tsp->ts_ports == 1);
		temp_matrix[0][0] = keysight_short(vcdp, f, 50.0);
		break;

	    case VNACAL_CALKIT_OPEN:
		assert(tsp->ts_ports == 1);
		temp_matrix[0][0] = keysight_open(vcdp, f, 50.0);
		break;

	    case VNACAL_CALKIT_LOAD:
		assert(tsp->ts_ports == 1);
		temp_matrix[0][0] = keysight_load(vcdp, f, 50.0);
		break;

	    case VNACAL_CALKIT_THROUGH:
		assert(tsp->ts_ports == 2);
		keysight_through(vcdp, f, 50.0, temp_matrix);
		break;

	    default:
		abort();
	    }
	    break;

	case T_DATA:
	    if (tsp->ts_has_fz0) {
		double fmin = tp->t_fmin;
		double fmax = tp->t_fmax;
		double fc = (fmin + fmax) / 2.0;

		for (int i = 0; i < tsp->ts_ports; ++i) {
		    temp_z1[i] = eval_fz0(tsp->ts_z0_vector[i], fc, f);
		}
	    } else {
		(void)memcpy((void *)temp_z1, (void *)tsp->ts_z0_vector,
			tsp->ts_ports * sizeof(double complex));
	    }
	    for (int r = 0; r < tsp->ts_ports; ++r) {
		for (int c = 0; c < tsp->ts_ports; ++c) {
		    int cell = tsp->ts_ports * r + c;

		    temp_matrix[r][c] = tf2_eval(&tsp->ts_matrix[cell], f);
		}
	    }
	    break;

	default:
	    abort();
	}

	/*
	 * Renormalize to the requested z0 vector.
	 */
	for (int i = 0; i < tsp->ts_ports; ++i) {
	    temp_z2[i] = z0_vector[tsp->ts_port_map[i]];
	}
	vnaconv_stosrn(*temp_matrix, *temp_matrix,
		temp_z1, temp_z2, tsp->ts_ports);

	/*
	 * Copy result.
	 */
	for (int r = 0; r < tsp->ts_ports; ++r) {
	    int row = tsp->ts_port_map[r];

	    for (int c = 0; c < tsp->ts_ports; ++c) {
		int column = tsp->ts_port_map[c];

		if (row < tp->t_rows && column < tp->t_columns) {
		    int cell = row * tp->t_columns + column;

		    result_matrix[cell] = temp_matrix[r][c];
		}
	    }
	}
    }
}

/*
 * check_result: check the actual matrix against expected at freuqency
 *   @tp: test common struct
 *   @findex: frequency index
 *   @frequency: frequency at which to evaluate
 *   @z0_vector: reference impedances
 *   @actual_matrix: actual values to compare against expected
 */
static libt_result_t check_result(test_t *tp, int findex, double frequency,
	double complex *z0_vector, const double complex *actual_matrix)
{
    int rows = tp->t_rows;
    int columns = tp->t_columns;
    double complex expected_matrix[rows][columns];
    libt_result_t result = T_FAIL;

    if (opt_v > 1) {
	(void)printf("findex %d frequency %e\n", findex, frequency);
	libt_print_cmatrix("actual", actual_matrix, rows, columns);
    }
    eval_standards(tp, frequency, z0_vector, &expected_matrix[0][0]);
    if (opt_v > 1) {
	libt_print_cmatrix("expected", &expected_matrix[0][0],
		rows, columns);
    }
    for (int row = 0; row < rows; ++row) {
	for (int column = 0; column < columns; ++column) {
	    int cell = row * columns + column;
	    char label[1 + 2 * 3 * sizeof(int) + 1];

	    (void)sprintf(label, "s%d%d", row + 1, column + 1);
	    TEST_EQUAL(actual_matrix[cell],
		    expected_matrix[row][column], label);
	}
    }

    /*
     * For single-ported standards, also check vnacal_eval_parameter.
     */
    if (opt_v > 1) {
	(void)printf("test vnacal_eval_parameter:\n");
    }
    for (test_standard_t *tsp = tp->t_head; tsp != NULL; tsp = tsp->ts_next) {
	int r, c, cell;
	double complex z0;
	double complex value;
	char label[1 + 2 * 3 * sizeof(int) + 1];

	if (tsp->ts_ports != 1) {
	    continue;
	}
	if (tsp->ts_type == T_TRADITIONAL_SCALAR ||
	    tsp->ts_type == T_TRADITIONAL_VECTOR) {
	    cell = tsp->ts_cell;
	    r = cell / columns;
	    c = cell % columns;
	    z0 = 1.0;  /* ignored */
	} else {
	    int port = tsp->ts_port_map[0];

	    r = port;
	    c = port;
	    cell = (columns + 1) * port;
	    z0 = z0_vector[port];
	}
	if (r >= rows || c >= columns) {
	    continue;
	}
	value = vnacal_eval_parameter(tp->t_vcp,
		tp->t_parameter_matrix[cell], frequency, z0);
	(void)sprintf(label, "s%d%d", cell / columns + 1, cell % columns + 1);
        TEST_EQUAL(value, expected_matrix[r][c], label);
    }
    if (opt_v > 1) {
	(void)printf("\n");
    }
    result = T_PASS;

out:
    return result;
}

/*
 * run_trial: run a test trial
 *   @trial: trial number
 *   @rows: number of rows in parameter matrix
 *   @columns: number of columns in parameter matrix
 */
static libt_result_t run_trial(int trial, int rows, int columns)
{
    test_t t;
    int ports = MAX(rows, columns);
    int frequencies;
    double complex z0_vector[ports];
    vnadata_t *vdp = NULL;
    libt_result_t result = T_FAIL;

    /*
     * If -v, print the test header.
     */
    if (opt_v != 0) {
	(void)printf("Test vnacal parameter matrix: trial %3d size %d x %d\n",
	    trial, rows, columns);
    }

    /*
     * Init the test common arguments.
     */
    if ((result = test_init(&t, rows, columns, FMIN, FMAX)) != T_PASS)
	goto out;
    if ((result = test_add_standards(&t)) != T_PASS)
	goto out;

    /*
     * Choose target reference impedances.
     */
    for (int i = 0; i < ports; ++i) {
	z0_vector[i] = 50.0 * libt_crand_nsmmra(0.832557, 0.5,
		0.3, 1000.0, 0.0, -140.0);
    }

    /*
     * Choose frequency points, evaulate at each frequency and compare.
     * On odd trials, use vnacal_eval_parameter_matrix; on eval trials
     * use vnacal_parameter_matrix_to_data.
     */
    frequencies = MIN_FPOINTS + random() % (MAX_FPOINTS - MIN_FPOINTS);
    if (trial & 1) {
	for (int findex = 0; findex < frequencies; ++findex) {
	    double f = FMIN + (FMAX - FMIN) * findex /
	               (double)(frequencies - 1);
	    double complex actual_matrix[rows][columns];

	    for (int r = 0; r < rows; ++r) {
		for (int c = 0; c < columns; ++c) {
		    actual_matrix[r][c] = NAN;
		}
	    }
	    if (vnacal_eval_parameter_matrix(t.t_vcp, t.t_parameter_matrix,
			rows, columns, f, z0_vector,
			&actual_matrix[0][0]) == -1) {
		(void)fprintf(stderr, "%s: vnacal_eval_parameter_matrix: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    result = check_result(&t, findex, f, z0_vector,
		    &actual_matrix[0][0]);
	    if (result != T_PASS)
	        goto out;
	}
    } else {
	vnadata_parameter_type_t type = (rows == columns) ? VPT_S : VPT_UNDEF;

	if ((vdp = vnadata_alloc_and_init(error_fn, (void *)NULL, type,
			rows, columns, frequencies)) == NULL) {
	    (void)fprintf(stderr, "%s: vnadata_alloc: %s\n",
		    progname, strerror(errno));
	    goto out;
	}
	for (int findex = 0; findex < frequencies; ++findex) {
	    double f = FMIN + (FMAX - FMIN) * findex /
	               (double)(frequencies - 1);

	    if (vnadata_set_frequency(vdp, findex, f) == -1) {
		(void)fprintf(stderr, "%s: vnadata_set_frequency: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	}
	if (vnadata_set_z0_vector(vdp, z0_vector) == -1) {
	    (void)fprintf(stderr, "%s: vnadata_set_z0_vector: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	if (vnacal_parameter_matrix_to_data(t.t_vcp, t.t_parameter_matrix,
		    rows, columns, vdp) == -1) {
	    (void)fprintf(stderr, "%s: vnacal_parameter_matrix_to_data: %s\n",
		    progname, strerror(errno));
	    result = T_FAIL;
	    goto out;
	}
	for (int findex = 0; findex < frequencies; ++findex) {
	    double f;
	    double complex *actual_matrix = vnadata_get_matrix(vdp, findex);

	    if ((f = vnadata_get_frequency(vdp, findex)) == HUGE_VAL) {
		(void)fprintf(stderr, "%s: vnadata_get_frequency: %s\n",
			progname, strerror(errno));
		result = T_FAIL;
		goto out;
	    }
	    result = check_result(&t, findex, f, z0_vector, actual_matrix);
	    if (result != T_PASS)
	        goto out;
	}
	vnadata_free(vdp);
	vdp = NULL;
    }
    result = T_PASS;

out:
    vnadata_free(vdp);
    test_free(&t);
    return result;
}

/*
 * run_test: run the test
 */
static libt_result_t run_test()
{
    libt_result_t result = T_FAIL;

    for (int trial = 1; trial <= NTRIALS; ++trial) {
	if (trial <= NRECTANGULAR) {
	    for (int rows = 0; rows <= 7; ++rows) {
		for (int columns = 0; columns <= 7; ++columns) {
		    if ((result = run_trial(trial, rows, columns)) != T_PASS) {
			goto out;
		    }
		}
	    }
	} else {
	    for (int ports = 1; ports <= 7; ++ports) {
		if ((result = run_trial(trial, ports, ports)) != T_PASS) {
		    goto out;
		}
	    }
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
    exit(99);
}

/*
 * main
 */
int
main(int argc, char **argv)
{
    /*
     * Parse Options
     */
    if ((char *)NULL == (progname = strrchr(argv[0], '/'))) {
	progname = argv[0];
    } else {
	++progname;
    }
    for (;;) {
	switch (getopt(argc, argv, options)) {
	case -1:
	    break;

	case 'a':
	    opt_a = true;
	    continue;

	case 'v':
	    ++opt_v;
	    continue;

	default:
	    print_usage();
	}
	break;
    }
    libt_isequal_init();
    exit(run_test());
}
