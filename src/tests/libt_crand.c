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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libt.h"
#include "libt_crand.h"


/*
 * rice_cdf_inner: inner summation
 */
static inline double rice_cdf_inner(int m, double u)
{
    double kf = 1.0;
    double s  = 1.0;
    double up = 1.0;

    if (u != 0.0) {
	for (int k = 1; k < m; ++k) {
	    double ds;

	    up *= u;
	    kf *= (double)k;
	    ds = up / kf;
	    if (!isnormal(ds)) {
		break;
	    }
	    s += ds;
	}
    }
    return 1.0 - exp(-u) * s;
}

/*
 * rice_cdf: compute the CDF of the Rice(ν, σ) distribution
 *   @nu: distance from the reference point to the origin
 *   @sigma: scale factor
 *
 *   We calculate the Marcum-Q function using the algoritm described
 *   by D.A. Shnidman (1989).  "The Calculation of the Probability of
 *   Detection and the Generalized Marcum Q-Function." IEEE Transactions
 *   on Information Theory, 35(2), 389-400.
 *
 * Return:
 *   CDF in [0.0 .. 1.0]
 */
static double rice_cdf(double nu, double sigma, double x)
{
    double ss;
    double u, v;
    double s;
    double mf;
    double vp;

    assert(nu >= 0.0 && sigma >= 0.0 && x >= 0.0);
    if (sigma == 0.0) {
	return x < nu ? 0.0 : x == nu ? 0.5 : 1.0;
    }
    ss = 2.0 * sigma * sigma;
    u  = nu * nu / ss;
    v  = x  * x  / ss;
    s  = 1.0;
    mf = 1.0;
    vp = 1.0;
    for (int m = 1; /*EMPTY*/; ++m) {
	double ds;

	vp *= v;
	mf *= (double)m;
	ds = rice_cdf_inner(m, u) * vp / mf;
	if (!isnormal(ds)) {
	    break;
	}
	s += ds;
	if (ds < 1.0e-6) {
	    break;
	}
    }
    return 1.0 - exp(-v) * s;
}

/*
 * rice_inverse_cdf: return the inverse of the Rice CDF
 *   @nu: distance from the reference point to the origin
 *   @sigma: scale factor
 *   @q: quantile [0.0 .. 1.0]
 */
static double rice_inverse_cdf(double nu, double sigma, double q)
{
    double a = 0.0;
    double b = 1.0;
    double fa = rice_cdf(nu, sigma, a) - q;
    double fb = rice_cdf(nu, sigma, b) - q;
    int side = 0;

    assert(0.0 <= nu && 0.0 <= sigma && 0 <= q && q <= 1.0);
    if (q == 1.0) {
#ifdef INFINITY
	return INFINITY;
#else
	return 1.0e+37;
#endif
    }

    /*
     * Surround the solution.
     */
    while (fb < 0.0) {
	assert(fa < 0.0);
	a = b;
	fa = fb;
	b *= 2.0;
	fb = rice_cdf(nu, sigma, b) - q;
    }

    /*
     * Use simplifed false position with Illinois method to find the root.
     */
    for (int k = 0; k < 50; ++k) {
	double c = (fa * b - fb * a) / (fa - fb);
	double fc;

	assert(fa <  0.0);
	assert(fb >= 0.0);
	if (c <= a || c >= b) {	/* If out of bounds, use bisection. */
	    c = (a + b) / 2.0;
	}
	if (fabs(b - a) < 1.0e-6) {
	    return c;
	}
	fc = rice_cdf(nu, sigma, c) - q;
	if (fc < 0.0) {
	    a = c;
	    fa = fc;
	    if (side == +1) {
		fb /= 2.0;
	    } else {
		side = +1;
	    }
	} else if (fc > 0.0) {
	    b = c;
	    fb = fc;
	    if (side == -1) {
		fa /= 2.0;
	    } else {
		side = -1;
	    }
	} else {
	    return c;
	}
    }
    errno = EDOM;
    return HUGE_VAL;
}

/*
 * libt_rand_nsmm: return a truncated Rice(nu, sigma) random number
 *   @nu:    distance from reference point to origin
 *   @sigma: scale
 *   @min:   minimum value to return
 *   @max:   maximum value to return
 */
double libt_rand_nsmm(double nu, double sigma, double min, double max)
{
    assert(0.0 <= nu && 0.0 <= sigma && 0.0 <= min && min <= max);

    /*
     * First try generating a few rice-distributed random numbers,
     * returning the first that we find in bounds.
     */
    for (int i = 0; i < 4; ++i) {
	double complex z = nu + libt_crandn_s(sigma);
	double r = cabs(z);

	if (r >= min && r <= max) {
	    return r;
	}
    }

    /*
     * If that failed, use the general approach of finding the quantiles
     * of the min and max, generating a quantile equally distributed
     * between those limits then finding the inverse Rice distribution.
     */
    {
	double q1 = rice_cdf(nu, sigma, min);
	double q2 = rice_cdf(nu, sigma, max);
	double u1 = (double)random() / RANDOM_MAX;
	double q = q1 + (q2 - q1) * u1;
	double r = rice_inverse_cdf(nu, sigma, q);

	/*
	 * r should already be in range, but if due to round-off error,
	 * it's slightly out of bounds, move it back in.
	 */
	if (r < min) {
	    r = min;
	} else if (r > max) {
	    r = max;
	}
	return r;
    }
}

/*
 * libt_crandn: return standard complex normal random numbers
 *     mode   of magnitude: sqrt(2)/2      ~~= 0.7071067811865475
 *     median of magnitude: sqrt(log(4)/2) ~~= 0.8325546111576978
 *     mean   of magnitude: sqrt(pi)/2     ~~= 0.8862269254527580
 */
double complex libt_crandn()
{
    double u1 = (random() + 1.0) / (RANDOM_MAX + 1.0); /* Box Muller method */
    double u2 = (double)random() / (RANDOM_MAX + 1.0);
    double r = sqrt(-log(u1));
    double a = 2 * M_PI * u2;

    return r * (cos(a) + I * sin(a));
}

/*
 * libt_crandn_s: return complex normal random numbers with scale factor
 *     mode   of magnitude: sigma
 *     median of magnitude: sigma * sqrt(2*log(2))
 *     mean   of magnitude: sigma * sqrt(pi/2)
 */
double complex libt_crandn_s(double sigma)
{
    double u1, u2, r, a;

    assert(sigma >= 0.0);
    u1 = (random() + 1.0) / (RANDOM_MAX + 1.0);    /* Box Muller method */
    u2 = (double)random() / (RANDOM_MAX + 1.0);
    r = sqrt(-2.0 * log(u1)) * sigma;
    a = 2 * M_PI * u2;
    return r * (cos(a) + I * sin(a));
}

/*
 * libt_crand_nsmm: return complex random numbers with magnitude
 *     distributed according to a truncated Rice(nu, sigma) distribution
 *   @nu:    distance from reference point to origin
 *   @sigma: scale
 *   @min:   minimum magnitude to return
 *   @max:   maximum magnitude to return
 */
double complex libt_crand_nsmm(double nu, double sigma, double min, double max)
{
    double r, u;

    assert(nu >= 0.0 && sigma >= 0.0 && min >= 0 && min <= max);
    r = libt_rand_nsmm(nu, sigma, min, max);
    u = (double)random() / (RANDOM_MAX + 1.0);
    return r * cexp(2 * M_PI * I * u);
}

/*
 * libt_crand_nsmmra: return complex random numbers with truncated mag and ang
 *   @nu:    distance from reference point to origin
 *   @sigma: scale
 *   @min:   minimum magnitude to return
 *   @max:   maximum magnitude to return
 *   @rotation: center of arc
 *   @angle: interior angle of arc (mirror if negative)
 */
double complex libt_crand_nsmmra(double nu, double sigma,
	double min, double max, double rotation, double angle)
{
    double r;

    assert(nu >= 0.0 && sigma >= 0.0 && min >= 0 && min <= max &&
	    -360.0 <= rotation && rotation <= 360.0 &&
	    -360.0 <= angle && angle <= 360.0);
    r = libt_rand_nsmm(nu, sigma, min, max);
    rotation *= M_PI / 180.0;
    if (angle != 0.0) {
	double u = (double)random() / RANDOM_MAX;

	if (angle < 0.0) {
	    u *= 2.0;
	    if (u >= 1.0) {
		u -= 1.0;
		rotation += M_PI;
	    }
	}
	rotation += (u - 0.5) * angle / 180.0 * M_PI;
    }
    return r * cexp(rotation * I);
}

/*
 * _cga: subclass of libt_crand_generator_t for _cga
 */
typedef struct _cga {
    libt_crand_generator_t cga_base;
    double		cga_sigma;
} _cga_t;

/*
 * _cga: complex random generator with sigma
 *   @cgp: instance of _cga_t
 */
static double complex _cga(libt_crand_generator_t *cgp)
{
    _cga_t *cgap = (_cga_t *)cgp;

    return libt_crandn_s(cgap->cga_sigma);
}

/*
 * _cga: subclass of libt_crand_generator_t for _cg1
 */
typedef struct _cg1 {
    libt_crand_generator_t cg_base;
    double		cg1_nu;
    double		cg1_sigma;
    double		cg1_min;
    double             	cg1_max;
    double		cg1_rotation;	/* in radians */
    double		cg1_angle;	/* in radians */
} _cg1_t;

/*
 * _cg1: iterative general purpose truncated magnitude and angle random function
 *   @cgp: an instance of _cg2_t
 */
static double complex _cg1(libt_crand_generator_t *cgp)
{
    _cg1_t *cg1p = (_cg1_t *)cgp;
    double r;
    double rotation = cg1p->cg1_rotation;
    double angle    = cg1p->cg1_angle;

    for (int i = 0; /*EMPTY*/; ++i) {
	double complex z = cg1p->cg1_nu + libt_crandn_s(cg1p->cg1_sigma);

	r = cabs(z);
	if (r >= cg1p->cg1_min && r <= cg1p->cg1_max) {
	    break;
	}
	if (i == 50) {
	    abort();
	}
    }
    if (angle != 0.0) {
	double u = (double)random() / RANDOM_MAX;

	if (angle < 0.0) {
	    u *= 2.0;
	    if (u >= 1.0) {
		u -= 1.0;
		rotation += M_PI;
	    }
	}
	rotation += (u - 0.5) * angle;
    }
    return r * cexp(rotation * I);
}

/*
 * _cga: subclass of libt_crand_generator_t for _cg2
 */
typedef struct _cg2 {
    libt_crand_generator_t cg_base;
    double		cg2_nu;
    double		cg2_sigma;
    double		cg2_min;
    double             	cg2_max;
    double		cg2_q1;
    double		cg2_q2;
    double		cg2_rotation;	/* in radians */
    double		cg2_angle;	/* in radians */
} _cg2_t;

/*
 * _cg2: general purpose truncated magnitude and angle random function
 *   @cgp: an instance of _cg2_t
 */
static double complex _cg2(libt_crand_generator_t *cgp)
{
    _cg2_t *cg2p = (_cg2_t *)cgp;
    double u1 = (double)random() / RANDOM_MAX;
    double q = cg2p->cg2_q1 + (cg2p->cg2_q2 - cg2p->cg2_q1) * u1;
    double r = rice_inverse_cdf(cg2p->cg2_nu, cg2p->cg2_sigma, q);
    double rotation = cg2p->cg2_rotation;
    double angle    = cg2p->cg2_angle;

    /*
     * r should already be in range, but if due to round-off error,
     * it's slightly out of bounds, move it back in.
     */
    if (r < cg2p->cg2_min) {
	r = cg2p->cg2_min;
    } else if (r > cg2p->cg2_max) {
	r = cg2p->cg2_max;
    }

    /*
     * Generate the random angle.
     */
    if (angle != 0.0) {
	double u2 = (double)random() / (RANDOM_MAX + 1.0);

	if (angle < 0.0) {
	    u2 *= 2.0;
	    if (u2 >= 1.0) {
		u2 -= 1.0;
		rotation += M_PI;
	    }
	}
	rotation += (u2 - 0.5) * angle;
    }
    return r * cexp(rotation * I);
}

/*
 * libt_crand_generator: return a random generator for the given parameters
 *   @nu:    distance from reference point to origin
 *   @sigma: scale
 *   @min:   minimum magnitude to return
 *   @max:   maximum magnitude to return
 *   @rotation: center of arc
 *   @angle: interior angle of arc (mirror if negative)
 *
 * This function returns a random number generator that produces complex
 * random numbers with magnitude between between min and max, following
 * a truncated Rice(ν, σ) distribution, and angle between rotation
 * plus or minus angle over two.  If angle is negative, then the range
 * of angles also includes minus rotation plus or minus angle over two.
 *
 * Caller is responsible for freeing the returned structure by a
 * call to free().
 */
libt_crand_generator_t *libt_crand_generator(double nu, double sigma,
	double min, double max, double rotation, double angle)
{
    libt_crand_generator_t *cgp = NULL;

    assert(nu >= 0.0 && sigma >= 0.0 && min >= 0 && min <= max &&
	    -360.0 <= rotation && rotation <= 360.0 &&
	    -360.0 <= angle && angle <= 360.0);

    /*
     * Handle the simple non-truncated cases first.
     */
    if (nu == 0.0 && min == 0.0 && rice_cdf(nu, sigma, max) >= 0.9999 &&
	    angle == 360.0) {
	_cga_t *cgap;

	if (sigma == M_SQRT1_2) {
	    if ((cgp = malloc(sizeof(libt_crand_generator_t))) == NULL) {
		(void)fprintf(stderr, "%s: malloc: %s\n",
			progname, strerror(errno));
		exit(99);
	    }
	    cgp->cg_crand = libt_crandn;

	} else {
	    if ((cgp = malloc(sizeof(_cga_t))) == NULL) {
		(void)fprintf(stderr, "%s: malloc: %s\n",
			progname, strerror(errno));
		exit(99);
	    }
	    cgap = (_cga_t *)cgp;
	    cgap->cga_sigma = sigma;
	    cgp->cg_crand = _cga;
	}

    } else {
	double q1 = rice_cdf(nu, sigma, min);
	double q2 = rice_cdf(nu, sigma, max);

	/*
	 * If there is a reasonable probability that we'll find a solution
	 * with magnitude within min and max in a few iterations, use
	 * the iterative generator.  The iterative generator is quite
	 * a bit faster than the general one, so making a few tries is
	 * faster.
	 */
	if (q2 - q1 >= 0.25) {
	    _cg1_t *cg1p;

	    if ((cgp = malloc(sizeof(_cg1_t))) == NULL) {
		(void)fprintf(stderr, "%s: malloc: %s\n",
			progname, strerror(errno));
		exit(99);
	    }
	    cg1p = (_cg1_t *)cgp;
	    cg1p->cg1_nu       = nu;
	    cg1p->cg1_sigma    = sigma;
	    cg1p->cg1_min      = min;
	    cg1p->cg1_max      = max;
	    cg1p->cg1_rotation = rotation * M_PI / 180.0;
	    cg1p->cg1_angle    = angle    * M_PI / 180.0;
	    cgp->cg_crand = _cg1;

	} else {
	    _cg2_t *cg2p;

	    if ((cgp = malloc(sizeof(_cg2_t))) == NULL) {
		(void)fprintf(stderr, "%s: malloc: %s\n",
			progname, strerror(errno));
		exit(99);
	    }
	    cg2p = (_cg2_t *)cgp;
	    cg2p->cg2_nu       = nu;
	    cg2p->cg2_sigma    = sigma;
	    cg2p->cg2_min      = min;
	    cg2p->cg2_max      = max;
	    cg2p->cg2_q1       = q1;
	    cg2p->cg2_q2       = q2;
	    cg2p->cg2_rotation = rotation * M_PI / 180.0;
	    cg2p->cg2_angle    = angle    * M_PI / 180.0;
	    cgp->cg_crand = _cg2;
	}
    }
    return cgp;
}
