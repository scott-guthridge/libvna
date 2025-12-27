#ifndef _R_H
#define _R_H

#include <complex.h>

/*
 * LIBT_IRLOG4: 1.0 / sqrt(log(4.0))
 *   The Rayleigh distribution has a median of 1 with this sigma value.
 */
#define LIBT_IRLOG4	0.84932180028801904272

/*
 * libt_crand_generator_t: abstract complex random number generator
 */
typedef struct libt_crand_generator libt_crand_generator_t;
typedef double complex libt_crand_function_t(libt_crand_generator_t *cgp);
struct libt_crand_generator {
    libt_crand_function_t *cg_crand;
};

/* libt_crandn: return standard complex normal random numbers */
extern double complex libt_crandn();

/* libt_crandn_s: return complex normal random numbers with scale factor */
extern double complex libt_crandn_s(double sigma);

/* libt_crand_nsmm: return complex random numbers with magnitude distributed
   according to a truncated Rice(nu, sigma) distribution */
extern double complex libt_crand_nsmm(double nu, double sigma,
	double min, double max);

/* libt_crand_nsmmra: complex random numbers with truncated mag and ang */
extern double complex libt_crand_nsmmra(double nu, double sigma,
	double min, double max, double rotation, double angle);

/* libt_crand_generator: return a random generator for the given parameters */
extern libt_crand_generator_t *libt_crand_generator(double nu, double sigma,
	double min, double max, double rotation, double angle);

#endif /* _R_H */
