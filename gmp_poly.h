#ifndef GMP_POLY
#define GMP_POLY


#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


// Index of coeff equals its power. x[0] = x ^ 0
typedef struct Fq_poly
{
	int degree;
	mpz_t *coeffs;
} Fq_poly;

mpz_t q;


#include "gmp_poly.c"
#include "gmp_utils.c"
#include "gmp_poly_testing.c"

#include "lagrangeInterpolation.h"

#endif