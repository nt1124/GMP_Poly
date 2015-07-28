#ifndef GMP_POLY
#define GMP_POLY


// Ronseal
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


// Index of coeff equals its power. x[0] = coeffs of x ^ 0 factor
// Therefore an Fq_Poly of degree N has N + 1 Coeffs
typedef struct Fq_poly
{
	int degree;
	mpz_t *coeffs;
} Fq_poly;


typedef struct PointwiseRep
{
	int numPoints;
	mpz_t *evalPoints;
} PointwiseRep;


#include "gmp_poly.c"
#include "gmp_utils.c"


#include "interpolationUtils.c"
#include "PointwiseRepUtils.c"


#include "Testing/gmp_poly_testing.c"


#endif