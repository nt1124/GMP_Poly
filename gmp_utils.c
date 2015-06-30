#include <ctime>

gmp_randstate_t *seedRandGen()
{
	gmp_randstate_t *state = (gmp_randstate_t*) calloc(1, sizeof(gmp_randstate_t));
	unsigned long int seed = time(NULL);
	gmp_randinit_default(*state);
	gmp_randseed_ui(*state, seed);

	return state;
}
