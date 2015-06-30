// Testing Functions.



void test_evaluation()
{
	struct Fq_poly *polyToEval;

	mpz_t *coeffs, x, *y;
	gmp_randstate_t *state = seedRandGen();
	int i, degree = 6;



	coeffs = (mpz_t *) calloc(degree, sizeof(mpz_t));
	for(i = 0; i < degree; i ++)
	{
		mpz_init(coeffs[i]);
		mpz_urandomm(coeffs[i], *state, q);
	}

	polyToEval = setPolyWithArray(coeffs, q, degree);


	// printPoly(polyToEval);
	// printf("\n");
	// printPolyReverse(polyToEval);
	// printf("\n");

	for(i = 0; i < 10; i ++)
	{
		mpz_init_set_ui(x, i);
		y = evalutePoly(polyToEval, x, q);
		gmp_printf("-- %Zd\n", y);
	}
}