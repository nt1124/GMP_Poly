// Testing Functions.



void test_evaluation()
{
	struct Fq_poly *polyToEval;

	mpz_t *coeffs, x, *y, yAlt;
	gmp_randstate_t *state = seedRandGen();
	int i, degree = 6, numCoeffs = degree + 1;


	mpz_init(yAlt);

	// Note the difference between the degree and the number of coeffients.
	// Number of coeffients is the degree plus one.
	coeffs = (mpz_t *) calloc(numCoeffs, sizeof(mpz_t));
	for(i = 0; i < numCoeffs; i ++)
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

		mpz_set_ui(yAlt, 0);
		evalutePolyAlt(yAlt, polyToEval, x, q);
		gmp_printf("++ %Zd\n\n", yAlt);
	}
}