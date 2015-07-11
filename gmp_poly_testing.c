// Testing Functions.



void test_evaluation(mpz_t q)
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



void test_interpolation(mpz_t q)
{
	struct Fq_poly *originalPoly, *interpolatedPoly;
	struct PointwiseRep *pointwiseVersion;

	mpz_t *coeffs;
	gmp_randstate_t *state = seedRandGen();
	int i, degree = 6, numCoeffs = degree + 1;


	coeffs = (mpz_t *) calloc(numCoeffs, sizeof(mpz_t));
	for(i = 0; i < numCoeffs; i ++)
	{
		mpz_init(coeffs[i]);
		mpz_urandomm(coeffs[i], *state, q);
	}
	originalPoly = setPolyWithArray(coeffs, q, degree);


	printPoly(originalPoly);
	printf("\n");

	pointwiseVersion = convertPolyToPointRep(originalPoly, q, degree);
	interpolatedPoly = interpolatePointwiseRep(pointwiseVersion, q);

	printPoly(interpolatedPoly);
	printf("\n");
}