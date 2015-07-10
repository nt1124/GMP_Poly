struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly, mpz_t q, int maxDegree)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint, yPoint;
	int i;


	outputPoly -> numPoints = maxDegree;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));

	mpz_init_set_ui(xPoint, 1);
	mpz_init_set(outputPoly -> evalPoints[0], inputPoly -> coeffs[0]);

	for(i = 1; i < outputPoly -> numPoints; i ++)
	{
		mpz_set_ui(xPoint, i);

		mpz_init(outputPoly -> evalPoints[i]);
		evalutePolyAlt(outputPoly -> evalPoints[i], inputPoly, xPoint, q);
	}


	return outputPoly;
}


struct PointwiseRep *addTwoPointwiseReps(struct PointwiseRep *polyA, struct PointwiseRep *polyB, int maxDegree)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t tempMPZ;
	int i;


	mpz_init(tempMPZ);
	outputPoly -> numPoints = maxDegree;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));

	for(i = 0; i < outputPoly -> numPoints; i ++)
	{
		mpz_init(outputPoly -> evalPoints[i]);
		mpz_add(tempMPZ, polyA -> evalPoints[i], polyB -> evalPoints[i]);
		mpz_mod(outputPoly -> evalPoints[i], tempMPZ, q);
	}

	mpz_clear(tempMPZ);


	return outputPoly;
}


struct PointwiseRep *multiplyTwoPointwiseReps(struct PointwiseRep *polyA, struct PointwiseRep *polyB, int maxDegree)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t tempMPZ;
	int i;


	mpz_init(tempMPZ);
	outputPoly -> numPoints = maxDegree;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));

	for(i = 0; i < outputPoly -> numPoints; i ++)
	{
		mpz_init(outputPoly -> evalPoints[i]);
		mpz_mul(tempMPZ, polyA -> evalPoints[i], polyB -> evalPoints[i]);
		mpz_mod(outputPoly -> evalPoints[i], tempMPZ, q);
	}

	mpz_clear(tempMPZ);


	return outputPoly;
}


