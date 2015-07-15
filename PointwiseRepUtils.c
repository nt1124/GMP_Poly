struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly, mpz_t q, int numPointToEval)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint, yPoint, *temp;
	int i;


	temp = (mpz_t *) calloc(1, sizeof(mpz_t));

	outputPoly -> numPoints = numPointToEval + 1;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));

	mpz_init_set_ui(xPoint, 1);
	mpz_init_set(outputPoly -> evalPoints[0], inputPoly -> coeffs[0]);


	for(i = 1; i <= outputPoly -> numPoints; i ++)
	{
		mpz_set_ui(xPoint, i);

		mpz_init(outputPoly -> evalPoints[i-1]);
		evalutePolyAlt(temp, inputPoly, xPoint, q);

		mpz_set(outputPoly -> evalPoints[i-1], *temp);
	}


	return outputPoly;
}


struct PointwiseRep *internalAddTwoPointwiseReps(struct PointwiseRep *polyA, struct PointwiseRep *polyB, int maxDegree, mpz_t q)
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


struct PointwiseRep *internalMultiplyTwoPointwiseReps(struct PointwiseRep *polyA, struct PointwiseRep *polyB, int maxDegree, mpz_t q)
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



struct Fq_poly *nLogN_MultiplyPolys(struct Fq_poly *rawPolyA, struct Fq_poly *rawPolyB, mpz_t q)
{
	struct PointwiseRep *processedPolyA, *processedPolyB, *outputPointwisePoly;
	struct Fq_poly *outputPoly;
	int maxDegreeDoubled;


	maxDegreeDoubled = 2 * MAX(rawPolyA -> degree, rawPolyB -> degree) + 1;
	processedPolyA = convertPolyToPointRep(rawPolyA, q, maxDegreeDoubled);
	processedPolyB = convertPolyToPointRep(rawPolyB, q, maxDegreeDoubled);

	outputPointwisePoly = internalMultiplyTwoPointwiseReps(processedPolyA, processedPolyB, maxDegreeDoubled, q);

	outputPoly = interpolatePointwiseRepMultiply(outputPointwisePoly, rawPolyA -> degree, rawPolyB -> degree, q);

	// outputPoly -> degree = rawPolyA -> degree + rawPolyB -> degree;
}

