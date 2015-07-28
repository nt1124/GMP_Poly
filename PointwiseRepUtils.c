// This function will convert a standard Coeffient representation poly to a pointwise representation.
struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly, mpz_t q, int numPointToEval)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint, yPoint, *temp;
	int i;


	// We need degree + 1 many points to represent the poly.
	// So initialise a pointwise rep struct of the required size.
	outputPoly -> numPoints = numPointToEval + 1;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));

	// Housekeeping
	temp = (mpz_t *) calloc(1, sizeof(mpz_t));
	mpz_init(*temp);
	mpz_init(xPoint);

	// Then for 1 up to numPoints (inclusive) eval the poly and store the result.
	// This gives us numPoints many values 
	for(i = 1; i <= outputPoly -> numPoints; i ++)
	{
		// Set X in an mpz_t variable
		mpz_set_ui(xPoint, i);

		// Evaluate the poly and store the result in a temp mpz_t
		// mpz_init(outputPoly -> evalPoints[i-1]);
		evalutePoly(temp, inputPoly, xPoint, q);

		// Put the output in the 
		mpz_init_set(outputPoly -> evalPoints[i-1], *temp);
	}


	mpz_clear(xPoint);
	mpz_clear(yPoint);
	mpz_clear(*temp);
	free(temp);

	return outputPoly;
}



void freePointwiseRep_Poly(struct PointwiseRep *polyToFree)
{
	int i;

	for(i = 0; i < polyToFree -> numPoints; i ++)
	{
		mpz_clear(polyToFree -> evalPoints[i]);
	}

	free(polyToFree -> evalPoints);
	free(polyToFree);
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

	freePointwiseRep_Poly(processedPolyA);
	freePointwiseRep_Poly(processedPolyB);
	freePointwiseRep_Poly(outputPointwisePoly);

	return outputPoly;
}
