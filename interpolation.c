struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly, int maxDegree)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint;
	int i;

	outputPoly -> numPoints = maxDegree;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));


	mpz_init(xPoint);
	mpz_init_set(outputPoly -> evalPoints[0], inputPoly -> coeffs[0]);
	for(i = 1; i < outputPoly -> numPoints; i ++)
	{
		mpz_set_ui(xPoint, i);

	}

	return outputPoly;
}