struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly, mpz_t q, int maxDegree)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint, yPoint;
	int i;

	outputPoly -> numPoints = maxDegree;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));


	mpz_init_set_ui(xPoint, 1);
	mpz_init(yPoint);
	mpz_init_set(outputPoly -> evalPoints[0], inputPoly -> coeffs[0]);

	for(i = 1; i < outputPoly -> numPoints; i ++)
	{
		mpz_set_ui(xPoint, i);

		evalutePolyAlt(yPoint, inputPoly, xPoint, q);
	}

	return outputPoly;
}