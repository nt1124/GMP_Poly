struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint;
	int i;

	outputPoly -> numPoints = inputPoly -> degree;
	outputPoly -> evalPoints = (mpz_t *) calloc(outputPoly -> numPoints, sizeof(mpz_t));

	mpz_init_set_ui(xPoint, 1);
	mpz_init_set_ui(outputPoly -> evalPoints[0], inputPoly -> coeffs[0]);


	for(i = 1; i < outputPoly -> degree; i ++)
	{

	}

	return outputPoly;
}