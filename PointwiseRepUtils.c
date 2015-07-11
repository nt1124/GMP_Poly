struct PointwiseRep *convertPolyToPointRep(struct Fq_poly *inputPoly, mpz_t q, int numPointToEval)
{
	struct PointwiseRep *outputPoly = (struct PointwiseRep *) calloc(1, sizeof(struct PointwiseRep));
	mpz_t xPoint, yPoint;
	int i;


	outputPoly -> numPoints = numPointToEval + 1;
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


struct PointwiseRep *addTwoPointwiseReps(struct PointwiseRep *polyA, struct PointwiseRep *polyB, int maxDegree, mpz_t q)
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


struct PointwiseRep *multiplyTwoPointwiseReps(struct PointwiseRep *polyA, struct PointwiseRep *polyB, int maxDegree, mpz_t q)
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




/*
struct Fq_poly *generateLagrangePoly(int n, unsigned int i, mpz_t q)
{
	struct Fq_poly **factors, **intermediates, *output;
	mpz_t divFactor, divFactorInv;
	int j;


	intermediates = (struct Fq_poly **) calloc(n - 1, sizeof(struct Fq_poly*));
	mpz_init(divFactorInv);


	factors = getLagrangeFactorsDivProduct(divFactor, n, delta_i, i, q);

	if(n > 1)
	{
		intermediates[0] = mulPolys(factors[0], factors[1], q);


		for(j = 1; j < n - 1; j ++)
		{
			intermediates[j] = mulPolys(intermediates[j - 1], factors[j + 1], q);
		}
	}
	else
	{
		intermediates[0] = factors[0];
	}

	mpz_invert(divFactorInv, divFactor, q);


	output = scalarMulti(intermediates[n - 2], divFactorInv, q);

	return output;
}



struct Fq_poly **generateFixedLagrangePolys(int length, mpz_t q)
{
	struct Fq_poly **lagrangePolys;
	int i;


	lagrangePolys = (struct Fq_poly **) calloc(length, sizeof(struct Fq_poly *));

	for(i = 1; i <= length; i ++)
	{
		lagrangePolys[i-1] = generateLagrangePoly(length, delta_i, i, q);
	}

	return lagrangePolys;
}


struct Fq_poly *interpolatePointwiseRep(struct PointwiseRep *polyToInterpolate, mpz_t q)
{
	struct Fq_poly *outputPoly = (struct Fq_poly *) calloc(1, sizeof(struct Fq_poly));
	struct Fq_poly **lagrangePolys;
	int *xPoints, i;


	xPoints = (int *) calloc(polyToInterpolate -> numPoints, sizeof(int));
	lagrangePolys = generateAllLagrangePolys(polyToInterpolate -> numPoints, xPoints, q);
	outputPoly = initPolyWithDegree(polyToInterpolate -> numPoints); 

	for(i = 1; i <= polyToInterpolate -> numPoints; i ++)
	{
		scalarMultiInPlace(lagrangePolys[i-1], polyToInterpolate -> evalPoints[i-1], q);
		outputPoly = addPolys(lagrangePolys[i-1], outputPoly, q);
	}


	return outputPoly;
}

/*
struct Fq_poly *getPolyFromCodewords(mpz_t *codewords, int *delta_i, int length, mpz_t q)
{
	struct Fq_poly **lagrangePolys, *output;
	int i;


	lagrangePolys = generateAllLagrangePolys(length, delta_i, q);
	output = initPolyWithDegree(length); 


	for(i = 1; i <= length; i ++)
	{
		scalarMultiInPlace(lagrangePolys[i-1], codewords[i-1], q);
		output = addPolys(lagrangePolys[i-1], output, q);
	}

	return output;
}
*/




