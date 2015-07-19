struct Fq_poly **getLagrangeFactorsDivProduct(mpz_t divFactorProduct, int n, unsigned int i, mpz_t q)
{
	struct Fq_poly **factors = (struct Fq_poly **) calloc(n, sizeof(struct Fq_poly *));
	mpz_t *inputArray = (mpz_t *) calloc(2, sizeof(mpz_t));
	mpz_t *divFactors = (mpz_t *) calloc(n, sizeof(mpz_t));
	mpz_t temp, x_i;
	unsigned int j;


	mpz_init(temp);
	mpz_init(inputArray[0]);
	mpz_init_set_ui(inputArray[1], 1);
	mpz_init_set_ui(x_i, i);

	for(j = 1; j <= n; j ++)
	{
		mpz_init(divFactors[j-1]);

		if(i != j)
		{
			mpz_sub_ui(temp, x_i, j);
			mpz_mod(divFactors[j-1], temp, q);
			mpz_sub_ui(inputArray[0], q, j);

			factors[j-1] = setPolyWithArray(inputArray, q, 1);
		}
		else
		{
			mpz_set_ui(divFactors[j-1], 1);
			mpz_set_ui(inputArray[0], 1);
			factors[j-1] = setPolyWithArray(inputArray, q, 0);
		}
	}


	productOfMPZs(divFactorProduct, divFactors, q, n);

	for(i = 0; i < n; i ++)
	{
		mpz_clear(divFactors[i]);
	}

	mpz_clear(inputArray[0]);
	mpz_clear(inputArray[1]);
	mpz_clear(temp);
	mpz_clear(x_i);

	free(inputArray);
	free(divFactors);

	return factors;
}


void tidyUpLagrangePolys(struct Fq_poly **factors, struct Fq_poly **intermediates, int n)
{
	int i;

	for(i = 0; i < n; i ++)
	{
		freeFq_Poly(factors[i]);
	}
	free(factors);

	for(i = 0; i < n - 1; i ++)
	{
		freeFq_Poly(intermediates[i]);
	}
	free(intermediates);
}


struct Fq_poly *generateLagrangePoly(int n, unsigned int i, mpz_t q)
{
	struct Fq_poly **factors, **intermediates, *output;
	mpz_t divFactor, divFactorInv;
	int j;


	intermediates = (struct Fq_poly **) calloc(n - 1, sizeof(struct Fq_poly*));
	mpz_init(divFactorInv);

	factors = getLagrangeFactorsDivProduct(divFactor, n, i, q);

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

	tidyUpLagrangePolys(factors, intermediates, n);
	mpz_clear(divFactor);
	mpz_clear(divFactorInv);

	return output;
}


struct Fq_poly **generateLagrangePolys(int length, mpz_t q)
{
	struct Fq_poly **lagrangePolys;
	int i;


	lagrangePolys = (struct Fq_poly **) calloc(length, sizeof(struct Fq_poly *));

	for(i = 1; i <= length; i ++)
	{
		lagrangePolys[i-1] = generateLagrangePoly(length, i, q);
	}

	return lagrangePolys;
}


struct Fq_poly *interpolatePointwiseRep(struct PointwiseRep *polyToInterpolate, mpz_t q)
{
	struct Fq_poly *outputPoly = (struct Fq_poly *) calloc(1, sizeof(struct Fq_poly));
	struct Fq_poly **lagrangePolys;
	int i;


	lagrangePolys = generateLagrangePolys(polyToInterpolate -> numPoints, q);
	outputPoly = initPolyWithDegree(polyToInterpolate -> numPoints - 1); 

	for(i = 0; i < polyToInterpolate -> numPoints; i ++)
	{
		scalarMultiInPlace(lagrangePolys[i], polyToInterpolate -> evalPoints[i], q);
		outputPoly = addPolys(lagrangePolys[i], outputPoly, q);

		freeFq_Poly(lagrangePolys[i]);
	}


	return outputPoly;
}


struct Fq_poly *interpolatePointwiseRepMultiply(struct PointwiseRep *polyToInterpolate, int degreeA, int degreeB, mpz_t q)
{
	struct Fq_poly *outputPoly, *tempPoly;// = (struct Fq_poly *) calloc(1, sizeof(struct Fq_poly));
	struct Fq_poly **lagrangePolys;
	int i;


	lagrangePolys = generateLagrangePolys(polyToInterpolate -> numPoints, q);
	outputPoly = initPolyWithDegree(degreeA + degreeB); 


	for(i = 0; i < polyToInterpolate -> numPoints; i ++)
	{
		scalarMultiInPlace(lagrangePolys[i], polyToInterpolate -> evalPoints[i], q);
		tempPoly = addPolys(lagrangePolys[i], outputPoly, q);

		freeFq_Poly(lagrangePolys[i]);
		freeFq_Poly(outputPoly);
		outputPoly = tempPoly;
	}


	trimLeadingZeroes(outputPoly);
	free(lagrangePolys);

	return outputPoly;
}
