// Function to print out a poly, least significant coeff first.
void printPoly(struct Fq_poly *polyToPrint)
{
	int i;

	// For each Coeff
	for(i = 0; i <= polyToPrint -> degree; i ++)
	{
		// If the coeff is not 0... (if it is no need to print it)
		if(0 != mpz_cmp_ui(polyToPrint -> coeffs[i], 0))
		{
			if(i == 0)
			{
				gmp_printf("%Zd", polyToPrint -> coeffs[i]);
			}
			else
			{
				gmp_printf("%Zd.X^%d", polyToPrint -> coeffs[i], i);
			}

			if(i < polyToPrint -> degree)
			{
				printf(" + ");
			}
		}
	}

	printf("\n");
}


// Function to print out a poly, most significant coeff first.
void printPolyReverse(struct Fq_poly *polyToPrint)
{
	int i;

	for(i = polyToPrint -> degree; i >= 0; i --)
	{

		if(0 != mpz_cmp_ui(polyToPrint -> coeffs[i], 0))
		{
			if(i == 0)
			{
				gmp_printf("%Zd", polyToPrint -> coeffs[i]);
			}
			else
			{
				gmp_printf("%Zd.X^%d", polyToPrint -> coeffs[i], i);
			}

			if(i != 0)
			{
				printf(" + ");
			}
		}
	}

	printf("\n");
}


// Initialise a polynomial with a given degree.
struct Fq_poly *initPolyWithDegree(int degree)
{
	struct Fq_poly *toReturn = (struct Fq_poly*) calloc(1, sizeof(struct Fq_poly));
	int i;


	toReturn -> degree = degree;
	toReturn -> coeffs = (mpz_t*) calloc(degree + 1, sizeof(mpz_t));

	for(i = 0; i <= degree; i ++)
	{
		mpz_init(toReturn -> coeffs[i]);
	}

	return toReturn;
}


// Clear the memory used for a polynomials coeffs.
void clearFqPoly(struct Fq_poly *toFree)
{
	int i;

	for(i = 0; i <= toFree -> degree; i ++)
	{
		mpz_clear(toFree -> coeffs[i]);
	}
}



struct Fq_poly *setPolyWithArray(mpz_t *inputArray, mpz_t q, int degree)
{
	struct Fq_poly *toReturn = initPolyWithDegree(degree);
	int i;


	for(i = 0; i <= degree; i ++)
	{
		mpz_mod(toReturn -> coeffs[i], inputArray[i], q);
	}

	return toReturn;
}


struct Fq_poly *addPolys_Internal(struct Fq_poly *x, struct Fq_poly *y, mpz_t q)
{
	struct Fq_poly *toReturn;
	mpz_t temp;
	int i, xDegree, yDegree;


	mpz_init(temp);
	xDegree = x -> degree;
	yDegree = y -> degree;
	toReturn = initPolyWithDegree(xDegree);

	for(i = 0; i <= yDegree; i ++)
	{
		mpz_add(temp, x -> coeffs[i], y -> coeffs[i]);
		mpz_mod(toReturn -> coeffs[i], temp, q);
	}

	for(; i <= xDegree; i ++)
	{
		mpz_mod(toReturn -> coeffs[i], x -> coeffs[i], q);
	}

	mpz_clear(temp);


	return toReturn;
}


struct Fq_poly *scalarMulti(struct Fq_poly *poly, mpz_t scalar, mpz_t q)
{
	struct Fq_poly *toReturn = initPolyWithDegree(poly -> degree);
	mpz_t temp;
	int i;

	mpz_init(temp);

	for(i = 0; i <= poly -> degree; i ++)
	{
		mpz_mul(temp, poly -> coeffs[i], scalar);
		mpz_mod(toReturn -> coeffs[i], temp, q);
	}

	mpz_clear(temp);

	return toReturn;
}


// 
void scalarMultiInPlace(struct Fq_poly *poly, mpz_t scalar, mpz_t q)
{
	mpz_t temp;
	int i;

	mpz_init(temp);

	for(i = 0; i <= poly -> degree; i ++)
	{
		mpz_mul(temp, poly -> coeffs[i], scalar);
		mpz_mod(poly -> coeffs[i], temp, q);
	}

	mpz_clear(temp);
}


// Function to wrap the functions for adding polynomials together.
struct Fq_poly *addPolys(struct Fq_poly *x, struct Fq_poly *y, mpz_t q)
{
	struct Fq_poly *toReturn;

	if(x -> degree >= y -> degree)
	{
		toReturn = addPolys_Internal(x, y, q);
	}
	else
	{
		toReturn = addPolys_Internal(y, x, q);
	}

	return toReturn;
}


// Performs naive multiplication of two Fq_polys.
struct Fq_poly *mulPolys(struct Fq_poly *x, struct Fq_poly *y, mpz_t q)
{
	struct Fq_poly *toReturn;
	mpz_t *tempArray, temp;
	int i, j, outputDegree = x -> degree + y -> degree;


	// Create a temporary mpz_t array
	tempArray = (mpz_t*) calloc(outputDegree + 1, sizeof(mpz_t));
	for(i = 0; i <= outputDegree; i ++)
	{
		mpz_init_set_ui(tempArray[i], 0);
	}

	// For each Coeff in the x Poly, multiply by...
	for(i = 0; i <= x -> degree; i ++)
	{
		// ...every coeff of the y Poly.
		for(j = 0; j <= y -> degree; j ++)
		{
			mpz_addmul(tempArray[i + j], x -> coeffs[i], y -> coeffs[j]);
		}
	}


	// Cut off anyleading zero coeffs.
	while(0 == tempArray[outputDegree])
	{
		outputDegree --;
	}

	// Put the tempArray into a poly.
	toReturn = setPolyWithArray(tempArray, q, outputDegree);

	for(i = 0; i <= x -> degree + y -> degree; i ++)
	{
		mpz_clear(tempArray[i]);
	}
	free(tempArray);

	return toReturn;
}


// Function to alter a single coeff of an Fq_poly struct.
void setPolyCoeff(struct Fq_poly *poly, mpz_t newCoeff, int coeffIndexToSet, mpz_t q)
{
	mpz_mod(poly -> coeffs[coeffIndexToSet], newCoeff, q);
}


// Evaluate an Fq_poly at a given point. Store output in the result pointer.
// The result pointer should be calloced and inited already.
void evalutePoly(mpz_t *result, struct Fq_poly *polyToEval, mpz_t x, mpz_t q)
{
	mpz_t temp, tempResult;
	int i;


	mpz_init(temp);
	mpz_set_ui(*result, 0);
	mpz_init(tempResult);

	// Standard Naive Polynomial evaluation. Could be improved upon?
	for(i = polyToEval -> degree; i >= 0; i--)
	{
		mpz_mul(temp, *result, x);
		mpz_add(tempResult, temp, polyToEval -> coeffs[i]);
		mpz_mod(*result, tempResult, q);
	}

	mpz_clear(temp);
	mpz_clear(tempResult);
}


// Function to take an list of mpz_t's and output the product of them.
void productOfMPZs(mpz_t output, mpz_t *inputArray, mpz_t q, int length)
{
	mpz_t unmoddedOutput;
	int i;


	// 
	mpz_init(output);
	mpz_init_set(unmoddedOutput, inputArray[0]);

	for(i = 1; i < length; i ++)
	{
		mpz_mul(unmoddedOutput, unmoddedOutput, inputArray[i]);
	}

	mpz_mod(output, unmoddedOutput, q);

	mpz_clear(unmoddedOutput);
}


// Function that takes a polynomial and trims the coeff array of leading zeroes.
void trimLeadingZeroes(struct Fq_poly *inputPoly)
{
	mpz_t *trimedCoeffs;
	int i, j;


	i = inputPoly -> degree;

	while(0 == mpz_cmp_ui(inputPoly -> coeffs[i], 0))
	{
		i --;
	}

	// Degree of the poly is actually i, thus we have i + 1 coeffs.
	trimedCoeffs = (mpz_t *) calloc(i + 1, sizeof(mpz_t));

	// For each valid coeff, copy into the new mpz_t array. Clear the old one as we go.
	for(j = 0; j <= i; j ++)
	{
		mpz_init_set(trimedCoeffs[j], inputPoly -> coeffs[j]);
		mpz_clear(inputPoly -> coeffs[j]);
	}
	// Clear the rest of the old coeff array.
	for(; j < inputPoly -> degree; j ++)
	{
		mpz_clear(inputPoly -> coeffs[j]);
	}

	// Housekeeping
	free(inputPoly -> coeffs);

	// Give the pointer to the new coeff array etc. 
	inputPoly -> degree = i;
	inputPoly -> coeffs = trimedCoeffs;
}


// Function to clear/free the memory used for an Fq_poly structure.
void freeFq_Poly(struct Fq_poly *polyToFree)
{
	int i;

	// For each coeff, clear the mpz_t's
	for(i = 0; i <= polyToFree -> degree; i ++)
	{
		mpz_clear(polyToFree -> coeffs[i]);
	}

	free(polyToFree -> coeffs);
	free(polyToFree);
}