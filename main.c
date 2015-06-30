#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include "gmp_poly.h"


// Primarily here in order for debugging purposes.
int main()
{
	const char *qStr = "A9FB57DBA1EEA9BC3E660A909D838D726E3BF623D52620282013481D1F6E5377";
	gmp_sscanf(qStr, "%ZX", q);


	test_evaluation();


	return 0;
}