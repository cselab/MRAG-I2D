#include "binomial.h"

//Binomial Coefficients for the 2D Vortex Case:
unsigned long long binomial(unsigned n, unsigned k) {
    if (k > n)
        return 0;
	
    if (k > n/2)
        k = n-k; // Take advantage of symmetry
	
    long double accum = 1;
    for (unsigned i = 1; i <= k; i++)
		accum = accum * (n-k+i) / i;
	
    return accum + 0.5; // avoid rounding error
}
