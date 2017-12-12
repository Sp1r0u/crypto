#ifndef _POLY_H_
#define _POLY_H_

#include <gmp.h>

#include "circuit.h"

void evaluate_Lagrange_polynomial( struct circuit_t*, mpz_t );

#endif /* _POLY_H_ */
