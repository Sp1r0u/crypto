#ifndef _FIELD_H_
#define _FIELD_H_

#include <gmp.h>

#include "BN.h"

/******************************************
 * p: prime number - order of the field
 *    i.e. field has p elements (0,...,p-1)
 *****************************************/

struct field_t {
  mpz_t p; //order of G1, G2, GT
};

void init_field( struct field_t*, struct curve_t* );
void rnd_field_elt( mpz_t, struct field_t*, gmp_randstate_t); //returns a randomly selected non-zero field element

#endif /* _FIELD_H_ */
