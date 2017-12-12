#ifndef _EC_H_
#define _EC_H_

#include "/home/elmagnifico/pbc-0.5.14/include/pbc.h"

struct curve_t {
  mpz_t x;
  mpz_t tx; // t(x): trace of the curve
  mpz_t nx; // n(x): curve order
  mpz_t px; // p(x): characteristic of base field
  mpz_t r;  // order of G1, G2, GT (see include/pbc_pairing.h line 18)

  element_t P; // <P>=G1
  element_t Q; // <Q>=G2
  
  pbc_param_t pbc_params;

  pairing_t   pairing;

  void (*init) (struct curve_t*);
  void (*view) (struct curve_t*);
  void (*free) (struct curve_t*);

};

//================================
//================================

void initEC (struct curve_t*); 
void viewEC (struct curve_t*);
void freeEC (struct curve_t*);

#endif /* _EC_H_ */
