#ifndef _BN_H_
#define _BN_H_

#include "/home/elmagnifico/pbc-0.5.14/include/pbc.h"

/**************************************
 * type a, d, e, f, g, a1, i
 * q order of the finite field Fq 
 * r order of the curve 
 * b E: y^2 = x^3 + b
 * beta quadratic nonresidue in Fq
 **************************************/

struct curve_t {
  mpz_t x;
  mpz_t tx; //t(x): trace of the curve
  mpz_t nx; //n(x): curve order
  mpz_t px; //p(x): characteristic of base field
  mpz_t r; //order of G1, G2, GT (see include/pbc_pairing.h line 18)

  element_t P; //<P>=G1
  element_t Q; //<Q>=G2
  
  pbc_param_t pbc_params;

  pairing_t   pairing;
 
};

void init_curve( struct curve_t* );

void free_curve( struct curve_t* );

#endif /* _BN_H_ */
