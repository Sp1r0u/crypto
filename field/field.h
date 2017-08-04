#ifndef _FIELD_H_
#define _FIELD_H_

#include <stdbool.h>

#include "element.h"

/************************************
 * p: prime 
 * k: degree of the field extension
 *    k=1 if F is the base field  
 * n: n=p^k the order/dimension of the field
 *    i.e. F has p^k elements
 * base: base field                 
 ************************************/

struct field_t {
  mpz_t p;
  mpz_t n;
  unsigned int k;

  struct field_t *base;

  struct element_t *irr_poly;
  
};

bool initialize( struct field_t*, const char*, unsigned int, struct field_t* );
bool is_field_order_prime( mpz_t );

#endif /* _FIELD_H_ */
