#ifndef _POLY_H_
#define _POLY_H_

#include <gmp.h>

#include "field.h"

struct poly_node_t {
  unsigned int exponent;
  struct field_elt_t *coefficient;
  struct poly_node_t *next;
};

//================================
//================================

struct poly_t {
  unsigned int length;
  
  struct poly_node_t *head;
  struct poly_node_t *tail;

  void (*init) (struct poly_t*);
  void (*view) (struct poly_t*);
  void (*free) (struct poly_t*);

  void (*reset) (struct poly_t*);
  
  void (*rndPoly) (struct poly_t*, struct field_t*, unsigned int, gmp_randstate_t);

  void (*insertNode) (struct poly_t*, mpz_t, unsigned int, struct field_t*);
  void (*removeNode) (struct poly_t*, struct poly_node_t*); 
  
  void (*mulPoly) (struct poly_t*, struct poly_t*, struct poly_t*);
  void (*divPoly) (struct poly_t*, struct poly_t*, struct poly_t*);
  void (*addPoly) (struct poly_t*, struct poly_t*, struct poly_t*);
  void (*subPoly) (struct poly_t*, struct poly_t*, struct poly_t*);

  void (*cstMulPoly) (struct poly_t*, struct field_elt_t*, struct poly_t*);
  
  void (*simplifyPoly) (struct poly_t*);
  void (*sort) (struct poly_t*);
  
  void (*copyPoly) (struct poly_t*, struct poly_t*);

  unsigned int (*hiExp) (struct poly_t*);

  void (*eval) (struct poly_t*, struct field_elt_t*, struct field_elt_t*);
  
};

//================================
//================================

void initPoly (struct poly_t*);
void viewPoly (struct poly_t*);
void freePoly (struct poly_t*);

void resetPoly (struct poly_t*);
  
// generates a random polynomial
void randomPoly (struct poly_t*, struct field_t*, unsigned int, gmp_randstate_t);

void insertPolyNode (struct poly_t*, mpz_t, unsigned int, struct field_t*);
void removePolyNode (struct poly_t*, struct poly_node_t*);

// p(x)<-q(x)*r(x)
void multiplyPoly (struct poly_t*, struct poly_t*, struct poly_t*);

// p(x)<-q(x)+r(x)
void addPoly (struct poly_t*, struct poly_t*, struct poly_t*);

// p(x)<-q(x)-r(x)
void subtractPoly (struct poly_t*, struct poly_t*, struct poly_t*);

// p(x)<-q(x)/r(x)
void dividePoly (struct poly_t*, struct poly_t*, struct poly_t*);

// q(x)<-c*p(x) c:constant
void cstMultiplyPoly (struct poly_t*, struct field_elt_t*, struct poly_t*);

void simplifyPoly (struct poly_t*);

void sortPoly (struct poly_t*);

void copyPoly (struct poly_t*, struct poly_t*);

// returns poly. highest exp. 
unsigned int getHighestExp (struct poly_t*);

// evaluates y = p(x)
void evaluatePoly (struct poly_t*, struct field_elt_t*, struct field_elt_t*);

#endif /* _POLY_H_ */
