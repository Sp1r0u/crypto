#ifndef _FIELD_H_
#define _FIELD_H_

#include <stdbool.h>
#include <gmp.h>

#include "EC.h"

extern bool _debug;
extern unsigned int _order;

struct field_t;
struct field_elt_t;

//================================
//================================

struct field_t {
  struct curve_t *EC;
  
  void (*init) (struct field_t*, struct curve_t*);
  void (*free) (struct field_t*);

  // field operations
  void (*rndElt) (struct field_elt_t*, gmp_randstate_t);
  void (*modElt) (struct field_elt_t*);
  void (*negElt) (struct field_elt_t*, struct field_elt_t*);
  void (*invElt) (struct field_elt_t*, struct field_elt_t*);
  void (*mulElt) (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);
  void (*divElt) (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);
  void (*addElt) (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);
  void (*subElt) (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);

  void (*powuiElt) (struct field_elt_t*, struct field_elt_t*, unsigned int); 
  
};

//================================
//================================

struct field_elt_t {
  
  mpz_t value;
  struct field_t *field;

  void (*init) (struct field_elt_t*, struct field_t*);
  
  void (*view) (struct field_elt_t*);
  void (*free) (struct field_elt_t*);
};

//================================
//================================

// initialize algebraic field
void initField (struct field_t*, struct curve_t*);

// initialize element of the field
void initFieldElt (struct field_elt_t*, struct field_t*);

void viewFieldElt (struct field_elt_t*);

// returns a randomly selected non-zero field element
void rndFieldElt (struct field_elt_t*, gmp_randstate_t);

// returns the modulo operation
void modFieldElt (struct field_elt_t*);

// returns the negative a<-(-a)
void negFieldElt (struct field_elt_t*, struct field_elt_t*);

// b<-a^-1 [mod r]
void invFieldElt (struct field_elt_t*, struct field_elt_t*);

// c<-a*b [mod r]
void mulFieldElt (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);

// c<-a/b [mod r]
void divFieldElt (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);

// c<-a+b [mod r]
void addFieldElt (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);

// c<-a-b [mod r]
void subFieldElt (struct field_elt_t*, struct field_elt_t*, struct field_elt_t*);

// b<-a^i [mod r], i:integer
void powuiFieldElt (struct field_elt_t*, struct field_elt_t*, unsigned int);

// free field
void freeField (struct field_t*);

// free field elt
void freeFieldElt (struct field_elt_t*);

#endif /* _FIELD_H_ */
