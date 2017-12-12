#ifndef _FIELD_H_
#define _FIELD_H_

#include <stdbool.h>
#include <gmp.h>

#include "BN.h"

extern bool fDebug;
extern unsigned int fOrder;

/******************************************
 * p: prime number - order of the field
 *    i.e. field has p elements (0,...,p-1)
 *****************************************/

struct field_t {
  mpz_t p; //order of G1, G2, GT
};

struct field_elt_t {
  mpz_t value;
  struct field_t *field;
};

struct poly_ring_node_t {
  unsigned int exponent;
  struct field_elt_t *coefficient;
  struct poly_ring_node_t *next;
};

struct poly_ring_t {
  int length;
  struct poly_ring_node_t *head;
  struct poly_ring_node_t *tail;
};

// initialize algebraic field
void init_field( struct field_t*, struct curve_t* );

// returns a randomly selected non-zero field element
void nonzero_rnd_field_elt( struct field_elt_t*, gmp_randstate_t );

// returns the modulo operation
void mod_field_elt( struct field_elt_t* );

// c<-a*b [mod p]
void mul_field_elt( struct field_elt_t*, struct field_elt_t*, struct field_elt_t* );

// c<-a/b [mod p]
void div_field_elt( struct field_elt_t*, struct field_elt_t*, struct field_elt_t* );

// c<-a+b [mod p]
void add_field_elt( struct field_elt_t*, struct field_elt_t*, struct field_elt_t* );

// c<-a-b [mod p]
void sub_field_elt( struct field_elt_t*, struct field_elt_t*, struct field_elt_t* );

// b<-a^-1 [mod p]
void inv_field_elt( struct field_elt_t*, struct field_elt_t* );

void insert_poly_ring_node( mpz_t, unsigned int, struct field_t*, struct poly_ring_t* );

void delete_poly_ring_node( struct poly_ring_t*, struct poly_ring_node_t* );

void print_poly( struct poly_ring_t* );

// p(x)<-q(x)*r(x)
void mul_poly( struct poly_ring_t*, struct poly_ring_t*, struct poly_ring_t* );

// p(x)<-q(x)+r(x)
void add_poly( struct poly_ring_t*, struct poly_ring_t*, struct poly_ring_t* );

// p(x)<-q(x)-r(x)
void sub_poly( struct poly_ring_t*, struct poly_ring_t*, struct poly_ring_t* );

// p(x)<-q(x)*(r(x)^-1)
void div_poly( struct poly_ring_t*, struct poly_ring_t*, struct poly_ring_t* );

// q(x)<-c*p(x) c:constant
void cst_mul_poly( struct poly_ring_t*, struct field_elt_t*, struct poly_ring_t* );

void evaluate_poly( struct poly_ring_t*, struct field_elt_t*, struct field_elt_t* );

void poly_ordering( struct poly_ring_t* );

unsigned int get_largest_exp( struct poly_ring_t* );

// generates a random monic polynomial
void random_monic_poly( struct field_t*, struct poly_ring_t*, unsigned int, gmp_randstate_t );

// generates a random polynomial
void random_poly( struct field_t*, struct poly_ring_t*, unsigned int, gmp_randstate_t );

// testing poly for irreducibility
bool is_poly_irreducible( struct field_t*, struct poly_ring_t* );

// freeing memory
void free_field( struct field_t* );

void free_poly_ring( struct poly_ring_node_t* );

void free_field_elt( struct field_elt_t* );

#endif /* _FIELD_H_ */
