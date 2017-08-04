#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <gmp.h>

#include "field.h"
#include "toolbox.h"

bool is_field_order_prime( mpz_t x ) {

  /* call mpz_probab_prime_p to check whether x is prime
   * reps=2 if n is definitely prime 
   * reps=1 if n is probably prime (without being certain)
   * reps=0 if n is definitely composite */

  int reps = -99;

  reps = mpz_probab_prime_p( x, reps );

  if( reps==2 ) { return true; }

  else { return false; }

}

//==================================================
//==================================================

bool initialize( struct field_t *F, const char *b, unsigned int k, struct field_t *Fb) {

  mpz_init( F->p );
  mpz_init( F->n );

  F->k = k;
  
  if( b != NULL ) {
      mpz_set_str( F->p, b, 10 );
      mpz_pow_ui( F->n, F->p, k );

      if( is_field_order_prime( F->n ) != true ) {
	gmp_printf("the order (%Zd) of the finite field is not a prime or a power of a prime... exiting now\n", F->p);
	exit(0);
      }
  }

  else {
    mpz_set( F->p, Fb->p );
    mpz_pow_ui( F->n, F->p, k );
  }
  
  //gmp_printf("p: %Zd\n", F->p);
  gmp_printf("n: %Zd\n", F->n);
  
  F->base = Fb;

  //if (F->base != NULL ) { printf("%d\n", (F->base)->k ); }

  return true;
}

//==================================================
//==================================================

/*void generate_random_monic_polynomial( struct field_t *F,
				       struct element_t **temp,
				       unsigned int p,
				       gmp_randstate_t GMPRandState ) {

  int i; 
  mpz_t coeff;
  mpz_init( coeff );
  
  /* bound for mpz_urandomm */
/*  mpz_t upBnd;

  if (F->base == NULL) {

    mpz_init_set( upBnd, F->n );
     
    for (i=p; i>=0; i--) {
      if (i==p) { mpz_init_set_ui( coeff, 1 ); }
      else      { mpz_urandomm( coeff, GMPRandState, upBnd ); }
      
      mpz_out_str( stdout, 10, coeff );
      
      printf("\n");
      
      //E->exp = p;
      
      //if (i==0) { E->next = NULL; }

      struct element_t *r;
      struct element_t *z;
      z = *temp;

      if (z==NULL) {
	r = (struct element_t*)malloc(sizeof(struct element_t));
	mpz_set( r->coeff, coeff );
	r->exp = exp;
	*temp = r;
	r->next = (struct element_t*)malloc(sizeof(struct element_t));
	r = r->next;
	r->next = NULL;
      }
      else {
	mpz_set( r->coeff, coeff );
	r->exp = exp;
	r->next = (struct element_t*)malloc(sizeof(struct element_t));
	r = r->next;
	r->next = NULL;
      }
    }
  }

  mpz_clear( upBnd );

}

//==================================================
//==================================================

void create_element( mpz_t coeff, unsigned int exp, struct element_t **temp ) {
  /*source: http://www.geeksforgeeks.org/adding-two-polynomials-using-linked-list/ */

/*  struct element_t *r;
  struct element_t *z;
  z = *temp;

  if (z==NULL) {
    r = (struct element_t*)malloc(sizeof(struct element_t));
    mpz_set( r->coeff, coeff );
    r->exp = exp;
    *temp = r;
    r->next = (struct element_t*)malloc(sizeof(struct element_t));
    r = r->next;
    r->next = NULL;
  }
  else {
    mpz_set( r->coeff, coeff );
    r->exp = exp;
    r->next = (struct element_t*)malloc(sizeof(struct element_t));
    r = r->next;
    r->next = NULL;
  }
  
}
*/
