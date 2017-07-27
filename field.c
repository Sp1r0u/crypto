#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gmp.h>

/******************************************
 * gcc field.c -L. -lpbc -lgmp -lm -o field 
 ******************************************/

void isPrime(mpz_t);
void getFieldOrder(mpz_t, unsigned int);

//===============================
//===============================

/************************************
 * p: prime
 * k: degree of the field extension
 *    k=1 if F is the base field  
 * dim: p^k is the order/dimension of the field
 *      i.e. F has p^k elements
 * base: base field                 
 ************************************/

struct field_t {
  mpz_t p, dim;
  unsigned int k;

  struct field_t *base;
  
  void (*isPrime)(mpz_t);
  void (*getFieldOrder)(mpz_t, unsigned int);
};

//===============================
//===============================

/*******************************************
 * an element a of a field is represented as
 *    a = a0 + a1*X + ... + a_{k-1}*X^{k-1}
 * n: degree of element x
 *    n is at most k-1
 *******************************************/

struct element_t {
  size_t n;

  mpz_t *x;
  
  struct field_t *F;
};

//===============================
//===============================

int main( int argc, char **argv ) {
  
  int i;
  size_t n = 10;
  gmp_randstate_t GMPRandState;

  char str[] = "313";
  
  struct field_t F;
  F.k = 1;

  mpz_init( F.p );
  
  mpz_set_str( F.p, str, 10 );

  isPrime( F.p );
  getFieldOrder( F );
  
  gmp_printf("field characteristic p:%Zd\n", F.p);
  
  struct element_t A;

  A.n = 3;

  A.x = malloc( A.n * sizeof(mpz_t) );

  A.F = &F;
  
  printf("size of x:%lu\n", A.n);
  gmp_printf("size base field:%Zd\n", A.F->p);
 
  /* initialize the random state with default algorithm */
  gmp_randinit_default( GMPRandState );

  /* seed the state with an unsigned long int */
  gmp_randseed_ui( GMPRandState, 1234567890 );

  /* bound for mpz_urandomm */
  mpz_t upBnd;
  mpz_init_set_str( upBnd, "18", 10 );
  
  for( i=0; i<A.n; i++ ) {
    mpz_init( A.x[i] );
    mpz_urandomm( A.x[i], GMPRandState, upBnd );
    mpz_out_str( stdout, 10, A.x[i] );
    printf("\n");
  }

  // freeing
  free( A.x );

  gmp_randclear( GMPRandState );
 
  return 0;
}

//===============================
//===============================

void getFieldOrder( mpz_t, unsigned int ) {
}

//===============================
//===============================

void isPrime( mpz_t x ) {
  /* call mpz_probab_prime_p to check 
   * whether x is prime
   * reps=2 if n is definitely prime 
   * reps=1 if n is probably prime (without being certain)
   * reps=0 if n is definitely composite */
  int reps;
  reps = mpz_probab_prime_p( x, reps );
  
  if (reps==2) { gmp_printf("%Zd is prime\n", x); }
}
