#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gmp.h>
#include "../../pbc-0.5.14/include/pbc.h"
#include "../../pbc-0.5.14/include/pbc_test.h"

/*****************************************
 * gcc test.c -L. -lpbc -lgmp -lm -o test 
 ****************************************/

struct curve_params {
  mpz_t x;
  mpz_t tx; //t(x): trace of the curve
  mpz_t nx; //n(x): curve order
  mpz_t px; //p(x): characteristic of base field
  pbc_param_t params;
  pairing_t pairing;
};

void set_params( struct curve_params* );
void init_pairing( struct curve_params* );
void get_rnd_element_of_Fp( element_t*, struct curve_params* );
void get_rnd_element_of_G1( element_t*, struct curve_params* );
void get_rnd_element_of_G2( element_t*, struct curve_params* );

/**************************************
 * type a, d, e, f, g, a1, i
 * q order of the finite field Fq 
 * r order of the curve 
 * b E: y^2 = x^3 + b
 * beta quadratic nonresidue in Fq
 **************************************/

int main( int argc, char **argv ) {

  struct curve_params BN; 

  set_params( &BN );
  
  init_pairing( &BN );

  double time1, time2;

  time1 = pbc_get_time();
  
  element_t a, b;
  get_rnd_element_of_Fp( &a, &BN );
  get_rnd_element_of_Fp( &b, &BN );

  element_t P, Q, R;
  get_rnd_element_of_G1( &P, &BN );
  get_rnd_element_of_G2( &Q, &BN );
  element_init_GT( R, BN.pairing );
  element_random( R );
  element_printf("R = %B\n", R);
  
  element_printf("a = %B\n", a);
  element_printf("b = %B\n", b);
  element_printf("P = %B\n", P);
  element_printf("Q = %B\n", Q);

  element_t aP, bQ;
  element_init_G1( aP, BN.pairing );
  element_init_G2( bQ, BN.pairing );
  element_mul_zn( aP, P, a ); //[a]P = P + P + ... + P (a times)
  element_mul_zn( bQ, Q, a ); //[b]Q = Q + Q + ... + Q (b times)
  
  element_printf("aP = %B\n", aP);
  element_printf("bQ = %B\n", bQ);

  element_pairing(R, aP, bQ);

  element_printf("R = %B\n", R);
  
  time2 = pbc_get_time();
  printf("execution time =%fs\n", time2 - time1);
  
  /* freeing memory */
  element_clear( a );
  element_clear( b );
  element_clear( P );
  element_clear( Q );
  element_clear( aP );
  element_clear( bQ );
  pairing_clear( BN.pairing );

  return 0;
}

//================================================
//================================================

void get_rnd_element_of_G2( element_t *a, struct curve_params* c ) {
  element_init_G2( *a, c->pairing );
  element_random( *a );
  //element_printf("elt = %B\n", *a);
}

//================================================
//================================================

void get_rnd_element_of_G1( element_t *a, struct curve_params* c ) {
  element_init_G1( *a, c->pairing );
  element_random( *a );
  //element_printf("elt = %B\n", *a);
}

//================================================
//================================================

void get_rnd_element_of_Fp( element_t *a, struct curve_params* c ) {
  element_init_Zr( *a, c->pairing );
  element_random( *a );
  //element_printf("elt = %B\n", *a);
}

//================================================
//================================================

void init_pairing( struct curve_params* c ) {
  printf("--- calling init_pairing ---\n");
  pairing_init_pbc_param( c->pairing, c->params );
  if( pairing_is_symmetric( c->pairing ) == 0 ) { printf("pairing function is not symmetric\n"); }
  else { printf("pairing function is symmetric..., exiting now!\n"); exit(1); }

  printf("length in bytes to represent an element of G1: %d\n", pairing_length_in_bytes_G1(c->pairing));
  printf("length in bytes to represent an element of G2: %d\n", pairing_length_in_bytes_G2(c->pairing));

}

//================================================
//================================================

void set_params( struct curve_params* c ) {
  const char * foo;
  mpz_init( c->x );
  mpz_init( c->tx );
  mpz_init( c->nx );
  mpz_init( c->px );

  foo = "1";
  
  if( mpz_set_str( c->x, foo, 10 ) != 0 ) { printf("error... exiting now!\n"); exit(1); }
  mpz_neg( c->x, c->x );
  if( mpz_set_str( c->tx, foo, 10 ) != 0 ) { printf("error... exiting now!\n"); exit(1); }
  mpz_neg( c->tx, c->tx );
  
  foo = "205523667896953300194896352429254920972540065223";
  if( mpz_set_str( c->px, foo, 10 ) != 0 ) { printf("error... exiting now!\n"); exit(1); }

  foo = "205523667896953300194895899082072403858390252929";
  if( mpz_set_str( c->nx, foo, 10 ) != 0 ) { printf("error... exiting now!\n"); exit(1); }

  gmp_printf("x  = %Zd\n", c->x);
  gmp_printf("tx = %Zd\n", c->tx);
  gmp_printf("px = %Zd\n", c->px);
  gmp_printf("nx = %Zd\n", c->nx);

  char str1[1024];
  char str2[1024];
  /* type of the ECC */
  strcat( str1, "type f\n"); 
  /* field characteristic */
  strcat( str1, "q ");
  mpz_get_str( str2, 10, c->px ); 
  strcat( str1, str2);
  strcat( str1, "\n");
  /* order of the curve */
  strcat( str1, "r ");
  mpz_get_str( str2, 10, c->nx );
  strcat( str1, str2);
  strcat( str1, "\n");
  /* b */
  strcat( str1, "b 40218105156867728698573668525883168222119515413\n");
  /* beta */
  strcat( str1, "beta 115334401956802802075595682801335644058796914268\n");
  /* alpha0 */
  strcat( str1, "alpha0 191079354656274778837764015557338301375963168470\n");
  /* alpha1 */
  strcat( str1, "alpha1 71445317903696340296199556072836940741717506375");
  printf("%s\n", str1);

  if( pbc_param_init_set_str(c->params, str1) == 1) {printf("error during init... exiting now\n"); exit(1);}
  FILE *fp = fopen("params.txt", "w");
  if( fp == NULL ) {printf("error opening file... exiting now\n"); exit(1);}
  pbc_param_out_str( fp, c->params );

  fclose( fp );
}

//================================================
//================================================

void set_params_bak( struct curve_params* c ) {
  /* x  = - ( 2^62 + 2^55 + 1 )
     tx = 6*x^2  + 1 
     px = 36*x^4 + 36*x^3 + 24*x^2 + 6*x + 1
     nx = 36*x^4 + 36*x^3 + 18*x^2 + 6*x + 1
   */
  printf("--- calling set_params ---\n");
  const char * const foo = "4647714815446351873"; //( pow(2, 62) + pow(2, 55) + 1 );
  mpz_init( c->x );
  mpz_init( c->tx );
  mpz_init( c->nx );
  mpz_init( c->px );

  if( mpz_set_str( c->x, foo, 10 ) != 0 ) { printf("error... exiting now!\n"); exit(1); }
  mpz_neg( c->x, c->x );

  mpz_t dummy1, dummy2, dummy3;
  mpz_init( dummy1 );
  mpz_init( dummy2 );
  mpz_init( dummy3 );
  mpz_pow_ui( dummy1, c->x, 2 ); //dummy1=x^2
  mpz_pow_ui( dummy2, c->x, 3 ); //dummy2=x^3
  mpz_pow_ui( dummy3, c->x, 4 ); //dummy3=x^4
  
  mpz_mul_ui( c->tx, dummy1, 6 );
  mpz_add_ui( c->tx, c->tx,  1 );

  mpz_mul_ui( dummy2, dummy2, 36 ); //dummy2=36*x^3
  mpz_mul_ui( dummy3, dummy3, 36 ); //dummy3=36*x^4

  mpz_add( dummy3, dummy3, dummy2 ); //dummy3=36*x^4+36*x^3

  mpz_mul_ui( dummy2, c->x, 6   ); //dummy2=6*x
  mpz_add_ui( dummy2, dummy2, 1 ); //dummy2=6*x+1

  mpz_add( dummy3, dummy3, dummy2 ); //dummy3=36*x^4+36*x^3+6*x+1

  mpz_mul_ui( dummy2, dummy1, 24 ); //dummy2=24*x^2

  mpz_add( c->px, dummy3, dummy2 );

  mpz_mul_ui( dummy2, dummy1, 18 ); //dummy2=18*x^2

  mpz_add( c->nx, dummy3, dummy2 );
  
  gmp_printf("x  = %Zd\n", c->x);
  gmp_printf("tx = %Zd\n", c->tx);
  gmp_printf("px = %Zd\n", c->px);
  gmp_printf("nx = %Zd\n", c->nx);

  char str1[1024];
  char str2[1024];
  /* type of the ECC */
  strcat( str1, "type f\n"); 
  /* field characteristic */
  strcat( str1, "q ");
  mpz_get_str( str2, 10, c->px ); 
  strcat( str1, str2);
  strcat( str1, "\n");
  /* order of the curve */
  strcat( str1, "r ");
  mpz_get_str( str2, 10, c->nx );
  strcat( str1, str2);
  strcat( str1, "\n");
  /* b */
  strcat( str1, "b 2\n");
  /* beta */
  strcat( str1, "beta -2\n");
  /* alpha0 */
  strcat( str1, "alpha0 1\n");
  /* alpha1 */
  strcat( str1, "alpha1 1");
  printf("%s\n", str1);

  if( pbc_param_init_set_str(c->params, str1) == 1) {printf("error during init... exiting now\n"); exit(1);}
  FILE *fp = fopen("params.txt", "w");
  if( fp == NULL ) {printf("error opening file... exiting now\n"); exit(1);}
  pbc_param_out_str( fp, c->params );


  fclose( fp );
  mpz_clear( dummy1 );
  mpz_clear( dummy2 );
  mpz_clear( dummy3 );
}

//================================================
//================================================

