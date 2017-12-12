#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "BN.h"

void init_curve( struct curve_t *c ) {
  printf ("\n**** EC Settings ****\n");
  mpz_init( c->x );
  mpz_init( c->tx );
  mpz_init( c->nx );
  mpz_init( c->px );
  
  mpz_set_str( c->x, "1", 10 );
  mpz_neg( c->x, c->x );

  mpz_set_str( c->tx, "1", 10 );
  mpz_neg( c->tx, c->tx );
  
  mpz_set_str( c->px, "205523667896953300194896352429254920972540065223", 10 );
  
  mpz_set_str( c->nx, "205523667896953300194895899082072403858390252929", 10 );
  
  gmp_printf("\nx  = %Zd\n", c->x);
  gmp_printf("tx = %Zd\n", c->tx);
  gmp_printf("px = %Zd\n", c->px);
  gmp_printf("nx = %Zd\n", c->nx);
  
  char str1[1024] = "";
  char str2[1024] = "";
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

  pbc_param_init_set_str(c->pbc_params, str1);

  FILE *fp = fopen("pbc_params.txt", "w");

  if (fp==NULL) {
    printf("error while opening file... exiting now\n");
    exit(1);
  }
    
  pbc_param_out_str( fp, c->pbc_params );

  pairing_init_pbc_param( c->pairing, c->pbc_params );

  if( pairing_is_symmetric( c->pairing ) == 0 ) {
    printf("\nSanety check: asymmetric pairing function\n");
  }
  else {
    printf("\nSanety check: symmetric pairing function..., exiting now!\n");
    exit(1);
  }

  printf("length in bytes to represent an element of G1: %d\n", pairing_length_in_bytes_G1(c->pairing));
  printf("length in bytes to represent an element of G2: %d\n", pairing_length_in_bytes_G2(c->pairing));
  
  element_init_G1( c->P, c->pairing );
  element_init_G2( c->Q, c->pairing );

  element_random( c->P );
  element_random( c->Q );

  element_printf("P = %B\n", c->P);
  element_printf("Q = %B\n", c->Q);

  // read pbc_pairing.h L18
  mpz_init( c->r );
  mpz_set( c->r, c->pairing->r );
  gmp_printf("r = %Zd\n", c->r);
  
  // call mpz_probab_prime_p to check whether r is prime
  // reps=2 if n is definitely prime 
  // reps=1 if n is probably prime (without being certain)
  // reps=0 if n is definitely composite
  
  int dummy = -99;
  dummy = mpz_probab_prime_p( c->r, dummy );
  if(dummy==2 ) {
    printf("Sanety check: G1, G2, GT order are prime\n");
  }
  else if(dummy==1 ) {
    printf("Sanety check: G1, G2, GT order are PROBABLY prime\n");
  }
  else {
    printf("Sanety check: G1, G2, GT order are NOT prime,... exiting now\n");
    exit(1);
  }
  
  fclose( fp );

}

//==============================================
//==============================================

void free_curve( struct curve_t *c ) {

  mpz_clear (c->x);
  mpz_clear (c->tx);
  mpz_clear (c->nx);
  mpz_clear (c->px);
  mpz_clear (c->r);

  element_clear (c->P);
  element_clear (c->Q);

  pbc_param_clear (c->pbc_params);

  pairing_clear (c->pairing);  
    
  free( c );

}
