#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gmp.h>

#include "field.h"
#include "element.h"
#include "toolbox.h"

int main( int argc, char **argv ) {

  /* characteristic of the field */
  char str[] = "313";

  gmp_randstate_t state;

  /* initialize the random state with default algorithm */
  gmp_randinit_default( state );

  /* seed the state with an unsigned long int */
  gmp_randseed_ui( state, 1234567890 );
  
  /* Fp: base field */
  struct field_t Fp;
  if( initialize( &Fp, str, 1, NULL ) != true ) {
    printf("error during field initialization... exiting now!\n");
    exit(0);
  }

  /* p: monic polynomial */
  struct element_t *p = malloc( sizeof( struct element_t ));
  generate_random_monic_polynomial( &Fp, p, 3, state );
  print_element( p );

  struct element_t *q = malloc( sizeof( struct element_t ));
  generate_random_monic_polynomial( &Fp, q, 4, state );
  print_element( q );

  struct element_t *r = malloc( sizeof( struct element_t ));
  multiply_elements( p, q, r );
  
  //is_poly_irreducible( p );
  
  /* Fp2: quadratic extension of Fp */
  struct field_t Fp2;
  if( initialize( &Fp2, NULL, 2, &Fp ) != true ) {
    printf("error during field initialization... exiting now!\n");
    exit(0);
  }

  /* freeing memory */
  gmp_randclear( state );
  free( p );
  
  return 0;
}
