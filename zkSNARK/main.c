#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gmp.h>

#include "/home/elmagnifico/pbc-0.5.14/include/pbc.h"
#include "field.h"
#include "BN.h"
//#include "element.h"
#include "circuit.h"
//#include "wire.h"
//#include "gate.h"

int main( int argc, char **argv ) {

  FILE *fptr;
  fptr = fopen( argv[1], "r" );
  if( fptr==NULL ) {
    printf("file %s does not exist... exiting now!\n", argv[1]);
    exit(1);
  }
  
  gmp_randstate_t RND_STATE;
  
  /* initialize the random state with default algorithm */
  //gmp_randinit_default( state );
  gmp_randinit( RND_STATE, 0, 128 );

  /* seed the state with an unsigned long int */
  long seed;
  time( &seed );
  //gmp_randseed_ui( state, 1234567890 );
  gmp_randseed_ui( RND_STATE, seed );
  
  struct curve_t *BN = malloc( sizeof( struct curve_t ));

  init_curve( BN );

  struct field_t *Fp = malloc( sizeof( struct field_t ));

  init_field( Fp, BN );
  
  struct circuit_t *C = malloc( sizeof( struct circuit_t ));
  
  init_circuit( C, Fp, fptr, RND_STATE );

  //display_circuit( C );
  //display_mul_gate( C );

  set_random_input_values( C, Fp, RND_STATE );
  display_circuit( C );
  display_mul_gate( C );

  eval_circuit( C );
  display_circuit( C );
  /*  struct element_t *e1 = malloc( sizeof( struct element_t ));
  struct element_t *e2 = malloc( sizeof( struct element_t ));
  struct element_t *e3 = malloc( sizeof( struct element_t ));
  struct element_t *e4 = malloc( sizeof( struct element_t ));
  struct element_t *e5 = malloc( sizeof( struct element_t ));
  struct element_t *e6 = malloc( sizeof( struct element_t ));

  init_elt_rnd( e1, Fp, RND_STATE );
  init_elt_ui( e2, Fp, 1234567890 );
  mul_elt( e1, e2, e3 );
  div_elt( e1, e2, e4 );
  add_elt( e1, e2, e5 );
  sub_elt( e1, e2, e6 );
  //printf("%p %p\n", e1->field, Fp);
  
  gmp_printf("       e1: %Zd\n", e1->value);
  gmp_printf("       e2: %Zd\n", e2->value);
  gmp_printf(" e3=e1*e2: %Zd\n", e3->value);
  gmp_printf(" e4=e1/e2: %Zd\n", e4->value);
  gmp_printf(" e5=e1+e2: %Zd\n", e5->value);
  gmp_printf(" e6=e1-e2: %Zd\n", e6->value);
  */

  
  /* freeing memory */
  free( BN );
  free( Fp );
  gmp_randclear( RND_STATE );

  fclose( fptr );
  
  return 0;
}
