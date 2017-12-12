#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <stdbool.h>

#include "/home/elmagnifico/pbc-0.5.14/include/pbc.h"
#include "field.h"
#include "BN.h"
#include "circuit.h"
#include "CRS.h"
#include "proof.h"

/* global var. declaration */
bool fDebug=false; //true;

unsigned int fOrder=101;

//==============================================
//==============================================

int main( int argc, char **argv ) {

  //char *str = malloc( 4 * sizeof( struct curve_t ));
  //strcpy (str, "papa"); 
  //printf ("%s\n", str);
  //free (str);

  //int *tmp;
  //tmp = NULL;
  //tmp = malloc( 10 * sizeof( int ));

  //int i;
  //for (i=0; i<10; i++) {
  //tmp[i] = i;
  //}
  //for (i=0; i<10; i++) {
  //printf(" %d\n",tmp[i]);
  //}
  
  //exit(0);
  
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

  display_circuit( C );

  struct crs_t *CRS = malloc( sizeof( struct crs_t ));
  
  generate_crs( CRS, C, BN, Fp, RND_STATE );

  printf("outside generate_crs\n");
  
  set_random_input_values( C, RND_STATE );

  display_circuit( C );
  //display_mul_gate( C );

  eval_circuit( C );

  display_circuit( C );

  struct proof_t *PI = malloc( sizeof( struct proof_t ));
  PI->genPI = generate_proof;
  PI->genPI (PI, BN, CRS, C, Fp);
  
  printf("--- testing here ---\n");

  /*int i;
  //unsigned int hi_exp;
  //for (i=0; i<C->nMul_gate; i++) {
  /*struct poly_ring_t *irr_poly = malloc( sizeof( struct poly_ring_t ));
    irr_poly->length=0;
    printf(" poly: ");
    print_poly( C->gate[C->index[i]]->lagrange_poly );
    hi_exp = get_largest_exp( C->gate[C->index[i]]->lagrange_poly );
    printf(" highest exp: %d\n", hi_exp);
    random_monic_poly( Fp, irr_poly, hi_exp, RND_STATE );
    printf(" irr. poly: ");
    print_poly( irr_poly );
    printf("\n");
    if (is_poly_irreducible( Fp, irr_poly )==true ) {printf("irreducible\n");}*/
  
  struct poly_ring_t *p1 = malloc( sizeof( struct poly_ring_t ));
  p1->length=0;
  
  struct poly_ring_t *p2 = malloc( sizeof( struct poly_ring_t ));
  p2->length=0;

  random_poly( Fp, p1, 2, RND_STATE );
  random_poly( Fp, p2, 2, RND_STATE );

  printf(" p1:");
  print_poly( p1 );

  printf(" p2:");
  print_poly( p2 );

  struct poly_ring_t *p3 = malloc( sizeof( struct poly_ring_t ));
  p3->length=0;
 
  //mul_poly( p1,  p2, p3);
  //print_poly( p3 );
  //poly_ordering( p3 );
  div_poly( p1,  p2, p3 );
  printf(" p3:");
  print_poly( p3 );
  //}
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
  printf ("\n**** Freeing Memory ****\n\n");
  printf ("--- freeing proof ---\n");
  PI->free (PI);
  printf ("--- freeing CRS ---\n");
  free_crs( CRS, C );
  printf ("--- freeing circuit ---\n");
  free_circuit( C );
  printf ("--- freeing field ---\n");
  free_field( Fp );
  printf ("--- freeing curve ---\n");
  free_curve( BN );
  //free_circuit( C );
  //free_crs( CRS );
  
  gmp_randclear( RND_STATE );

  fclose( fptr );
  
  return 0;
}
