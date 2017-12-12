#include <time.h>
#include <gmp.h>
#include <stdbool.h>

#include "EC.h"
#include "field.h"
#include "poly.h"
#include "circuit.h"
#include "CRS.h"
#include "proof.h"

// global var. declaration
bool _debug=true;
unsigned int _order=101;

//================================
//================================

int main (int argc, char **argv) {
  
  if (argc!=3) {
    printf (" Please specify 1 or 0 for debug %i\n", argc);
    exit (1);
  }
    
  FILE *fptr;
  fptr = fopen( argv[1], "r" );
  if( fptr==NULL ) {
    printf("file %s does not exist... exiting now!\n", argv[1]);
    exit(1);
  }
  
  _debug = (atoi(argv[2]) == 0) ? false : true;

  printf("\n");
  printf (" **** To Do List ****\n");
  printf (" 1. implementing polynomial operations using FFT\n");
  printf (" 2. rewriting BN.c\n");
  printf (" 3. rewriting CRS.c\n");
  printf (" 4. rewriting simplify/sorting in poly.c\n");
  
  printf("\n");
  printf(" **** Debug Mode ****\n");
  printf(" set to %s\n", _debug?"true":"false");
  
  gmp_randstate_t RND_STATE;

  // initialize the random state with default algorithm
  gmp_randinit (RND_STATE, 0, 128);

  // seed the state with an unsigned long int
  long seed;
  time (&seed);
  gmp_randseed_ui (RND_STATE, seed);

  // set the elliptic curve
  struct curve_t *EC = malloc (sizeof (struct curve_t) );
  EC->init = initEC;
  EC->init (EC);

  // set the field
  struct field_t *Fp = malloc (sizeof (struct field_t) );
  Fp->init = initField;
  Fp->init (Fp, EC);

  // set the arithmetic circuit
  struct circuit_t *C = malloc (sizeof (struct circuit_t) );
  C->init = initCircuit;
  C->init (C, Fp, fptr, RND_STATE);

  // generates CRS
  struct crs_t *CRS = malloc (sizeof (struct crs_t) );
  CRS->init = initCRS;
  CRS->init (CRS, C, RND_STATE);

  // evaluates the circuit & generates the proof
  struct proof_t *PI = malloc (sizeof (struct proof_t) );
  PI->init = initProof;
  PI->init (PI, CRS);
  
  // free memory
  gmp_randclear (RND_STATE);
  PI->free (PI);
  CRS->free (CRS); 
  C->free (C);
  Fp->free (Fp);
  EC->free (EC);
  
  return 0;
}
