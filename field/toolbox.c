#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <gmp.h>

#include "field.h"
#include "toolbox.h"
#include "element.h"

//==================================================
//==================================================

bool is_poly_irreducible( struct element_t *p ) {
  
  mpz_t coeff;
  mpz_init( coeff );
  mpz_set_ui( coeff, 1 );
  
  struct element_t *u = malloc( sizeof( struct element_t ));
  set_element( u, coeff, 1 );
  print_element( u );
  
  mpz_clear( coeff );
  free( u );
  
  return false;
}

//==================================================
//==================================================

void str2mpz( char *b, mpz_t m ) {
  mpz_set_str( m, b, 10 );
}

//==================================================
//==================================================

void mpz2str( mpz_t m, char *b ) {
  mpz_get_str( b, 10, m );
}

//==================================================
//==================================================

void generate_random_monic_polynomial( struct field_t *f, struct element_t *p, unsigned int n, gmp_randstate_t s ) {
  
  /* bound for mpz_urandomm */
  mpz_t upBnd;
  mpz_init_set( upBnd, f->n );

  mpz_t coeff;
  mpz_init( coeff );

  p->length = 0;
  
  int i;
  for( i=0; i<=n; i++ ) {

    if( i==n ) {mpz_set_ui( coeff, 1 );}
    else {mpz_urandomm( coeff, s, upBnd );} 

    set_element( p, coeff, i );
    
  }

  p->field = f;
  
  mpz_clear( upBnd );
  mpz_clear( coeff );
}
