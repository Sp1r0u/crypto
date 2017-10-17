#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <gmp.h>

#include "field.h"

void init_field( struct field_t *f, struct curve_t *c ) {
  mpz_init( f->p );
  mpz_set( f->p, c->r );
}

void rnd_field_elt( mpz_t elt, struct field_t *f, gmp_randstate_t rnd_state ) {
  mpz_init( elt );
  while ( mpz_cmp_ui(elt, 0)==0 ) {
    mpz_urandomm( elt, rnd_state, f->p );
  }
}
