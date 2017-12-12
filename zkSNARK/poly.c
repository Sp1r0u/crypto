#include "poly.h"

/* evaluate Lagrange poly. at x */
/*void evaluate_Lagrange_polynomial( struct circuit_t *c, mpz_t x ) {
  gmp_printf("\n   Evaluating Lagrange poly. at %Zd\n", x);

  int i, j;
  mpz_t num; 
  mpz_t den;

  mpz_init_set_ui( num, 1 );
  mpz_init_set_ui( den, 1 );
  
  for (i=0; i<c->nMul_gate; i++) {  
    for (j=0; j<c->nMul_gate; j++) {
      if (i!=j) {
	num *= x - c->label[j];
	den *= c->label[i] - c->label[j]  
      }
    }
    mpz_init_set_ui( num, 1 );
    mpz_init_set_ui( den, 1 );
  }
}
*/
