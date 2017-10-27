#include <string.h>

#include "CRS.h"

//==============================================
//==============================================

void generate_VKey( struct VKey_t *vk, struct circuit_t *c, struct curve_t *ec, struct field_t *f, mpz_t param[] ) {
  int i;
  mpz_t target; //target polynomial evaluated at the secret point s
  mpz_init( target );
  evaluate_target_poly( c, param[3], target, f );
  gmp_printf("target: %Zd\n", target);

  vk->rvVP = malloc( ( c->nIOWire ) * sizeof( element_t ));
  vk->rwWQ = malloc( ( c->nIOWire ) * sizeof( element_t ));
  vk->ryYP = malloc( ( c->nIOWire ) * sizeof( element_t ));
   
  vk->P = &(ec->P);
  vk->Q = &(ec->Q);
  
  element_init_G2( vk->alpha_vQ, ec->pairing );
  element_init_G2( vk->alpha_wQ, ec->pairing );
  element_init_G1( vk->alpha_wP, ec->pairing );
  element_init_G1( vk->betaP,    ec->pairing );
  element_init_G2( vk->betaQ,    ec->pairing );
  element_init_G1( vk->rytP,     ec->pairing );
  for (i=0; i<c->nIOWire; i++) {
    element_init_G1( vk->rvVP[i], ec->pairing );
    element_init_G2( vk->rwWQ[i], ec->pairing );
    element_init_G1( vk->ryYP[i], ec->pairing );
  }
    
  element_mul_mpz( vk->alpha_vQ, *(vk->Q), param[4] );
  element_mul_mpz( vk->alpha_wQ, *(vk->Q), param[5] );
  element_mul_mpz( vk->alpha_wP, *(vk->P), param[5] );
  element_mul_mpz( vk->betaP,    *(vk->P), param[7] );
  element_mul_mpz( vk->betaQ,    *(vk->Q), param[7] );
  element_mul_mpz( vk->rytP,     *(vk->P), param[2] ); //rytp <- ry   * P
  element_mul_mpz( vk->rytP,     vk->rytP, target   ); //rytp <- t(s) * rytp

  //  for (i=0; i<c->nWire; i++) {
  // if ( strcmp(c->wire[i].type, "input")==0 || strcmp(c->wire[i].type, "output")==0 ) {
      
  // }
  //}
  
  printf("in generate_VKey\n");
  element_printf("P = %B\n", vk->P);
  element_printf("Q = %B\n", vk->Q); 
  element_printf("alpha_vQ = %B\n", vk->alpha_vQ);
  element_printf("alpha_wQ = %B\n", vk->alpha_wQ);
  element_printf("alpha_wP = %B\n", vk->alpha_wP);
  element_printf("betaP    = %B\n", vk->betaP);
  element_printf("betaQ    = %B\n", vk->betaQ);
  element_printf("rytP     = %B\n", vk->rytP);
  exit(0);
}

//==============================================
//==============================================

void generate_EKey( struct EKey_t *ek, struct circuit_t *c ) {
  
}
  
//==============================================
//==============================================

void generate_crs( struct crs_t *crs, struct circuit_t *c, struct curve_t *ec, struct field_t *f, gmp_randstate_t rnd_state ) {
  int i;
  printf ("\n**** Generating CRS ****\n\n");
    
  for (i=0; i<NPARAMS; i++) {
    mpz_init( crs->param[i] );
    if (i!=2) {rnd_field_elt( crs->param[i], f, rnd_state);}
    else {
      mpz_mul( crs->param[2], crs->param[0], crs->param[1] ); //ry=rv*rw [mod p]
      mpz_mod( crs->param[2], crs->param[2], f->p );
    }
  }

  for (i=0; i<NPARAMS; i++) {gmp_printf("params[%d]: %Zd\n", i, crs->param[i]);}
  
  //exit(0);
  
  crs->EK = malloc( sizeof( struct EKey_t ));
  crs->VK = malloc( sizeof( struct VKey_t ));
  
  generate_EKey( crs->EK, c );
  generate_VKey( crs->VK, c, ec, f, crs->param );
}

//==============================================
//==============================================

void evaluate_target_poly( struct circuit_t *c, mpz_t x, mpz_t t, struct field_t *f) {

  int i;
  mpz_set_ui( t, 1 );
 
  mpz_t foo;
  mpz_init( foo );

  for (i=0; i<c->nMul_gate; i++) {
    mpz_sub( foo, x, c->gate[c->index[i]]->label );//foo=(x-x_g) [mod p]
    mpz_mod( foo, foo, f->p );
    mpz_mul( t, t, foo );//t=t*foo [mod p]
    mpz_mod( t, t, f->p );
  }
}

//==============================================
//==============================================
