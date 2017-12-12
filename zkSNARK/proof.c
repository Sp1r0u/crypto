#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "proof.h"

void generate_proof (struct proof_t *this, struct curve_t *ec,
		     struct crs_t *crs, struct circuit_t *c, struct field_t *f) {
  printf ("\n**** Generating proof ****\n\n");
  
  int i;
  
  this->free = free_proof;
  this->init = init_proof;
  this->view = view_proof;

  this->genHPoly = generate_HPoly;
  
  this->init (this, ec);
  
  element_t tmpG1;
  element_t tmpG2;
  
  element_init_G1 (tmpG1, ec->pairing);
  element_init_G2 (tmpG2, ec->pairing);
  
  for (i=0; i<*(crs->nInt_V); i++) {
    gmp_printf(" V[%d]=%Zd\n",crs->EK->rvVP[i]->wire->id,crs->EK->rvVP[i]->wire->value->value);
    //a_i * r_v * v_i(s) * P
    element_mul_mpz (tmpG1, crs->EK->rvVP[i]->elt, crs->EK->rvVP[i]->wire->value->value);
    element_add (this->arvVP, this->arvVP, tmpG1);
    //a_i * alpha_v * r_v * v_i(s) * P
    element_mul_mpz (tmpG1, crs->EK->alpha_vrvVP[i]->elt, crs->EK->alpha_vrvVP[i]->wire->value->value);
    element_add (this->aalpha_vrvVP, this->aalpha_vrvVP, tmpG1);
  }

  for (i=0; i<*(crs->nInt_W); i++) {
    gmp_printf(" W[%d]=%Zd\n",crs->EK->rwWQ[i]->wire->id,crs->EK->rwWQ[i]->wire->value->value);
    //a_i * r_w * w_i(s) * Q
    element_mul_mpz (tmpG2, crs->EK->rwWQ[i]->elt, crs->EK->rwWQ[i]->wire->value->value);
    element_add (this->arwWQ, this->arwWQ, tmpG2);
    //a_i * alpha_w * r_w * w_i(s) * Q
    element_mul_mpz (tmpG2, crs->EK->alpha_wrwWQ[i]->elt, crs->EK->alpha_wrwWQ[i]->wire->value->value);
    element_add (this->aalpha_wrwWQ, this->aalpha_wrwWQ, tmpG2);
  }

  for (i=0; i<*(crs->nInt_Y); i++) {
    gmp_printf(" Y[%d]=%Zd\n",crs->EK->ryYP[i]->wire->id,crs->EK->ryYP[i]->wire->value->value);
    //a_i * r_v * v_i(s) * P
    element_mul_mpz (tmpG1, crs->EK->ryYP[i]->elt, crs->EK->ryYP[i]->wire->value->value);
    element_add (this->aryYP, this->aryYP, tmpG1);
    //a_i * alpha_v * r_v * v_i(s) * P
    element_mul_mpz (tmpG1, crs->EK->alpha_yryYP[i]->elt, crs->EK->alpha_yryYP[i]->wire->value->value);
    element_add (this->aalpha_yryYP, this->aalpha_yryYP, tmpG1);
  }

  for (i=0; i<*(crs->nInt_Y); i++) {
    element_mul_mpz (tmpG1, crs->EK->r_xbetaVP[i]->elt, crs->EK->r_xbetaVP[i]->wire->value->value);
    element_add (this->ar_xbetaVP, this->ar_xbetaVP, tmpG1);
  }

  this->genHPoly (this, c, f);
  
  this->view (this);
  
  element_clear (tmpG1);
  element_clear (tmpG2);
}

//==============================================
//==============================================

void generate_HPoly (struct proof_t *this, struct circuit_t *c, struct field_t *f) {
  printf (" inside generate_HPoly\n");

  this->genTPoly = generate_target_poly;
  this->genPPoly = generate_PPoly;

  struct poly_ring_t *T = malloc (sizeof (struct poly_ring_t)); // T(x): univariate polynomial
  struct poly_ring_t *P = malloc (sizeof (struct poly_ring_t)); // P(x): univariate polynomial

  this->genTPoly (c, f, T);
  print_poly (T);

  this->genPPoly ();
  print_poly (P);

  if (T->length>0) {
    free_poly_ring (T->head);
    free (T);
  }

  if (P->length>0) {
    free_poly_ring (P->head);
    free (P);
  }
  
  printf (" leaving generate_HPoly\n");
}

//==============================================
//==============================================

void generate_PPoly () {
  int i;

  struct poly_ring_t *V = malloc (sizeof (struct poly_ring_t)); // V(x)
  struct poly_ring_t *W = malloc (sizeof (struct poly_ring_t)); // W(x)
  struct poly_ring_t *Y = malloc (sizeof (struct poly_ring_t)); // Y(x)
  

  
  
}

//==============================================
//==============================================
  
void free_proof (struct proof_t *this) {
  
  element_clear (this->arvVP);
  element_clear (this->aalpha_vrvVP);
  
  element_clear (this->arwWQ);
  element_clear (this->aalpha_wrwWQ);
  
  element_clear (this->aryYP);
  element_clear (this->aalpha_yryYP);

  element_clear (this->ar_xbetaVP);

  if (this->H->length>0) {
    free_poly_ring (this->H->head);
    free (H);
  }
  
  free (this);
}
 
//==============================================
//==============================================
 
void init_proof (struct proof_t *this, struct curve_t *ec) {
  
  element_init_G1 (this->arvVP, ec->pairing);
  element_init_G1 (this->aalpha_vrvVP, ec->pairing);

  element_init_G2 (this->arwWQ, ec->pairing);
  element_init_G2 (this->aalpha_wrwWQ, ec->pairing);

  element_init_G1 (this->aryYP, ec->pairing);
  element_init_G1 (this->aalpha_yryYP, ec->pairing);

  element_init_G1 (this->ar_xbetaVP, ec->pairing);
  
  element_set0 (this->arvVP);
  element_set0 (this->aalpha_vrvVP);
  
  element_set0 (this->arwWQ);
  element_set0 (this->aalpha_wrwWQ);
  
  element_set0 (this->aryYP);
  element_set0 (this->aalpha_yryYP);

  element_set0 (this->ar_xbetaVP);
  
  H = malloc (sizeof (struct poly_ring_t)); // H(x): univariate polynomial

}
 
//==============================================
//==============================================

void view_proof (struct proof_t *this) {
  
  element_printf ("%B\n", this->arvVP);
  element_printf ("%B\n", this->aalpha_vrvVP);
  
  element_printf ("%B\n", this->arwWQ);
  element_printf ("%B\n", this->aalpha_wrwWQ);

  element_printf ("%B\n", this->aryYP);
  element_printf ("%B\n", this->aalpha_yryYP);

  element_printf ("%B\n", this->ar_xbetaVP);

}
