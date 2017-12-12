#include "proof.h"

//================================
//================================

void initProof (struct proof_t *this, struct crs_t *CRS) {

  this->view = viewProof;
  this->free = freeProof;

  this->verif = verifyProof;
  
  this->CRS = CRS;

  this->genP = genPPoly;
  this->genH = genHPoly;
  
  // evaluate the arith. circuit C
  this->CRS->C->eval (CRS->C);

  element_init_G1 (this->aRvVP, CRS->C->field->EC->pairing);
  element_init_G1 (this->aALPHAvRvVP, CRS->C->field->EC->pairing);

  element_init_G2 (this->aRwWQ, CRS->C->field->EC->pairing);
  element_init_G1 (this->aRwWP, CRS->C->field->EC->pairing);
  element_init_G2 (this->aALPHAwRwWQ, CRS->C->field->EC->pairing);
  element_init_G1 (this->aALPHAwRwWP, CRS->C->field->EC->pairing);
  
  element_init_G1 (this->aRyYP, CRS->C->field->EC->pairing);
  element_init_G1 (this->aALPHAyRyYP, CRS->C->field->EC->pairing);

  element_init_G1 (this->aRxBETAXP, CRS->C->field->EC->pairing);

  element_init_G2 (this->hSQ, CRS->C->field->EC->pairing);
  
  element_set0 (this->aRvVP);
  element_set0 (this->aALPHAvRvVP);
  
  element_set0 (this->aRwWQ);
  element_set0 (this->aRwWP);
  element_set0 (this->aALPHAwRwWQ);
  element_set0 (this->aALPHAwRwWP);
  
  element_set0 (this->aRyYP);
  element_set0 (this->aALPHAyRyYP);

  element_set0 (this->aRxBETAXP);

  element_set0 (this->hSQ);
  
  // temporary elements
  element_t tmp_elt1;
  element_t tmp_elt2;

  element_init_G1 (tmp_elt1, CRS->C->field->EC->pairing);
  element_init_G2 (tmp_elt2, CRS->C->field->EC->pairing);
  
  int i;
  for (i=0; i<*(CRS->nInt_V); i++) {
    // sum_i ( a_i * r_v * v_i(s) * P ) (i=Int wire IDs)
    element_set0 (tmp_elt1);
    element_mul_mpz (tmp_elt1, CRS->EK->RvVP[i]->elt, CRS->EK->RvVP[i]->wire->value->value);
    element_add (this->aRvVP, this->aRvVP, tmp_elt1);

    // sum_i ( a_i * alpha_v * r_v * v_i(s) * P ) (i=Int wire IDs)
    element_set0 (tmp_elt1);
    element_mul_mpz (tmp_elt1, CRS->EK->ALPHAvRvVP[i]->elt, CRS->EK->ALPHAvRvVP[i]->wire->value->value);
    element_add (this->aALPHAvRvVP, this->aALPHAvRvVP, tmp_elt1);
  }

  for (i=0; i<*(CRS->nInt_W); i++) {
    // sum_i ( a_i * r_w * w_i(s) * Q )
    element_set0 (tmp_elt2);
    element_mul_mpz (tmp_elt2, CRS->EK->RwWQ[i]->elt, CRS->EK->RwWQ[i]->wire->value->value);
    element_add (this->aRwWQ, this->aRwWQ, tmp_elt2);
    
    // sum_i ( a_i * alpha_w * r_w * w_i(s) * Q )
    element_set0 (tmp_elt2);
    element_mul_mpz (tmp_elt2, CRS->EK->ALPHAwRwWQ[i]->elt, CRS->EK->ALPHAwRwWQ[i]->wire->value->value);
    element_add (this->aALPHAwRwWQ, this->aALPHAwRwWQ, tmp_elt2);

    // sum_i ( a_i * alpha_w * r_w * w_i(s) * P )
    element_set0 (tmp_elt1);
    element_mul_mpz (tmp_elt1, CRS->EK->ALPHAwRwWP[i]->elt, CRS->EK->ALPHAwRwWP[i]->wire->value->value);
    element_add (this->aALPHAwRwWP, this->aALPHAwRwWP, tmp_elt1);
  }

  for (i=0; i<*(CRS->nInt_Y); i++) {
    // sum_i ( a_i * r_y * y_i(s) * P )
    element_set0 (tmp_elt1);
    element_mul_mpz (tmp_elt1, CRS->EK->RyYP[i]->elt, CRS->EK->RyYP[i]->wire->value->value);
    element_add (this->aRyYP, this->aRyYP, tmp_elt1);
    
    // sum_i ( a_i * alpha_y * r_y * y_i(s) * P )
    element_set0 (tmp_elt1);
    element_mul_mpz (tmp_elt1, CRS->EK->ALPHAyRyYP[i]->elt, CRS->EK->ALPHAyRyYP[i]->wire->value->value);
    element_add (this->aALPHAyRyYP, this->aALPHAyRyYP, tmp_elt1);
  }

  for (i=0; i<(CRS->C->nWire - CRS->C->nIOWire); i++) {
    // sum_i ( a_i * ( r_v * beta * v_i(s) + r_w * beta * w_i(s) + r_y * beta * y_i(s) ) * P )
    element_set0 (tmp_elt1);
    element_mul_mpz (tmp_elt1, CRS->EK->RxBETAXP[i]->elt, CRS->EK->RxBETAXP[i]->wire->value->value);
    element_add (this->aRxBETAXP, this->aRxBETAXP, tmp_elt1);
  }
  
  genHPoly (this);
  
  this->view (this);

  this->verif (this);
  
  element_clear (tmp_elt1);
  element_clear (tmp_elt2);
}

//================================
//================================

void viewProof (struct proof_t *this) {
  int i;
  printf ("\n");
  printf (" **** Evaluation of arithm. circuit ****\n");

  printf (" wire ID | type(input, intermed., output) | value(dec, hex)\n");
  for (i=0; i<this->CRS->C->nWire; i++) {
    this->CRS->C->wire[i]->view (this->CRS->C->wire[i]);
  }

  printf ("\n");
  printf (" **** Proof ****\n");
  element_printf(" aRvVP = %B\n", this->aRvVP);
  element_printf(" aRwWQ = %B\n", this->aRwWQ);
  element_printf(" aRyYP = %B\n", this->aRyYP);

  element_printf(" aALPHAvRvVP = %B\n", this->aALPHAvRvVP);
  element_printf(" aALPHAwRwWQ = %B\n", this->aALPHAwRwWQ);
  element_printf(" aALPHAwRwWP = %B\n", this->aALPHAwRwWP);
  element_printf(" aALPHAyRyYP = %B\n", this->aALPHAyRyYP);

  element_printf(" aRxBETAXP  = %B\n", this->aRxBETAXP);

  element_printf(" hSQ = %B\n", this->hSQ);
}

//================================
//================================

void freeProof (struct proof_t *this) {
  element_clear (this->aRvVP);
  element_clear (this->aRwWQ);
  element_clear (this->aRwWP);
  element_clear (this->aRyYP);
  
  element_clear (this->aALPHAvRvVP);
  element_clear (this->aALPHAwRwWQ);
  element_clear (this->aALPHAwRwWP);
  element_clear (this->aALPHAyRyYP);

  element_clear (this->aRxBETAXP);

  this->P->free (this->P);
  this->H->free (this->H);

  element_clear (this->hSQ);
  
  free (this);
}  

//================================
//================================

void genPPoly (struct proof_t *this) {
  struct wire_t *ptr;
 
  // temporary variables
  struct poly_t *tmp = malloc (sizeof (struct poly_t) );
  tmp->init = initPoly;
  tmp->init (tmp);
  
  struct poly_t *tmp_V = malloc (sizeof (struct poly_t) );
  tmp_V->init = initPoly;
  tmp_V->init (tmp_V);

  struct poly_t *tmp_W = malloc (sizeof (struct poly_t) );
  tmp_W->init = initPoly;
  tmp_W->init (tmp_W);

  struct poly_t *tmp_Y = malloc (sizeof (struct poly_t) );
  tmp_Y->init = initPoly;
  tmp_Y->init (tmp_Y);

  struct field_elt_t *tmp_elt = malloc (sizeof (struct field_elt_t) );
  tmp_elt->init = initFieldElt;
  tmp_elt->init (tmp_elt, this->CRS->C->field);
  mpz_set_ui (tmp_elt->value, 0); // tmp_elt <- 0
  
  struct poly_t *V = malloc (sizeof (struct poly_t) );
  V->init = initPoly;
  V->init (V);
  V->insertNode (V, tmp_elt->value, 0, tmp_elt->field); // V(x) = 0*x^0
  
  struct poly_t *W = malloc (sizeof (struct poly_t) );
  W->init = initPoly;
  W->init (W);
  W->insertNode (W, tmp_elt->value, 0, tmp_elt->field); // W(x) = 0*x^0
 
  struct poly_t *Y = malloc (sizeof (struct poly_t) );
  Y->init = initPoly;
  Y->init (Y);
  Y->insertNode (Y, tmp_elt->value, 0, tmp_elt->field); // Y(x) = 0*x^0
    
  int i;
  for (i=0; i<this->CRS->C->nWire; i++) {
    ptr = this->CRS->C->wire[i];

    if (ptr->V->length>0) {
      tmp_V->cstMulPoly (ptr->V, ptr->value, tmp_V); // tmp_V(x) <- a_i * V_i(x)
      tmp->copyPoly (V, tmp); // tmp(x) <- V(x)
      V->reset (V);
      V->addPoly (tmp_V, tmp, V); // V(x) += tmp_V(x)
      tmp_V->reset (tmp_V);
      tmp->reset (tmp);
    }
    
    if (ptr->W->length>0) {
      tmp_W->cstMulPoly (ptr->W, ptr->value, tmp_W); // tmp_W(x) <- a_i * W_i(x)
      tmp->copyPoly (W, tmp); // tmp(x) <- W(x)
      W->reset (W);
      W->addPoly (tmp_W, tmp, W); // W(x) += tmp_W(x)
      tmp_W->reset (tmp_W);
      tmp->reset (tmp);
    }

    if (ptr->Y->length>0) {
      tmp_Y->cstMulPoly (ptr->Y, ptr->value, tmp_Y); // tmp_Y(x) <- a_i * Y_i(x)
      tmp->copyPoly (Y, tmp); // tmp(x) <- Y(x)
      Y->reset (Y);
      Y->addPoly (tmp_Y, tmp, Y); // Y(x) += tmp_Y(x)
      tmp_Y->reset (tmp_Y);
      tmp->reset (tmp);
    }
    
  }

  tmp->mulPoly (V, W, tmp);

  this->P->subPoly (tmp, Y, this->P);
  
  tmp_elt->free (tmp_elt);

  tmp->free (tmp);
  
  tmp_V->free (tmp_V);
  tmp_W->free (tmp_W);
  tmp_Y->free (tmp_Y);
  
  V->free (V);
  W->free (W);
  Y->free (Y);

}

//================================
//================================

void genHPoly (struct proof_t *this) {

  this->P = malloc (sizeof (struct poly_t) );
  this->P->init = initPoly;
  this->P->init (this->P);

  this->H = malloc (sizeof (struct poly_t) );
  this->H->init = initPoly;
  this->H->init (this->H);

  genPPoly (this);

  this->H->divPoly (this->P, this->CRS->T, this->H);

  if (this->CRS->T->head->exponent < this->CRS->T->tail->exponent) {
    printf (" error in genHPoly... exiting now\n");
    exit (1);
  }					     

  // temporary elements
  element_t tmp_elt;
  element_init_G2 (tmp_elt, this->CRS->C->field->EC->pairing);

  //this->H->view (this->H);

  if (this->H->length!=0) {
    struct poly_node_t *node = this->H->head;
    while (node!=NULL) {
      element_mul_mpz (tmp_elt, this->CRS->EK->SQ[node->exponent]->elt, node->coefficient->value);
      element_add (this->hSQ, this->hSQ, tmp_elt);
      node = node->next;
    }
  }

  element_clear (tmp_elt);

}

//================================
//================================

void verifyProof (struct proof_t *this) {
  printf ("\n");
  printf (" --- Proof Verification ---\n");

  element_t R;
  element_t L;

  element_init_GT (R, this->CRS->C->field->EC->pairing);
  element_init_GT (L, this->CRS->C->field->EC->pairing);
  
  // EQ1
  element_pairing(R, this->aALPHAvRvVP, *(this->CRS->VK->Q));
  element_pairing(L, this->aRvVP,       this->CRS->VK->ALPHAvQ);
  
  if (element_cmp (R,L)==0) {
    printf (" EQ1 OK\n");
  }

  // EQ2
  element_pairing(R, this->aALPHAwRwWP, *(this->CRS->VK->Q));
  element_pairing(L, this->CRS->VK->ALPHAwP, this->aRwWQ);
  
  if (element_cmp (R,L)==0) {
    printf (" EQ2 OK\n");
  }

  // EQ3
  element_pairing(R, this->aALPHAyRyYP, *(this->CRS->VK->Q));
  element_pairing(L, this->aRyYP,       this->CRS->VK->ALPHAyQ);
  
  if (element_cmp (R,L)==0) {
    printf (" EQ3 OK\n");
  }

  // EQ4
  element_t tmp_elt1;
  element_t tmp_elt2;

  element_init_G1 (tmp_elt1, this->CRS->C->field->EC->pairing);
  element_init_G2 (tmp_elt2, this->CRS->C->field->EC->pairing);

  element_t R1, R2;
  element_t L1, L2;

  element_init_GT (R1, this->CRS->C->field->EC->pairing);
  element_init_GT (R2, this->CRS->C->field->EC->pairing);
  element_init_GT (L1, this->CRS->C->field->EC->pairing);
  element_init_GT (L2, this->CRS->C->field->EC->pairing);
  
  element_add (tmp_elt1, this->aRvVP, this->aRyYP); 
  element_pairing(R1, tmp_elt1, this->CRS->VK->BETAQ);

  element_pairing(R2, this->CRS->VK->BETAP, this->aRwWQ);

  element_mul (R, R1, R2); 

  element_pairing(L, this->aRxBETAXP, *(this->CRS->VK->Q));

  if (element_cmp (R,L)==0) {
    printf (" EQ4 OK\n");
  }
  
  // EQ5
  element_t tmp_elt3;
  element_t tmp_elt4;

  element_init_G1 (tmp_elt3, this->CRS->C->field->EC->pairing);
  element_init_G2 (tmp_elt4, this->CRS->C->field->EC->pairing);

  element_set0 (tmp_elt3);
  element_set0 (tmp_elt4);

  int i;
  // sum_i ( a_i * r_v * v_i(s) * P ) (i=IO+Int wire IDs)
  for (i=0; i<*(this->CRS->nIO_V); i++) {
    element_mul_mpz
      (tmp_elt1, this->CRS->VK->RvVP[i]->elt, this->CRS->VK->RvVP[i]->wire->value->value); // a_i * r_v * v_i(s) * P (i=IO wire IDs)
    element_add (tmp_elt3, tmp_elt3, tmp_elt1); // sum_i ( a_i * r_v * v_i(s) * P ) (i=IO wire IDs)
  }
  element_add (tmp_elt1, tmp_elt3, this->aRvVP);
  element_set0 (tmp_elt3);
  
  // sum_i ( a_i * r_w * w_i(s) * Q ) (i=IO+Int wire IDs)
  for (i=0; i<*(this->CRS->nIO_W); i++) {
    element_mul_mpz
      (tmp_elt2, this->CRS->VK->RwWQ[i]->elt, this->CRS->VK->RwWQ[i]->wire->value->value); // a_i * r_w * w_i(s) * Q (i=IO wire IDs)
    element_add (tmp_elt4, tmp_elt4, tmp_elt2); // sum_i ( a_i * r_w * w_i(s) * Q ) (i=IO wire IDs)
  }
  element_add (tmp_elt2, tmp_elt4, this->aRwWQ);
  element_set0 (tmp_elt4);
  
  // R=e(RvVP,RwWQ) 
  element_pairing(R, tmp_elt1, tmp_elt2);
  element_set0 (tmp_elt1);
  element_set0 (tmp_elt2);
  
  // sum_i ( a_i * r_y * y_i(s) * P ) (i=IO+Int wire IDs)
  for (i=0; i<*(this->CRS->nIO_Y); i++) {
    element_mul_mpz
      (tmp_elt1, this->CRS->VK->RyYP[i]->elt, this->CRS->VK->RyYP[i]->wire->value->value); // a_i * r_y * y_i(s) * P (i=IO wire IDs)
    element_add (tmp_elt3, tmp_elt3, tmp_elt1); // sum_i ( a_i * r_y * y_i(s) * P ) (i=IO wire IDs)
  }
  element_add (tmp_elt1, tmp_elt3, this->aRyYP);
  element_set0 (tmp_elt3);
  
  // L1=e(RyYP,Q) 
  element_pairing(L1, tmp_elt1, *(this->CRS->VK->Q));
  element_set0 (tmp_elt1);
  
  // L2=e(RyTP,hSQ) 
  element_pairing(L2, this->CRS->VK->RyTP, this->hSQ);

  element_mul (L, L1, L2); 

  if (element_cmp (R,L)==0) {
    printf (" EQ5 OK\n");
  }

  /*
  struct field_elt_t *y0 = malloc (sizeof (struct field_elt_t) );
  struct field_elt_t *y1 = malloc (sizeof (struct field_elt_t) );
  struct field_elt_t *y2 = malloc (sizeof (struct field_elt_t) );
  y0->init = initFieldElt;
  y1->init = initFieldElt;
  y2->init = initFieldElt;
  y0->init (y0, this->CRS->C->field);
  y1->init (y1, this->CRS->C->field);
  y2->init (y2, this->CRS->C->field);
 
  this->P->eval (this->P, this->CRS->param[3], y0); //y0 = P(s)
  this->H->eval (this->H, this->CRS->param[3], y1); //y1 = H(s)

  element_mul_mpz
      (tmp_elt4, *(this->CRS->VK->Q), y1->value);

  printf ("\n");
  element_printf ("\n H(s)*Q:%B\n",tmp_elt4);
  element_printf ("\n hSQ:%B\n",this->hSQ);
  printf ("\n");
  
  this->CRS->T->eval (this->CRS->T, this->CRS->param[3], y2); //y2 = T(s)

  struct field_elt_t *y3 = malloc (sizeof (struct field_elt_t) );
  y3->init = initFieldElt;
  y3->init (y3, this->CRS->C->field);

  y3->field->mulElt (y1, y2, y3); //y3 = H(s) * T(s)

  gmp_printf (" P(s)      = %Zd\n", y0);
  gmp_printf (" T(s)*H(s) = %Zd\n", y3);

  y1->field->mulElt (y3, this->CRS->param[2], y1); //y1 = ry * H(s) * T(s)
  
  element_mul_mpz
      (tmp_elt1, *(this->CRS->VK->P), y1->value);

  element_pairing(L, tmp_elt1, *(this->CRS->VK->Q));
  element_printf("\n\n %B\n",L);
  
  element_printf("\n\n %B\n",L2);
  exit(0);
  



  */
  
  element_clear (R);
  element_clear (L);

  element_clear (tmp_elt1);
  element_clear (tmp_elt2);

  element_clear (R1);
  element_clear (R2);
  
  element_clear (L1);
  element_clear (L2);

  element_clear (tmp_elt3);
  element_clear (tmp_elt4);
}

//================================
//================================

