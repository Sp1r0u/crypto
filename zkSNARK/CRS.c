#include "CRS.h"

void initCRS (struct crs_t *this, struct circuit_t *C, gmp_randstate_t state) {

  this->view = viewCRS;
  this->free = freeCRS;

  this->Lagrange = generateLagrangePoly;

  this->genQAP  = generateQAP;
  this->genQAPP = generateQAPPoly;

  this->genEKVK = generateEKeyVKey;
  
  this->genT = generateTargetPoly;

  this->genRxXY = generateRxXY;
  
  this->param = malloc (NPARAMS * sizeof (struct field_elt_t*) );

  this->EK = malloc( sizeof( struct EKey_t ));
  this->VK = malloc( sizeof( struct VKey_t ));
    
  this->C = C;
  
  int i;
  for (i=0; i<NPARAMS; i++) {
    this->param[i] = malloc (sizeof (struct field_elt_t) );
    this->param[i]->init = initFieldElt;
    this->param[i]->init (this->param[i], this->C->field);

    //mpz_set_ui (this->param[i]->value, 1);

    if (i!=2) {
      this->param[i]->field->rndElt (this->param[i], state);
    }
    else {
      //Ry = Rv * Rw [mod r]
      this->param[i]->field->mulElt (this->param[0], this->param[1], this->param[2]); 
    }
    
  }

  this->Lagrange (this);

  this->nIO_V = malloc( sizeof( unsigned int ));
  this->nIO_W = malloc( sizeof( unsigned int ));
  this->nIO_Y = malloc( sizeof( unsigned int ));

  this->nInt_V = malloc( sizeof( unsigned int ));
  this->nInt_W = malloc( sizeof( unsigned int ));
  this->nInt_Y = malloc( sizeof( unsigned int ));
  
  *(this->nIO_V) = 0;
  *(this->nIO_W) = 0;
  *(this->nIO_Y) = 0;

  *(this->nInt_V) = 0;
  *(this->nInt_W) = 0;
  *(this->nInt_Y) = 0;

  this->genQAP (this);

  this->genEKVK (this);
  
  this->view (this);
}

//================================
//================================

void generateEKeyVKey (struct crs_t *this) {

  // T: target poly
  this->T = malloc (sizeof (struct poly_t) ); //univariate polynomial  
  this->T->init = initPoly;
  this->T->init (this->T);

  this->genT (this->T, this->C);

  struct field_elt_t *tmp_elt = malloc (sizeof (struct field_elt_t) );
  tmp_elt->init = initFieldElt;
  tmp_elt->init (tmp_elt, this->C->field);
  
  this->T->eval (this->T, this->param[3], tmp_elt);
  //tmp_elt->view (tmp_elt);

  this->VK->RvVP = malloc( *( this->nIO_V ) * sizeof( struct extraKey_t* ));
  this->VK->RwWQ = malloc( *( this->nIO_W ) * sizeof( struct extraKey_t* ));
  this->VK->RyYP = malloc( *( this->nIO_Y ) * sizeof( struct extraKey_t* ));

  this->EK->RvVP = malloc( *( this->nInt_V ) * sizeof( struct extraKey_t* ));
  this->EK->RwWQ = malloc( *( this->nInt_W ) * sizeof( struct extraKey_t* ));
  this->EK->RyYP = malloc( *( this->nInt_Y ) * sizeof( struct extraKey_t* ));

  this->EK->ALPHAvRvVP = malloc( *( this->nInt_V ) * sizeof( struct extraKey_t* ));
  this->EK->ALPHAwRwWQ = malloc( *( this->nInt_W ) * sizeof( struct extraKey_t* ));
  this->EK->ALPHAwRwWP = malloc( *( this->nInt_W ) * sizeof( struct extraKey_t* ));
  this->EK->ALPHAyRyYP = malloc( *( this->nInt_Y ) * sizeof( struct extraKey_t* ));

  this->EK->RxBETAXP =
    malloc( ( this->C->nWire - this->C->nIOWire ) * sizeof( struct extraKey_t* ));
  
  this->VK->P = &(this->C->field->EC->P);
  this->VK->Q = &(this->C->field->EC->Q);

  element_init_G2( this->VK->ALPHAvQ, this->C->field->EC->pairing );
  element_init_G2( this->VK->ALPHAwQ, this->C->field->EC->pairing );
  element_init_G1( this->VK->ALPHAwP, this->C->field->EC->pairing );
  element_init_G2( this->VK->ALPHAyQ, this->C->field->EC->pairing );
  element_init_G1( this->VK->BETAP,   this->C->field->EC->pairing );
  element_init_G2( this->VK->BETAQ,   this->C->field->EC->pairing );
  element_init_G1( this->VK->RyTP,    this->C->field->EC->pairing );
  
  int i;
  
  for (i=0; i<*(this->nIO_V); i++) {
    this->VK->RvVP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->VK->RvVP[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->VK->RvVP[i]->value);
  }
 
  for (i=0; i<*(this->nIO_W); i++) {
    this->VK->RwWQ[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G2 (this->VK->RwWQ[i]->elt, this->C->field->EC->pairing );
    mpz_init (this->VK->RwWQ[i]->value);
  }

  for (i=0; i<*(this->nIO_Y); i++) {
    this->VK->RyYP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->VK->RyYP[i]->elt, this->C->field->EC->pairing );
    mpz_init (this->VK->RyYP[i]->value);
  }

  for (i=0; i<*(this->nInt_V); i++) {
    this->EK->RvVP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->EK->RvVP[i]->elt, this->C->field->EC->pairing );
    mpz_init (this->EK->RvVP[i]->value);
    
    this->EK->ALPHAvRvVP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->EK->ALPHAvRvVP[i]->elt, this->C->field->EC->pairing );
    mpz_init (this->EK->ALPHAvRvVP[i]->value);
  }

  for (i=0; i<*(this->nInt_W); i++) {
    this->EK->RwWQ[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G2 (this->EK->RwWQ[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->RwWQ[i]->value);
    
    this->EK->ALPHAwRwWQ[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G2 (this->EK->ALPHAwRwWQ[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->ALPHAwRwWQ[i]->value);
    
    this->EK->ALPHAwRwWP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->EK->ALPHAwRwWP[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->ALPHAwRwWP[i]->value);
  }
  
  for (i=0; i<*(this->nInt_Y); i++) {
    this->EK->RyYP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->EK->RyYP[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->RyYP[i]->value);
    
    this->EK->ALPHAyRyYP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->EK->ALPHAyRyYP[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->ALPHAyRyYP[i]->value);
  }

  for (i=0; i<(this->C->nWire - this->C->nIOWire); i++) {
    this->EK->RxBETAXP[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G1 (this->EK->RxBETAXP[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->RxBETAXP[i]->value);
  }

  element_mul_mpz( this->VK->ALPHAvQ, *(this->VK->Q), this->param[4]->value );
  element_mul_mpz( this->VK->ALPHAwQ, *(this->VK->Q), this->param[5]->value );
  element_mul_mpz( this->VK->ALPHAwP, *(this->VK->P), this->param[5]->value );
  element_mul_mpz( this->VK->ALPHAyQ, *(this->VK->Q), this->param[6]->value );
  element_mul_mpz( this->VK->BETAP,   *(this->VK->P), this->param[7]->value );
  element_mul_mpz( this->VK->BETAQ,   *(this->VK->Q), this->param[7]->value );
  element_mul_mpz( this->VK->RyTP,    *(this->VK->P), this->param[2]->value );
  element_mul_mpz( this->VK->RyTP,     this->VK->RyTP, tmp_elt->value       );

  // generate RxXY and ALPHAxRxXY (x={v,w,y} y={P,Q})
  this->genRxXY (this);

  // generate s^j * Q (j=0,..,deg(t))
  this->EK->SQ = malloc( (this->T->hiExp (this->T) + 1) * sizeof( struct extraKey_t* ));
  for (i=0; i<=(int)(this->T->hiExp (this->T)); i++) {
    this->EK->SQ[i] = malloc( sizeof( struct extraKey_t ));
    element_init_G2 (this->EK->SQ[i]->elt, this->C->field->EC->pairing);
    mpz_init (this->EK->SQ[i]->value);
    tmp_elt->field->powuiElt (this->param[3], tmp_elt, (unsigned int)(i));
    element_mul_mpz (this->EK->SQ[i]->elt, *(this->VK->Q), tmp_elt->value);
    mpz_set (this->EK->SQ[i]->value, tmp_elt->value);
  }
  
  tmp_elt->free (tmp_elt);
}

//================================
//================================

void generateRxXY (struct crs_t *this) {

  unsigned int iIOV = 0;
  unsigned int iIOW = 0;
  unsigned int iIOY = 0;
  
  unsigned int iIntV = 0;
  unsigned int iIntW = 0;
  unsigned int iIntY = 0;

  unsigned int iInt = 0;

  struct field_elt_t *y0 = malloc( sizeof ( struct field_elt_t ) );
  y0->init = initFieldElt;
  y0->init (y0, this->C->field);

  struct field_elt_t *y1 = malloc( sizeof ( struct field_elt_t ) );
  y1->init = initFieldElt;
  y1->init (y1, this->C->field);

  struct field_elt_t *y2 = malloc( sizeof ( struct field_elt_t ) );
  y2->init = initFieldElt;
  y2->init (y2, this->C->field);

  struct field_elt_t *y3 = malloc( sizeof ( struct field_elt_t ) );
  y3->init = initFieldElt;
  y3->init (y3, this->C->field);

  int i;
  struct wire_t *ptr;
  
  for (i=0; i<this->C->nWire; i++) {

    ptr = this->C->wire[i];

    if ( strcmp(ptr->type,"input")==0 || strcmp(ptr->type,"output")==0 ) {

      if (ptr->left_ctrb->length>0) {
	this->VK->RvVP[iIOV]->wire = ptr;
	ptr->V->eval (ptr->V, this->param[3], y0);  // y0 = V(s)
	y0->field->mulElt (this->param[0], y0, y0); // y0 = Rv * V(s)
	element_mul_mpz (this->VK->RvVP[iIOV]->elt, *(this->VK->P), y0->value);
	mpz_set (this->VK->RvVP[iIOV]->value, y0->value);
	iIOV+=1;
      }
 
      if (ptr->right_ctrb->length>0) {
	this->VK->RwWQ[iIOW]->wire = ptr;
	ptr->W->eval (ptr->W, this->param[3], y0);  // y0 = W(s)
	y0->field->mulElt (this->param[1], y0, y0); // y0 = Rw * W(s)
	element_mul_mpz (this->VK->RwWQ[iIOW]->elt, *(this->VK->Q), y0->value);
	mpz_set (this->VK->RwWQ[iIOW]->value, y0->value);
	iIOW+=1;
      }
      if (ptr->output_ctrb->length>0) {
	this->VK->RyYP[iIOY]->wire = ptr;
	ptr->Y->eval (ptr->Y, this->param[3], y0);  // y0 = Y(s)
	y0->field->mulElt (this->param[2], y0, y0); // y0 = Ry * Y(s)
	element_mul_mpz (this->VK->RyYP[iIOY]->elt, *(this->VK->P), y0->value);
	mpz_set (this->VK->RyYP[iIOY]->value, y0->value);
	iIOY+=1;
      }
    }//end if(input/output)
   
    else if ( strcmp(ptr->type,"intermediate")==0 ) {
      
      if (ptr->left_ctrb->length>0) {
	this->EK->RvVP[iIntV]->wire = ptr;
	this->EK->ALPHAvRvVP[iIntV]->wire = ptr;

	ptr->V->eval (ptr->V, this->param[3], y0);  //y0 = V(s)
	y0->field->mulElt (this->param[0], y0, y0); //y0 = Rv * V(s)
	y1->field->mulElt (this->param[7], y0, y1); //y1 = BETA * Rv * V(s)
	element_mul_mpz (this->EK->RvVP[iIntV]->elt, *(this->VK->P), y0->value);
	mpz_set (this->EK->RvVP[iIntV]->value, y0->value);
	y0->field->mulElt (this->param[4], y0, y0); //y0 = ALPHAv * Rv * V(s)
	element_mul_mpz (this->EK->ALPHAvRvVP[iIntV]->elt, *(this->VK->P), y0->value);
	mpz_set (this->EK->ALPHAvRvVP[iIntV]->value, y0->value);
	iIntV+=1;
      }
      
      if (ptr->right_ctrb->length>0) {
	this->EK->RwWQ[iIntW]->wire = ptr;
	this->EK->ALPHAwRwWQ[iIntW]->wire = ptr;
	this->EK->ALPHAwRwWP[iIntW]->wire = ptr;
	
	ptr->W->eval (ptr->W, this->param[3], y0);  //y0 = W(s)
	y0->field->mulElt (this->param[1], y0, y0); //y0 = Rw * W(s)
	y2->field->mulElt (this->param[7], y0, y2); //y2 = BETA * Rw * W(s)
	element_mul_mpz (this->EK->RwWQ[iIntW]->elt, *(this->VK->Q), y0->value);
	mpz_set (this->EK->RwWQ[iIntW]->value, y0->value);
	y0->field->mulElt (this->param[5], y0, y0); //y0 = ALPHAw * Rw * W(s)
	element_mul_mpz (this->EK->ALPHAwRwWQ[iIntW]->elt, *(this->VK->Q), y0->value);
	element_mul_mpz (this->EK->ALPHAwRwWP[iIntW]->elt, *(this->VK->P), y0->value);
	mpz_set (this->EK->ALPHAwRwWQ[iIntW]->value, y0->value);
	iIntW+=1;
      }
      
      if (ptr->output_ctrb->length>0) {
	this->EK->RyYP[iIntY]->wire = ptr;
	this->EK->ALPHAyRyYP[iIntY]->wire = ptr;

	ptr->Y->eval (ptr->Y, this->param[3], y0);  //y0 = Y(s)
	y0->field->mulElt (this->param[2], y0, y0); //y0 = Ry * Y(s)
	y3->field->mulElt (this->param[7], y0, y3); //y3 = BETA * Ry * Y(s)
	element_mul_mpz (this->EK->RyYP[iIntY]->elt, *(this->VK->P), y0->value);
	mpz_set (this->EK->RyYP[iIntY]->value, y0->value);
	y0->field->mulElt (this->param[6], y0, y0); //y0 = ALPHAy * Ry * Y(s)
	element_mul_mpz (this->EK->ALPHAyRyYP[iIntY]->elt, *(this->VK->P), y0->value);
	mpz_set (this->EK->ALPHAyRyYP[iIntY]->value, y0->value);
	iIntY+=1;
      }

      // y0 = (BETA * Rv * V(s)) + (BETA * Rw * W(s))
      y0->field->addElt (y1, y2, y0);

      // y1 = y0 + (BETA * Ry * Y(s))
      y0->field->addElt (y0, y3, y1);

      //gmp_printf(" y1 = %Zd\n", y1->value);
      element_mul_mpz (this->EK->RxBETAXP[iInt]->elt, *(this->VK->P),  y1->value);
      mpz_set (this->EK->RxBETAXP[iInt]->value, y1->value);
      this->EK->RxBETAXP[iInt]->wire = ptr;
      iInt+=1;

      //mpz_init( y1->value );
      //mpz_init( y2->value );
      //mpz_init( y3->value );
      mpz_set_ui (y1->value, 0);
      mpz_set_ui (y2->value, 0);
      mpz_set_ui (y3->value, 0);
    }//end else if(intermediate)
      
    else {
      printf ("error in generateRxXY, type not defined... exiting now");
      exit (0);
    }
    
  }

  y0->free (y0);
  y1->free (y1);
  y2->free (y2);
  y3->free (y3);
  
}

//================================
//================================

void generateTargetPoly (struct poly_t *this, struct circuit_t *C) {

  // temporary variables
  struct poly_t *tmp_poly1 = malloc (sizeof (struct poly_t) );
  tmp_poly1->init = initPoly;
  tmp_poly1->init (tmp_poly1);

  struct poly_t *tmp_poly2 = malloc (sizeof (struct poly_t) );
  tmp_poly2->init = initPoly;
  tmp_poly2->init (tmp_poly2);
  
  struct field_elt_t *tmp_elt = malloc (sizeof (struct field_elt_t) );
  tmp_elt->init = initFieldElt;
  tmp_elt->init (tmp_elt, C->field);

  // T(x) = 1*x^0
  mpz_set_ui (tmp_elt->value, 1); 
  this->insertNode (this, tmp_elt->value, 0, tmp_elt->field);
  
  int i;
  for (i=0; i<C->nMul_gate; i++) {
    mpz_set_ui (tmp_elt->value, 1);
    
    // tmp_poly1(x) = 1*x
    tmp_poly1->insertNode (tmp_poly1, tmp_elt->value, 1, tmp_elt->field);

    // tmp_elt <- (-ri)
    tmp_elt->field->negElt (C->gate[C->index[i]]->root, tmp_elt);

    // tmp_poly1(x) = 1*x - ri *x^0
    tmp_poly1->insertNode (tmp_poly1, tmp_elt->value, 0, tmp_elt->field);
        
    tmp_poly2->mulPoly (this, tmp_poly1, tmp_poly2);
    tmp_poly1->reset (tmp_poly1);
    
    this->reset (this);
    this->copyPoly (tmp_poly2, this);
    tmp_poly2->reset (tmp_poly2);
  }

  tmp_poly1->free (tmp_poly1);
  tmp_poly2->free (tmp_poly2);
  tmp_elt->free (tmp_elt);
}

//================================
//================================

void generateQAP (struct crs_t *this) {
  
  if (this->C==NULL) {
    printf("error in generateQAP, circuit does not exist... exiting now!");
    exit(0);
  }

  int i;
  for (i=0; i<this->C->nWire; i++) {

    if (this->C->wire[i]==NULL) {
      printf("error in generateQAP, wire[%d] does not exist... exiting now!", i);
      exit(0);
    }
    
    this->C->wire[i]->V = malloc (sizeof (struct poly_t) );
    this->C->wire[i]->W = malloc (sizeof (struct poly_t) );
    this->C->wire[i]->Y = malloc (sizeof (struct poly_t) );

    this->C->wire[i]->V->init = initPoly;
    this->C->wire[i]->W->init = initPoly;
    this->C->wire[i]->Y->init = initPoly;

    this->C->wire[i]->V->init (this->C->wire[i]->V);
    this->C->wire[i]->W->init (this->C->wire[i]->W);
    this->C->wire[i]->Y->init (this->C->wire[i]->Y);
    
    if (this->C->wire[i]->left_ctrb->length>0) {
      this->genQAPP (this->C->wire[i]->V, this, this->C->wire[i]->left_ctrb->head, this->C->wire[i], this->nIO_V, this->nInt_V);
    }

    if (this->C->wire[i]->right_ctrb->length>0) {
      this->genQAPP (this->C->wire[i]->W, this, this->C->wire[i]->right_ctrb->head, this->C->wire[i], this->nIO_W, this->nInt_W);
    }

     if (this->C->wire[i]->output_ctrb->length>0) {
      this->genQAPP (this->C->wire[i]->Y, this, this->C->wire[i]->output_ctrb->head, this->C->wire[i], this->nIO_Y, this->nInt_Y);
    }
    
  }// end-loop
}

//================================
//================================

void generateQAPPoly (struct poly_t *this, struct crs_t *CRS, struct node_t *node,
		      struct wire_t *wire, unsigned int *nIO, unsigned int *nInt) {
  
  // temporary variables
  struct poly_t *tmp_poly = malloc (sizeof (struct poly_t) );
  tmp_poly->init = initPoly;
  tmp_poly->init (tmp_poly);
 
  struct field_elt_t *tmp_elt = malloc (sizeof (struct field_elt_t) );
  tmp_elt->init = initFieldElt;
  tmp_elt->init (tmp_elt, CRS->C->field);
  mpz_set_ui (tmp_elt->value, 0);
  
  this->insertNode (this, tmp_elt->value, 0, tmp_elt->field);
  
  while (node!=NULL) {
    tmp_poly->addPoly (this, node->gate->lagrange_poly, tmp_poly);
    this->reset (this);
    this->copyPoly (tmp_poly, this);
    tmp_poly->reset (tmp_poly);
    node = node->next;
  }

  if (strcmp(wire->type,"input")==0 || strcmp(wire->type,"output")==0) {
    *(nIO)+=1;
  }
  else if (strcmp(wire->type,"intermediate")==0) {
    *(nInt)+=1;
  }
  else {
    printf("error in generateQAPP, type not defined... exiting now");
    exit(0);
  }
  
  tmp_elt->free (tmp_elt);
  tmp_poly->free (tmp_poly);
      
}

//================================
//================================

void generateLagrangePoly (struct crs_t *this) {
  int i, j;

  // temporary variables
  struct field_elt_t *tmp_elt = malloc (sizeof (struct field_elt_t) );
  tmp_elt->init = initFieldElt;
  tmp_elt->init (tmp_elt, this->C->field);

  struct field_elt_t *tmp_den = malloc (sizeof (struct field_elt_t) );
  tmp_den->init = initFieldElt;
  tmp_den->init (tmp_den, this->C->field);
    
  struct poly_t *tmp_ptr;
  struct poly_t *tmp_poly1 = malloc (sizeof (struct poly_t) );
  struct poly_t *tmp_poly2 = malloc (sizeof (struct poly_t) );
  
  tmp_poly1->init = initPoly;
  tmp_poly1->init (tmp_poly1);

  tmp_poly2->init = initPoly;
  tmp_poly2->init (tmp_poly2);
  
  for (i=0; i<this->C->nMul_gate; i++) {
    this->C->gate[this->C->index[i]]->lagrange_poly = malloc (sizeof (struct poly_t) );
    this->C->gate[this->C->index[i]]->lagrange_poly->init = initPoly;
    this->C->gate[this->C->index[i]]->lagrange_poly->init (this->C->gate[this->C->index[i]]->lagrange_poly);

    // temporary pointer
    tmp_ptr = this->C->gate[this->C->index[i]]->lagrange_poly;

    // tmp_ptr(x) = 1*x^0
    mpz_set_ui (tmp_elt->value, 1);
    tmp_ptr->insertNode (tmp_ptr, tmp_elt->value, 0, tmp_elt->field);
  
    // tmp_den = 1
    mpz_set_ui (tmp_den->value, 1);
    
    for (j=0; j<this->C->nMul_gate; j++) {  

      if (i!=j) {
	// tmp_poly1(x) = 1*x^1-rj*x^0
	mpz_set_ui (tmp_elt->value, 1);
	tmp_poly1->insertNode (tmp_poly1, tmp_elt->value, 1, tmp_elt->field); //tmp_poly1(x) = 1*x^1
		
	tmp_elt->field->negElt (this->C->gate[this->C->index[j]]->root, tmp_elt);
	tmp_poly1->insertNode (tmp_poly1, tmp_elt->value, 0, tmp_elt->field); //tmp_poly1(x) = 1*x^1 - rj*x^0
	
	// tmp_poly2(x) = tmp_poly1(x) * tmp_ptr(x)
	tmp_poly2->mulPoly (tmp_ptr, tmp_poly1, tmp_poly2);
	tmp_poly1->reset (tmp_poly1);
 
	// tmp_ptr(x) <- tmp_poly2(x)
	tmp_ptr->reset (tmp_ptr);
	tmp_ptr->copyPoly (tmp_poly2, tmp_ptr); 
	tmp_poly2->reset (tmp_poly2);
	  
	// tmp_elt = (ri-rj)
	tmp_elt->field->subElt
	  (this->C->gate[this->C->index[i]]->root, this->C->gate[this->C->index[j]]->root, tmp_elt);

	// tmp_den *= tmp_elt
	tmp_den->field->mulElt (tmp_elt, tmp_den, tmp_den); 	
	
      }    
    } // end-for
    
    // tmp_den <- 1/tmp_den
    tmp_elt->field->invElt (tmp_den, tmp_den);

    // tmp_poly1(x) = tmp_ptr(x) * tmp_den
    tmp_ptr->cstMulPoly (tmp_ptr, tmp_den, tmp_poly1);
    tmp_ptr->reset (tmp_ptr);
    tmp_ptr->copyPoly (tmp_poly1, tmp_ptr); 
    tmp_poly1->reset (tmp_poly1);

  } // end-for

  tmp_elt->free (tmp_elt);
  tmp_den->free (tmp_den);
  
  tmp_poly1->free (tmp_poly1);
  tmp_poly2->free (tmp_poly2);
}

//================================
//================================

void freeCRS (struct crs_t *this) {
  int i;

  if (this!=NULL) {

    // param
    if (this->param!=NULL) {
      for (i=0; i<NPARAMS; i++) {
	if (this->param[i]!=NULL) {
	  this->param[i]->free (this->param[i]);
	}
      }
      free (this->param);
    }
    
    // EKey
    if (this->EK!=NULL) {

      for (i=0; i<*(this->nInt_V); i++) {
	if (this->EK->RvVP[i]!=NULL) {
	  mpz_clear (this->EK->RvVP[i]->value);
	  element_clear (this->EK->RvVP[i]->elt);
	  free (this->EK->RvVP[i]);
	}
	if (this->EK->ALPHAvRvVP[i]!=NULL) {
	  mpz_clear (this->EK->ALPHAvRvVP[i]->value);
	  element_clear (this->EK->ALPHAvRvVP[i]->elt);
	  free (this->EK->ALPHAvRvVP[i]);
	}
      }
      free (this->EK->RvVP);
      free (this->EK->ALPHAvRvVP);
      
      for (i=0; i<*(this->nInt_W); i++) {
	if (this->EK->RwWQ[i]!=NULL) {
	  mpz_clear (this->EK->RwWQ[i]->value);
	  element_clear (this->EK->RwWQ[i]->elt);
	  free (this->EK->RwWQ[i]);
	}
	if (this->EK->ALPHAwRwWQ[i]!=NULL) {
	  mpz_clear (this->EK->ALPHAwRwWQ[i]->value);
	  element_clear (this->EK->ALPHAwRwWQ[i]->elt);
	  free (this->EK->ALPHAwRwWQ[i]);
	}
	if (this->EK->ALPHAwRwWP[i]!=NULL) {
	  mpz_clear (this->EK->ALPHAwRwWP[i]->value);
	  element_clear (this->EK->ALPHAwRwWP[i]->elt);
	  free (this->EK->ALPHAwRwWP[i]);
	}
      }
      free (this->EK->RwWQ);
      free (this->EK->ALPHAwRwWQ);
      free (this->EK->ALPHAwRwWP);
      
      for (i=0; i<*(this->nInt_Y); i++) {
	if (this->EK->RyYP[i]!=NULL) {
	  mpz_clear (this->EK->RyYP[i]->value);
	  element_clear (this->EK->RyYP[i]->elt);
	  free (this->EK->RyYP[i]);
	}
	if (this->EK->RyYP[i]!=NULL) {
	  mpz_clear (this->EK->ALPHAyRyYP[i]->value);
	  element_clear (this->EK->ALPHAyRyYP[i]->elt);
	  free (this->EK->ALPHAyRyYP[i]);
	}
      }
      free (this->EK->RyYP);
      free (this->EK->ALPHAyRyYP);
	
      for (i=0; i<(this->C->nWire - this->C->nIOWire); i++) {
	if (this->EK->RxBETAXP[i]!=NULL) {
	  mpz_clear (this->EK->RxBETAXP[i]->value);
	  element_clear (this->EK->RxBETAXP[i]->elt);
	  free (this->EK->RxBETAXP[i]);
	}
      }
      free (this->EK->RxBETAXP);
      
      for (i=0; i<=(int)(this->T->hiExp (this->T)); i++) {
	mpz_clear (this->EK->SQ[i]->value);
	element_clear (this->EK->SQ[i]->elt);
	free (this->EK->SQ[i]);
      }
      free (this->EK->SQ);
      
      free (this->EK);
    }

    // VKey
    if (this->VK!=NULL) {

      element_clear (this->VK->ALPHAvQ);
      element_clear (this->VK->ALPHAwQ);
      element_clear (this->VK->ALPHAwP);
      element_clear (this->VK->ALPHAyQ);
      element_clear (this->VK->BETAP);
      element_clear (this->VK->BETAQ);
      element_clear (this->VK->RyTP);
      
      for (i=0; i<*(this->nIO_V); i++) {
	if (this->VK->RvVP[i]!=NULL) {
	  mpz_clear (this->VK->RvVP[i]->value);
	  element_clear (this->VK->RvVP[i]->elt);
	  free (this->VK->RvVP[i]);
	}
      }
      free (this->VK->RvVP);
      
      for (i=0; i<*(this->nIO_W); i++) {
	if (this->VK->RwWQ[i]!=NULL) {
	  mpz_clear (this->VK->RwWQ[i]->value);
	  element_clear (this->VK->RwWQ[i]->elt);
	  free (this->VK->RwWQ[i]);
	}
      }
      free (this->VK->RwWQ);
      	
      for (i=0; i<*(this->nIO_Y); i++) {
	if (this->VK->RyYP[i]!=NULL) {
	  mpz_clear (this->VK->RyYP[i]->value);
	  element_clear (this->VK->RyYP[i]->elt);
	  free (this->VK->RyYP[i]);
	}
      }
      free (this->VK->RyYP);
       
      free (this->VK);
    }
       
    // T(x)
    if (this->T!=NULL) {
      this->T->free (this->T);
    }
    
    // nIO_V
    if (this->nIO_V!=NULL) {
      free (this->nIO_V);
    }

    // nIO_W
    if (this->nIO_W!=NULL) {
      free (this->nIO_W);
    }

    // nIO_Y
    if (this->nIO_Y!=NULL) {
      free (this->nIO_Y);
    }

    // nInt_V
    if (this->nInt_V!=NULL) {
      free (this->nInt_V);
    }

    // nInt_W
    if (this->nInt_W!=NULL) {
      free (this->nInt_W);
    }

    // nInt_Y
    if (this->nInt_Y!=NULL) {
      free (this->nInt_Y);
    }
    
    free (this);
  }
}

//================================
//================================

void viewCRS (struct crs_t *this) {
  int i;
  printf ("\n");
  printf (" **** CRS ****\n");
  
  for (i=0; i<NPARAMS; i++) {
    gmp_printf(" param[%d]: %Zd\n", i, this->param[i]->value);
  }

  printf("\n --- Lagrange polynomial ---\n");
  for (i=0; i<this->C->nMul_gate; i++) {
    gmp_printf(" Gate %Zd: L(x) =", this->C->gate[this->C->index[i]]->root->value);
    this->C->gate[this->C->index[i]]->lagrange_poly->view
      (this->C->gate[this->C->index[i]]->lagrange_poly);
  }

  printf("\n --- Wire ctrb ---\n");
  for (i=0; i<this->C->nWire; i++) {
    this->C->wire[i]->viewCtrb (this->C->wire[i]);
  }

  printf("\n --- Target polynomial ---\n");
  printf (" T(x) =");
  this->T->view (this->T);

  printf("\n --- VKey ---\n");
  
  element_printf(" P = %B\n", this->VK->P);
  element_printf(" Q = %B\n", this->VK->Q); 
  element_printf(" ALPHAvQ = %B\n", this->VK->ALPHAvQ);
  element_printf(" ALPHAwQ = %B\n", this->VK->ALPHAwQ);
  element_printf(" ALPHAwP = %B\n", this->VK->ALPHAwP);
  element_printf(" ALPHAyQ = %B\n", this->VK->ALPHAyQ);
  element_printf(" BETAP   = %B\n", this->VK->BETAP);
  element_printf(" BETAQ   = %B\n", this->VK->BETAQ);
  element_printf(" RyTP    = %B\n", this->VK->RyTP);

  for (i=0; i<*(this->nIO_V); i++) {
    //gmp_printf(" RvV (%d) = %Zd\n", this->VK->RvVP[i]->wire->id, this->VK->RvVP[i]->value);
    element_printf(" RvVP(%d) = %B\n", this->VK->RvVP[i]->wire->id, this->VK->RvVP[i]->elt);
  }
  
  for (i=0; i<*(this->nIO_W); i++) {
    //gmp_printf(" RwW (%d) = %Zd\n", this->VK->RwWQ[i]->wire->id, this->VK->RwWQ[i]->value);
    element_printf(" RwWQ(%d) = %B\n", this->VK->RwWQ[i]->wire->id, this->VK->RwWQ[i]->elt);
  }
  
  for (i=0; i<*(this->nIO_Y); i++) {
    //gmp_printf(" RyY (%d) = %Zd\n", this->VK->RyYP[i]->wire->id, this->VK->RyYP[i]->value);
    element_printf(" RyYP(%d) = %B\n", this->VK->RyYP[i]->wire->id, this->VK->RyYP[i]->elt);
  }

  printf("\n --- EKey ---\n");

  for (i=0; i<*(this->nInt_V); i++) {
    //gmp_printf(" RvV (%d) = %Zd\n", this->EK->RvVP[i]->wire->id, this->EK->RvVP[i]->value);
    element_printf(" RvVP(%d) = %B\n", this->EK->RvVP[i]->wire->id, this->EK->RvVP[i]->elt);
    //gmp_printf(" ALPHAvRvV (%d) = %Zd\n", this->EK->ALPHAvRvVP[i]->wire->id, this->EK->ALPHAvRvVP[i]->value);
    element_printf(" ALPHAvRvVP(%d) = %B\n", this->EK->ALPHAvRvVP[i]->wire->id, this->EK->ALPHAvRvVP[i]->elt);
  }
  
  for (i=0; i<*(this->nInt_W); i++) {
    //gmp_printf(" RwW (%d) = %Zd\n", this->EK->RwWQ[i]->wire->id, this->EK->RwWQ[i]->value);
    element_printf(" RwWQ(%d) = %B\n", this->EK->RwWQ[i]->wire->id, this->EK->RwWQ[i]->elt);
    //gmp_printf(" ALPHAwRwW (%d) = %Zd\n", this->EK->ALPHAwRwWQ[i]->wire->id, this->EK->ALPHAwRwWQ[i]->value);
    element_printf(" ALPHAwRwWQ(%d) = %B\n", this->EK->ALPHAwRwWQ[i]->wire->id, this->EK->ALPHAwRwWQ[i]->elt);
    element_printf(" ALPHAwRwWP(%d) = %B\n", this->EK->ALPHAwRwWP[i]->wire->id, this->EK->ALPHAwRwWP[i]->elt);
  }
  
  for (i=0; i<*(this->nInt_Y); i++) {
    //gmp_printf(" RyY (%d) = %Zd\n", this->EK->RyYP[i]->wire->id, this->EK->RyYP[i]->value);
    element_printf(" RyYP(%d) = %B\n", this->EK->RyYP[i]->wire->id, this->EK->RyYP[i]->elt);
    //gmp_printf(" ALPHAyRyY (%d) = %Zd\n", this->EK->ALPHAyRyYP[i]->wire->id, this->EK->ALPHAyRyYP[i]->value);
    element_printf(" ALPHAyRyYP(%d) = %B\n", this->EK->ALPHAyRyYP[i]->wire->id, this->EK->ALPHAyRyYP[i]->elt);
  }

  for (i=0; i<(this->C->nWire - this->C->nIOWire); i++) {
    //gmp_printf(" ( RvBETAV(%d) + RwBETAW(%d) + RyBETAY(%d) )     = %Zd\n",
    //this->EK->RxBETAXP[i]->wire->id, this->EK->RxBETAXP[i]->wire->id, this->EK->RxBETAXP[i]->wire->id, this->EK->RxBETAXP[i]->value);
    element_printf(" ( RvBETAV(%d) + RwBETAW(%d) + RyBETAY(%d) ) * P = %B\n",
		   this->EK->RxBETAXP[i]->wire->id, this->EK->RxBETAXP[i]->wire->id, this->EK->RxBETAXP[i]->wire->id, this->EK->RxBETAXP[i]->elt);
  }

  for (i=0; i<=(int)(this->T->hiExp (this->T)); i++) {
    //gmp_printf(" s%d  = %Zd\n", i, this->EK->SQ[i]->value);
    element_printf(" s%dQ = %B\n", i, this->EK->SQ[i]->elt);
  }
  
}

//================================
//================================
