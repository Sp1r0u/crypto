#include <string.h>

#include "CRS.h"

//==============================================
//==============================================

void generate_EKVK( struct crs_t *crs, struct circuit_t *c, struct curve_t *ec, struct field_t *f ) {
  int i;

  struct poly_ring_t *target = malloc( sizeof( struct poly_ring_t )); //univariate polynomial ring
  generate_target_poly( c, f, target );

  struct field_elt_t *target_eval = malloc( sizeof ( struct field_elt_t ) );
  mpz_init( target_eval->value );

  target_eval->field = crs->param[3]->field;
  
  evaluate_poly( target, crs->param[3], target_eval);

  printf(" target poly.: ");
  print_poly( target );
  //gmp_printf(" t(%Zd) = %Zd\n", crs->param[3]->value, target_eval->value);

  //printf(" nIO_V %d nInt_V %d\n",*(crs->nIO_V), *(crs->nInt_V));
  //printf(" nIO_W %d nInt_W %d\n",*(crs->nIO_W), *(crs->nInt_W));
  //printf(" nIO_Y %d nInt_Y %d\n",*(crs->nIO_Y), *(crs->nInt_Y));
  
  crs->VK->rvVP = malloc( *( crs->nIO_V ) * sizeof( struct extrakey_t* ));
  crs->VK->rwWQ = malloc( *( crs->nIO_W ) * sizeof( struct extrakey_t* ));
  crs->VK->ryYP = malloc( *( crs->nIO_Y ) * sizeof( struct extrakey_t* ));

  crs->EK->rvVP = malloc( *( crs->nInt_V ) * sizeof( struct extrakey_t* ));
  crs->EK->rwWQ = malloc( *( crs->nInt_W ) * sizeof( struct extrakey_t* ));
  crs->EK->ryYP = malloc( *( crs->nInt_Y ) * sizeof( struct extrakey_t* ));
  
  //crs->EK->rvVP = malloc( ( c->nWire - c->nIOWire ) * sizeof( struct extrakey_t* ));
  //crs->EK->rwWQ = malloc( ( c->nWire - c->nIOWire ) * sizeof( struct extrakey_t* ));
  //crs->EK->ryYP = malloc( ( c->nWire - c->nIOWire ) * sizeof( struct extrakey_t* ));
  
  crs->EK->alpha_vrvVP = malloc( *( crs->nInt_V ) * sizeof( struct extrakey_t* ));
  crs->EK->alpha_wrwWQ = malloc( *( crs->nInt_W ) * sizeof( struct extrakey_t* ));
  crs->EK->alpha_yryYP = malloc( *( crs->nInt_Y ) * sizeof( struct extrakey_t* ));

  crs->EK->r_xbetaVP = malloc( ( c->nWire - c->nIOWire ) * sizeof( struct extrakey_t* ));
  
  crs->VK->P = &(ec->P);
  crs->VK->Q = &(ec->Q);
  
  element_init_G2( crs->VK->alpha_vQ, ec->pairing );
  element_init_G2( crs->VK->alpha_wQ, ec->pairing );
  element_init_G1( crs->VK->alpha_wP, ec->pairing );
  element_init_G2( crs->VK->alpha_yQ, ec->pairing );
  element_init_G1( crs->VK->betaP,    ec->pairing );
  element_init_G2( crs->VK->betaQ,    ec->pairing );
  element_init_G1( crs->VK->rytP,     ec->pairing );

  for (i=0; i<*(crs->nIO_V); i++) {
    crs->VK->rvVP[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G1( crs->VK->rvVP[i]->elt, ec->pairing );
  }

  for (i=0; i<*(crs->nIO_W); i++) {
    crs->VK->rwWQ[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G2( crs->VK->rwWQ[i]->elt, ec->pairing );
  }

  for (i=0; i<*(crs->nIO_Y); i++) {
    crs->VK->ryYP[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G1( crs->VK->ryYP[i]->elt, ec->pairing );
  }

  for (i=0; i<*(crs->nInt_V); i++) {
    crs->EK->rvVP[i] = malloc( sizeof( struct extrakey_t ));
    crs->EK->alpha_vrvVP[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G1( crs->EK->rvVP[i]->elt, ec->pairing );
    element_init_G1( crs->EK->alpha_vrvVP[i]->elt, ec->pairing );
  }

  for (i=0; i<*(crs->nInt_W); i++) {
    crs->EK->rwWQ[i] = malloc( sizeof( struct extrakey_t ));
    crs->EK->alpha_wrwWQ[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G2( crs->EK->rwWQ[i]->elt, ec->pairing );
    element_init_G2( crs->EK->alpha_wrwWQ[i]->elt, ec->pairing );
  }

  for (i=0; i<*(crs->nInt_Y); i++) {
    crs->EK->ryYP[i] = malloc( sizeof( struct extrakey_t ));
    crs->EK->alpha_yryYP[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G1( crs->EK->ryYP[i]->elt, ec->pairing );
    element_init_G1( crs->EK->alpha_yryYP[i]->elt, ec->pairing );
  }

  for (i=0; i<(c->nWire-c->nIOWire); i++) {
    crs->EK->r_xbetaVP[i] = malloc( sizeof( struct extrakey_t ));
    element_init_G1( crs->EK->r_xbetaVP[i]->elt, ec->pairing );
  }
  
  element_mul_mpz( crs->VK->alpha_vQ, *(crs->VK->Q), crs->param[4]->value );
  element_mul_mpz( crs->VK->alpha_wQ, *(crs->VK->Q), crs->param[5]->value );
  element_mul_mpz( crs->VK->alpha_wP, *(crs->VK->P), crs->param[5]->value );
  element_mul_mpz( crs->VK->alpha_yQ, *(crs->VK->Q), crs->param[6]->value );
  element_mul_mpz( crs->VK->betaP,    *(crs->VK->P), crs->param[7]->value );
  element_mul_mpz( crs->VK->betaQ,    *(crs->VK->Q), crs->param[7]->value );
  element_mul_mpz( crs->VK->rytP,     *(crs->VK->P), crs->param[2]->value ); //rytp=ry*P
  element_mul_mpz( crs->VK->rytP,     crs->VK->rytP, target_eval->value   ); //rytp=t(s)*rytp

  evaluate_rxXY( c, crs );
  
  print_EKVK( c, crs );

  free_poly_ring (target->head);
  free (target);
  free_field_elt (target_eval);

  target=NULL;
  target_eval=NULL;
  
}

//==============================================
//==============================================

void generate_crs( struct crs_t *crs, struct circuit_t *c, struct curve_t *ec, struct field_t *f, gmp_randstate_t rnd_state ) {
  int i;
  printf ("\n**** Generating CRS ****\n\n");

  crs->param = malloc( NPARAMS * sizeof( struct field_elt_t* ));
  
  crs->EK = malloc( sizeof( struct EKey_t ));
   
  crs->VK = malloc( sizeof( struct VKey_t ));

  for (i=0; i<NPARAMS; i++) {
    crs->param[i] = malloc( sizeof( struct field_elt_t ));
        
    crs->param[i]->field = f;
    mpz_init( crs->param[i]->value );

    if (i!=2) {
      nonzero_rnd_field_elt( crs->param[i], rnd_state);
    }
    else {
      mul_field_elt( crs->param[0], crs->param[1], crs->param[2] ); //ry=rv*rw [mod p]
    }
  }

  for (i=0; i<NPARAMS; i++) {
    gmp_printf(" param[%d]: %Zd\n", i, crs->param[i]->value);
  }

  generate_Lagrange_polynomial( c, f );

  printf("\n --- Lagrange polynomial ---\n");
  for (i=0; i<c->nMul_gate; i++) {
    gmp_printf(" Gate %Zd: ", c->gate[c->index[i]]->root->value);
    print_poly( c->gate[c->index[i]]->lagrange_poly );
  }
  

  crs->nIO_V = malloc( sizeof( unsigned int ));
  crs->nIO_W = malloc( sizeof( unsigned int ));
  crs->nIO_Y = malloc( sizeof( unsigned int ));

  crs->nInt_V = malloc( sizeof( unsigned int ));
  crs->nInt_W = malloc( sizeof( unsigned int ));
  crs->nInt_Y = malloc( sizeof( unsigned int ));
  
  *(crs->nIO_V) = 0;
  *(crs->nIO_W) = 0;
  *(crs->nIO_Y) = 0;

  *(crs->nInt_V) = 0;
  *(crs->nInt_W) = 0;
  *(crs->nInt_Y) = 0;

  generate_QAP_polynomial( c, f, crs );

  printf("\n --- QAP polynomial ---\n");
  for (i=0; i<c->nWire; i++) {
    gmp_printf(" wire %d:", c->wire[i]->id);
    printf(" v(x)=");
    print_poly( c->wire[i]->V );
    printf("         w(x)=");
    print_poly( c->wire[i]->W );
    printf("         y(x)=");
    print_poly( c->wire[i]->Y );
  }

  generate_EKVK( crs, c, ec, f );

}

//==============================================
//==============================================

void generate_Lagrange_polynomial( struct circuit_t *c, struct field_t *f ) {
  int i, j;

  mpz_t foo;
  mpz_init(foo);

  struct field_elt_t *inv = malloc( sizeof ( struct field_elt_t ) ); //inverse of denominator
  mpz_init( inv->value ); 
  //inv->field = f;
  
  struct field_elt_t *temp1 = malloc( sizeof ( struct field_elt_t ) );
  mpz_init( temp1->value );
  //temp1->field = f;
 
  struct field_elt_t *temp2 = malloc( sizeof ( struct field_elt_t ) );
  mpz_init_set_ui( temp2->value, 1 );
  temp2->field = f;
  
  for (i=0; i<c->nMul_gate; i++) {

    gmp_printf("mul_gate: %Zd\n", c->gate[c->index[i]]->root->value);
    
    c->gate[c->index[i]]->lagrange_poly = malloc( sizeof( struct poly_ring_t ));

    c->gate[c->index[i]]->lagrange_poly->length = 0;

    mpz_set_ui(foo, 1);
    insert_poly_ring_node( foo, 0, f, c->gate[c->index[i]]->lagrange_poly );

    //print_poly( c->gate[c->index[i]]->lagrange_poly );

    struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
    temp->length=0;
    
    for (j=0; j<c->nMul_gate; j++) {
      
      if (i!=j) {
	temp->length=0;
	
	mpz_set_ui(foo, 1);
	insert_poly_ring_node( foo, 1, f, temp );

	mpz_neg( foo, c->gate[c->index[j]]->root->value );

	if (fDebug==false) {
	  mpz_mod( foo, foo, c->gate[c->index[j]]->root->field->p );
	}
	
	else {
	  mpz_t mod;
	  mpz_init_set_ui( mod, fOrder);
	  mpz_mod( foo, foo, mod );
	  mpz_clear( mod );
	}
	insert_poly_ring_node( foo, 0, f, temp );

	mul_poly( c->gate[c->index[i]]->lagrange_poly, temp, c->gate[c->index[i]]->lagrange_poly );

	sub_field_elt( c->gate[c->index[i]]->root, c->gate[c->index[j]]->root, temp1 ); //temp1 = x_i - x_j

	mul_field_elt( temp1, temp2, temp2 ); //temp2 = temp2 * (x_i-x_j)
      }
      
    }//end j-loop

    poly_ordering( c->gate[c->index[i]]->lagrange_poly );

    //print_poly( c->gate[c->index[i]]->lagrange_poly );
    //gmp_printf("temp2 = %Zd\n", temp2);

    inv_field_elt( inv, temp2 );
    //gmp_printf("inv = %Zd\n", inv);

    cst_mul_poly( c->gate[c->index[i]]->lagrange_poly, inv, c->gate[c->index[i]]->lagrange_poly);
    //print_poly( c->gate[c->index[i]]->lagrange_poly );
  
    mpz_set_ui( inv->value, 1 );
    mpz_set_ui( temp2->value, 1 );

    if (temp->length>0) {
      free_poly_ring (temp->head);
    }
    free (temp);
    temp=NULL;
  }//end i-loop

  free_field_elt (inv);
  free_field_elt (temp1);
  free_field_elt (temp2);
  mpz_clear( foo );

  //exit(0);
}

//==============================================
//==============================================

void generate_target_poly( struct circuit_t *c, struct field_t *f, struct poly_ring_t *t ) {

  mpz_t foo;
  mpz_init(foo);

  t->length = 0;
  
  mpz_set_ui(foo, 1);
  insert_poly_ring_node( foo, 0, f, t );

  int i;
  for (i=0; i<c->nMul_gate; i++) {

    struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
    temp->length=0;
    
    mpz_set_ui(foo, 1);
    insert_poly_ring_node( foo, 1, f, temp );

    mpz_neg( foo, c->gate[c->index[i]]->root->value );
    insert_poly_ring_node( foo, 0, f, temp );
    
    mul_poly( t, temp, t );

    free_poly_ring( temp->head);
    free (temp);
    temp = NULL;
    
  }
  
  poly_ordering( t );
}

//==============================================
//==============================================

void generate_Tpoly( struct circuit_t *c, struct field_t *f, struct poly_ring_t *t ) {

  mpz_t foo;
  mpz_init (foo);

  t->length = 0;
  
  mpz_set_ui (foo, 1);
  insert_poly_ring_node (foo, 0, f, t);

  int i;
  for (i=0; i<c->nMul_gate; i++) {

    struct poly_ring_t *tmp = malloc (sizeof (struct poly_ring_t));
    tmp->length=0;
    
    mpz_set_ui (foo, 1);
    insert_poly_ring_node (foo, 1, f, tmp);

    mpz_neg (foo, c->gate[c->index[i]]->root->value);
    insert_poly_ring_node (foo, 0, f, tmp);
    
    mul_poly (t, tmp, t);

    free_poly_ring (tmp->head);
    free (tmp);
    
  }
  
  poly_ordering (t);

}

//==============================================
//==============================================

void print_EKVK( struct circuit_t *c, struct crs_t *crs ) {
  
  printf("\n --- VKey ---\n");

  int i;
  
  element_printf(" P = %B\n", crs->VK->P);
  element_printf(" Q = %B\n", crs->VK->Q); 
  element_printf(" alpha_vQ = %B\n", crs->VK->alpha_vQ);
  element_printf(" alpha_wQ = %B\n", crs->VK->alpha_wQ);
  element_printf(" alpha_wP = %B\n", crs->VK->alpha_wP);
  element_printf(" alpha_yQ = %B\n", crs->VK->alpha_yQ);
  element_printf(" betaP    = %B\n", crs->VK->betaP);
  element_printf(" betaQ    = %B\n", crs->VK->betaQ);
  element_printf(" rytP     = %B\n", crs->VK->rytP);

  for (i=0; i<*(crs->nIO_V); i++) {
    element_printf(" rvVP(%d) = %B\n", crs->VK->rvVP[i]->wire->id, crs->VK->rvVP[i]->elt);
  }
  
  for (i=0; i<*(crs->nIO_W); i++) {
    element_printf(" rwWQ(%d) = %B\n", crs->VK->rwWQ[i]->wire->id, crs->VK->rwWQ[i]->elt);
  }
  
  for (i=0; i<*(crs->nIO_Y); i++) {
    element_printf(" ryYP(%d) = %B\n", crs->VK->ryYP[i]->wire->id, crs->VK->ryYP[i]->elt);
  }

  printf("\n --- EKey ---\n");

  for (i=0; i<*(crs->nInt_V); i++) {
    element_printf(" rvVP(%d) = %B\n", crs->EK->rvVP[i]->wire->id, crs->EK->rvVP[i]->elt);
    element_printf(" alpha_vrvVP(%d) = %B\n", crs->EK->alpha_vrvVP[i]->wire->id, crs->EK->alpha_vrvVP[i]->elt);
  }
  
  for (i=0; i<*(crs->nInt_W); i++) {
    element_printf(" rwWQ(%d) = %B\n", crs->EK->rwWQ[i]->wire->id, crs->EK->rwWQ[i]->elt);
    element_printf(" alpha_wrwWQ(%d) = %B\n", crs->EK->alpha_wrwWQ[i]->wire->id, crs->EK->alpha_wrwWQ[i]->elt);
  }
  
  for (i=0; i<*(crs->nInt_Y); i++) {
    element_printf(" ryYP(%d) = %B\n", crs->EK->ryYP[i]->wire->id, crs->EK->ryYP[i]->elt);
    element_printf(" alpha_yryYP(%d) = %B\n", crs->EK->alpha_yryYP[i]->wire->id, crs->EK->alpha_yryYP[i]->elt);
  }

  for (i=0; i<(c->nWire-c->nIOWire); i++) {
    element_printf(" ( rvbetaV(%d) + rwbetaW(%d) + rybetaY(%d) ) * P = %B\n", crs->EK->r_xbetaVP[i]->wire->id,
		   crs->EK->r_xbetaVP[i]->wire->id, crs->EK->r_xbetaVP[i]->wire->id, crs->EK->r_xbetaVP[i]->elt);
  }
  
}

//==============================================
//==============================================

void generate_QAP_polynomial( struct circuit_t *c, struct field_t *f, struct crs_t *crs ) {

  mpz_t foo;
  mpz_init(foo); 
  
  int i;
  
  if (c==NULL) {
    printf("in generate_QAP_polynomial, c does not exist... exiting now!");
    exit(0);
  }
  
  for (i=0; i<c->nWire; i++) {

    if (c->wire[i]==NULL) {
      printf("in generate_QAP_polynomial, wire does not exist... exiting now!");
      exit(0);
    }
    
    c->wire[i]->V = malloc( sizeof( struct poly_ring_t ));
    c->wire[i]->W = malloc( sizeof( struct poly_ring_t ));
    c->wire[i]->Y = malloc( sizeof( struct poly_ring_t ));
    
    c->wire[i]->V->length = 0;
    c->wire[i]->W->length = 0;
    c->wire[i]->Y->length = 0;
    
    mpz_set_ui(foo, 0);
    
    // V(x)
    if (c->wire[i]->left_ctrb==NULL) {
      printf("in generate_QAP_polynomial, left_ctrb does not exist... exiting now!");
      exit(0);
    }

    //printf(" length %d\n", c->wire[i]->left_ctrb->length);
    
    if (c->wire[i]->left_ctrb->length>0) {
      insert_poly_ring_node( foo, 0, f, c->wire[i]->V );
      print_poly(c->wire[i]->V);
      
      struct node_t *node = c->wire[i]->left_ctrb->head;
      
      while( node != NULL ) {
	add_poly( c->wire[i]->V, node->gate->lagrange_poly, c->wire[i]->V );
	node = node->next;
      }
      if ( strcmp(c->wire[i]->type,"input")==0 || strcmp(c->wire[i]->type,"output")==0 ) {*(crs->nIO_V)+=1;}
      else if ( strcmp(c->wire[i]->type,"intermediate")==0 ) {*(crs->nInt_V)+=1;}
      else { printf("wire's type not defined... exiting now"); exit(0);}
    }

    // W(x)
    if (c->wire[i]->right_ctrb->length>0) {
      insert_poly_ring_node( foo, 0, f, c->wire[i]->W );
      
      struct node_t *node = c->wire[i]->right_ctrb->head;
      
      while( node != NULL ) {
	add_poly( c->wire[i]->W, node->gate->lagrange_poly, c->wire[i]->W );
	node = node->next;
      }
      if ( strcmp(c->wire[i]->type,"input")==0 || strcmp(c->wire[i]->type,"output")==0 ) {*(crs->nIO_W)+=1;}
      else if ( strcmp(c->wire[i]->type,"intermediate")==0 ) {*(crs->nInt_W)+=1;}
      else { printf("wire's type not defined... exiting now"); exit(0);}
    }
    
    // Y(x)
    if (c->wire[i]->output_ctrb->length>0) {
      insert_poly_ring_node( foo, 0, f, c->wire[i]->Y );
      
      struct node_t *node = c->wire[i]->output_ctrb->head;
      
      while( node != NULL ) {
	add_poly( c->wire[i]->Y, node->gate->lagrange_poly, c->wire[i]->Y );
	node = node->next;
      }
      if ( strcmp(c->wire[i]->type,"input")==0 || strcmp(c->wire[i]->type,"output")==0 ) {*(crs->nIO_Y)+=1;}
      else if ( strcmp(c->wire[i]->type,"intermediate")==0 ) {*(crs->nInt_Y)+=1;}
      else { printf("wire's type not defined... exiting now"); exit(0);}
    }
    
  }// end-loop
  mpz_clear (foo);  
}

//==============================================
//==============================================

void evaluate_rxXY( struct circuit_t *c, struct crs_t *crs) {

  unsigned int iIOV = 0;
  unsigned int iIOW = 0;
  unsigned int iIOY = 0;
  
  unsigned int iIntV = 0;
  unsigned int iIntW = 0;
  unsigned int iIntY = 0;

  unsigned int iInt = 0;

  struct field_elt_t *y = malloc( sizeof ( struct field_elt_t ) );
  y->field = crs->param[3]->field;

  struct field_elt_t *y1 = malloc( sizeof ( struct field_elt_t ) );
  y1->field = crs->param[3]->field;

  struct field_elt_t *y2 = malloc( sizeof ( struct field_elt_t ) );
  y2->field = crs->param[3]->field;

  struct field_elt_t *y3 = malloc( sizeof ( struct field_elt_t ) );
  y3->field = crs->param[3]->field;

  mpz_init( y1->value );
  mpz_init( y2->value );
  mpz_init( y3->value );
    
  int i;
  for (i=0; i<c->nWire; i++) {

    if ( strcmp(c->wire[i]->type,"input")==0 || strcmp(c->wire[i]->type,"output")==0 ) {
      if (c->wire[i]->left_ctrb->length>0) {
	crs->VK->rvVP[iIOV]->wire = c->wire[i];
	mpz_init( y->value );
	evaluate_poly( c->wire[i]->V, crs->param[3], y );
	mul_field_elt( crs->param[0], y, y );
	element_mul_mpz( crs->VK->rvVP[iIOV]->elt, *(crs->VK->P),  y->value);
	//printf(" wire %d:", crs->VK->rvVP[iIOV]->wire->id);
	//gmp_printf(" rv*v(s) = %Zd\n", y->value);
	//element_printf("         rvVP = %B\n", crs->VK->rvVP[iIOV]->elt);
	//printf (" %d ptr add %p\n", iIOV, crs->VK->rvVP[iIOV]);
	iIOV+=1;
      }
 
      if (c->wire[i]->right_ctrb->length>0) {
	crs->VK->rwWQ[iIOW]->wire = c->wire[i];
	mpz_init( y->value );
	evaluate_poly( c->wire[i]->W, crs->param[3], y );
	mul_field_elt( crs->param[1], y, y );
	element_mul_mpz( crs->VK->rwWQ[iIOW]->elt, *(crs->VK->Q),  y->value);
	//printf(" wire %d:", crs->VK->rwWQ[iIOW]->wire->id);
	//gmp_printf(" rw*w(s) = %Zd\n", y->value);
	//element_printf("         rwWQ = %B\n", crs->VK->rwWQ[iIOW]->elt);
	iIOW+=1;
      }
      if (c->wire[i]->output_ctrb->length>0) {
	crs->VK->ryYP[iIOY]->wire = c->wire[i];
	mpz_init( y->value );
	evaluate_poly( c->wire[i]->Y, crs->param[3], y );
	mul_field_elt( crs->param[2], y, y );
	element_mul_mpz( crs->VK->ryYP[iIOY]->elt, *(crs->VK->P),  y->value);
	//printf(" wire %d:", crs->VK->ryYP[iIOY]->wire->id);
	//gmp_printf(" ry*y(s) = %Zd\n", y->value);
	//element_printf("         ryYP = %B\n", crs->VK->ryYP[iIOY]->elt);
	iIOY+=1;
      }
    }//end if(input/output)
    
    else if ( strcmp(c->wire[i]->type,"intermediate")==0 ) {
      if (c->wire[i]->left_ctrb->length>0) {
	crs->EK->rvVP[iIntV]->wire = c->wire[i];
	crs->EK->alpha_vrvVP[iIntV]->wire = c->wire[i];
	mpz_init( y->value );
	evaluate_poly( c->wire[i]->V, crs->param[3], y );//y  = v(s)
	mul_field_elt( crs->param[0], y, y );            //y  = rv * v(s)
	mul_field_elt( crs->param[7], y, y1 );           //y1 = beta * rv * v(s)
	element_mul_mpz( crs->EK->rvVP[iIntV]->elt, *(crs->VK->P),  y->value);//y = rv * v(s) * P
	mul_field_elt( crs->param[4], y, y );            //y  = alpha_v * rv * v(s)
	element_mul_mpz( crs->EK->alpha_vrvVP[iIntV]->elt, *(crs->VK->P),  y->value);//y = alpha_v * rv * v(s) * P
	//printf(" wire %d:", crs->EK->rvVP[iIntV]->wire->id);
	//gmp_printf(" rv*v(s) = %Zd\n", y->value);
	//element_printf("         rvVP = %B\n", crs->EK->rvVP[iIntV]->elt);
	iIntV+=1;
      }
      if (c->wire[i]->right_ctrb->length>0) {
	crs->EK->rwWQ[iIntW]->wire = c->wire[i];
	crs->EK->alpha_wrwWQ[iIntW]->wire = c->wire[i];
	mpz_init( y->value );
	evaluate_poly( c->wire[i]->W, crs->param[3], y );//y  = w(s)
	mul_field_elt( crs->param[1], y, y );            //y  = rw * w(s)
	mul_field_elt( crs->param[7], y, y2 );           //y2 = beta * rw * w(s)
	element_mul_mpz( crs->EK->rwWQ[iIntW]->elt, *(crs->VK->Q),  y->value);//y = rw * w(s) * Q
	mul_field_elt( crs->param[5], y, y );            //y  = alpha_w * rw * w(s)
	element_mul_mpz( crs->EK->alpha_wrwWQ[iIntW]->elt, *(crs->VK->Q),  y->value);//y = alpha_w * rw * w(s) * Q
	//printf(" wire %d:", crs->EK->rwWQ[iIntW]->wire->id);
	//gmp_printf(" rw*w(s) = %Zd\n", y->value);
	//element_printf("         rwWQ = %B\n", crs->EK->rwWQ[iIntW]->elt);
	iIntW+=1;
      }
      if (c->wire[i]->output_ctrb->length>0) {
	crs->EK->ryYP[iIntY]->wire = c->wire[i];
	crs->EK->alpha_yryYP[iIntY]->wire = c->wire[i];
	mpz_init( y->value );
	evaluate_poly( c->wire[i]->Y, crs->param[3], y );//y  = y(s)
	mul_field_elt( crs->param[2], y, y );            //y  = ry * y(s)
	mul_field_elt( crs->param[7], y, y3 );           //y3 = beta * ry * y(s)
	element_mul_mpz( crs->EK->ryYP[iIntY]->elt, *(crs->VK->P),  y->value);//y = ry * y(s) * P
	mul_field_elt( crs->param[6], y, y );            //y  = alpha_y * ry * y(s)
	element_mul_mpz( crs->EK->alpha_yryYP[iIntY]->elt, *(crs->VK->P),  y->value);//y = alpha_y * ry * y(s) * P
	//printf(" wire %d:", crs->EK->ryYP[iIntY]->wire->id);
	//gmp_printf(" ry*y(s) = %Zd\n", y->value);
	//element_printf("         ryYP = %B\n", crs->EK->ryYP[iIntY]->elt);
	iIntY+=1;
      }
      //gmp_printf(" y1 = %Zd y2 = %Zd y3 = %Zd\n", y1->value, y2->value, y3->value);
      add_field_elt( y1, y2, y1 ); //y1 = ( beta * rv * v(s) ) + ( beta * rw * w(s) )
      add_field_elt( y1, y3, y1 ); //y1 = ( beta * rv * v(s) ) + ( beta * rw * w(s) ) + ( beta * ry * y(s) )
      //gmp_printf(" y1 = %Zd\n", y1->value);
      element_mul_mpz( crs->EK->r_xbetaVP[iInt]->elt, *(crs->VK->P),  y1->value);//y1 = ( ( beta * rv * v(s) ) + ( beta * rw * w(s) ) + ( beta * ry * y(s) ) ) * P
      crs->EK->r_xbetaVP[iInt]->wire = c->wire[i];
      iInt+=1;
      mpz_init( y1->value );
      mpz_init( y2->value );
      mpz_init( y3->value );
    }//end else if(intermediate)
      
    else {
      printf("wire's type not defined... exiting now");
      exit(0);
    }
    
  }
  free_field_elt (y);
  free_field_elt (y1);
  free_field_elt (y2);
  free_field_elt (y3);

}

//==============================================
//==============================================

void free_crs( struct crs_t *crs, struct circuit_t *c  ) {

  int i;

  printf (" freeing param\n");
  if (crs->param!=NULL) {
    for (i=0; i<NPARAMS; i++) {
      if (crs->param[i]!=NULL) {
	free_field_elt (crs->param[i]);
      }
    }
    free (crs->param);
  }
  
  /* EK */
  printf (" freeing EK\n");
  if (crs->EK!=NULL) {

    for (i=0; i<*(crs->nInt_V); i++) {
      if (crs->EK->rvVP[i]!=NULL) {
	element_clear( crs->EK->rvVP[i]->elt );
	free (crs->EK->rvVP[i]->wire);
	free (crs->EK->rvVP[i]);
      }
    }
    free (crs->EK->rvVP);

    for (i=0; i<*(crs->nInt_W); i++) {
      if (crs->EK->rwWQ[i]!=NULL) {
	element_clear( crs->EK->rwWQ[i]->elt );
	free (crs->EK->rwWQ[i]->wire);
	free (crs->EK->rwWQ[i]);
      }
    }
    free (crs->EK->rwWQ);

    for (i=0; i<*(crs->nInt_Y); i++) {
      if (crs->EK->ryYP[i]!=NULL) {
	element_clear( crs->EK->ryYP[i]->elt );
	free (crs->EK->ryYP[i]->wire);
	free (crs->EK->ryYP[i]);
      }
    }
    free (crs->EK->ryYP);

    for (i=0; i<*(crs->nInt_V); i++) {
      if (crs->EK->alpha_vrvVP[i]!=NULL) {
	element_clear( crs->EK->alpha_vrvVP[i]->elt );
	free (crs->EK->alpha_vrvVP[i]->wire);
	free (crs->EK->alpha_vrvVP[i]);
	crs->EK->alpha_vrvVP[i]->wire=NULL;
	crs->EK->alpha_vrvVP[i]=NULL;
      }
    }
    free (crs->EK->alpha_vrvVP);

     for (i=0; i<*(crs->nInt_W); i++) {
      if (crs->EK->alpha_wrwWQ[i]!=NULL) {
	element_clear( crs->EK->alpha_wrwWQ[i]->elt );
	free (crs->EK->alpha_wrwWQ[i]->wire);
	free (crs->EK->alpha_wrwWQ[i]);
      }
    }
    free (crs->EK->alpha_wrwWQ);

    for (i=0; i<*(crs->nInt_Y); i++) {
      if (crs->EK->alpha_yryYP[i]!=NULL) {
	element_clear( crs->EK->alpha_yryYP[i]->elt );
	free (crs->EK->alpha_yryYP[i]->wire);
	free (crs->EK->alpha_yryYP[i]);
      }
    }
    free (crs->EK->ryYP);

    for (i=0; i<(c->nWire-c->nIOWire); i++) {
      if (crs->EK->r_xbetaVP[i]!=NULL) {
	element_clear( crs->EK->r_xbetaVP[i]->elt );
	//free (crs->EK->r_xbetaVP[i]->wire);
	free (crs->EK->r_xbetaVP[i]);
      }
    }
    free (crs->EK->r_xbetaVP);
        
    free (crs->EK);
  }

  /* VK */
  printf (" freeing VK\n");
  if (crs->VK!=NULL) {

    element_clear( crs->VK->alpha_vQ );
    element_clear( crs->VK->alpha_wQ );
    element_clear( crs->VK->alpha_wP );
    element_clear( crs->VK->alpha_yQ );
    element_clear( crs->VK->betaP );
    element_clear( crs->VK->betaQ );
    element_clear( crs->VK->rytP  );
    
    for (i=0; i<*(crs->nIO_V); i++) {
      if (crs->VK->rvVP[i]!=NULL) {
	element_clear( crs->VK->rvVP[i]->elt );
	free (crs->VK->rvVP[i]->wire);
	free (crs->VK->rvVP[i]);
      }
    }
    free (crs->VK->rvVP);
    
    for (i=0; i<*(crs->nIO_W); i++) {
      if (crs->VK->rwWQ[i]!=NULL) {
	element_clear( crs->VK->rwWQ[i]->elt );
	free (crs->VK->rwWQ[i]->wire);
	free (crs->VK->rwWQ[i]);
      }
    }
    free (crs->VK->rwWQ);

    for (i=0; i<*(crs->nIO_Y); i++) {
      if (crs->VK->ryYP[i]!=NULL) {
	element_clear( crs->VK->ryYP[i]->elt );
	free (crs->VK->ryYP[i]->wire);
	free (crs->VK->ryYP[i]);
      }
    }
    free (crs->VK->ryYP);
    
    free (crs->VK);
  }

  printf (" freeing nIO_X\n");
  if (crs->nIO_V!=NULL) {
    free (crs->nIO_V);
  }

  if (crs->nIO_W!=NULL) {
    free (crs->nIO_W);
  }

  if (crs->nIO_Y!=NULL) {
    free (crs->nIO_Y);
  }

  printf (" freeing nIint_X\n");
  if (crs->nInt_V!=NULL) {
    free (crs->nInt_V);
  }
  
  if (crs->nInt_W!=NULL) {
    free (crs->nInt_W);
  }
  
  if (crs->nInt_Y!=NULL) {
    free (crs->nInt_Y);
  }
   
  free (crs);
  
}
