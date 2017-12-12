#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "field.h"

//==============================================
//==============================================

void init_field( struct field_t *f, struct curve_t *c ) {
  mpz_init( f->p );
  mpz_set( f->p, c->r );
}

//==============================================
//==============================================

void free_field( struct field_t *f ) {
  mpz_clear( f->p );
  free( f );
}

//==============================================
//==============================================

void nonzero_rnd_field_elt( struct field_elt_t *e, gmp_randstate_t s ) {
  mpz_init (e->value);
  
  if (fDebug==false) {
    while ( mpz_cmp_ui(e->value, 0)==0 ) {
      mpz_urandomm( e->value, s, e->field->p );
    }
  }
  
  /* bound for mpz_urandomm */
  else {
    mpz_t upBnd;
    mpz_init_set_ui( upBnd, fOrder);
    while ( mpz_cmp_ui(e->value, 0)==0 ) {
      mpz_urandomm( e->value, s, upBnd );
    }
    mpz_clear( upBnd );
  }
}

//==============================================
//==============================================

void mul_field_elt( struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out ) {

  if ( in1->field != in2->field ) {
    gmp_printf(" %Zd %Zd",in1->field->p,in2->field->p);
    printf("error in mul_elt: in1 and in2 are not elements of the same field... exiting now\n!");
    exit(1);
  }

  out->field = in1->field;
  //mpz_init( out->value );

  mpz_mul( out->value, in1->value, in2->value );

  mod_field_elt( out );

  //if (fDebug==false) {
  //mpz_mod( out->value, out->value, out->field->p );
  //}
  //else {
  //mpz_t mod;
  //mpz_init_set_ui( mod, fOrder);
  //mpz_mod( out->value, out->value, mod );
  //mpz_clear( mod );
  //}
}

//==============================================
//==============================================

void div_field_elt( struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out ) {
  if ( in1->field != in2->field ) {
    printf("error in div_elt: in1 and in2 are not elements of the same field... exiting now!\n");
    exit(1);
  }
  
  struct field_elt_t *temp = malloc( sizeof( struct field_elt_t ));
  mpz_init( temp->value );

  inv_field_elt (temp, in2);//temp<-in2^-1
  
  //temp->field = in2->field;

  //out->field = in1->field;
  
  //if (fDebug==false) {
  //mpz_invert( temp->value, in2->value, in2->field->p );
    //mpz_mul( out->value, in1->value, temp->value );
    //mpz_mod( out->value, out->value, out->field->p );
  //}
  //else {
  //mpz_t mod;
  //mpz_init_set_ui( mod, fOrder);
  //mpz_invert( temp->value, in2->value, mod );
    //mpz_mul( out->value, in1->value, temp->value );
    //mpz_mod( out->value, out->value, mod );
    //mpz_clear( mod );
  //}
  mul_field_elt( in1, temp, out );
  //mpz_div( out->value, in1->value, in2->value );

  free_field_elt( temp );

  //free( temp->field );
  //mpz_clear( temp->value );
  //free( temp );
}

//==============================================
//==============================================

void add_field_elt( struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out ) {
  if ( in1->field != in2->field ) {
    printf("error in add_elt: in1 and in2 are not elements of the same field... exiting now!\n");
    exit(1);
  }
  out->field = in1->field;
  //mpz_init( out->value );
  mpz_add( out->value, in1->value, in2->value );

  mod_field_elt( out );
  
  //if (fDebug==false) {
  //mpz_mod( out->value, out->value, out->field->p );
  //}
  //else {
  //mpz_t mod;
  //mpz_init_set_ui( mod, fOrder);
  //mpz_mod( out->value, out->value, mod );
  //mpz_clear( mod );
  //}
}

//==============================================
//==============================================

void sub_field_elt( struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out ) {
  if ( in1->field != in2->field ) {
    printf("error in sub_elt: in1 and in2 are not elements of the same field... exiting now!\n");
    exit(1);
  }
  out->field = in1->field;
  //mpz_init( out->value );
  mpz_sub( out->value, in1->value, in2->value );

  mod_field_elt( out );
  
  //if (fDebug==false) {
  //mpz_mod( out->value, out->value, out->field->p );
  //}
  //else {
  //mpz_t mod;
  //mpz_init_set_ui( mod, fOrder);
  //mpz_mod( out->value, out->value, mod );
  //mpz_clear( mod );
  //}
}

//==============================================
//==============================================

void inv_field_elt( struct field_elt_t *out, struct field_elt_t *in ) {
  out->field = in->field;
  if (fDebug==false) {
    if( mpz_invert( out->value, in->value, out->field->p ) == 0 ) {
      printf("inverse does not exist... exiting now!");
      exit(0);
    }
  }
  else {
    mpz_t mod;
    mpz_init_set_ui( mod, fOrder);
    if( mpz_invert( out->value, in->value, mod ) == 0 ) {
      printf("inverse does not exist... exiting now!");
      exit(0);
    }
    mpz_clear( mod );
  }
}

//==============================================
//==============================================

void neg_field_elt( struct field_elt_t *out, struct field_elt_t *in ) {
  out->field = in->field;
  mpz_neg( out->value, in->value );
  mod_field_elt( out );
  //if (fDebug==false) {
  //mpz_mod( out->value, out->value, out->field->p );
  //}
  //else {
  //mpz_t mod;
  //mpz_init_set_ui( mod, fOrder);
  //mpz_mod( out->value, out->value, mod );
  //mpz_clear( mod );
  //}
}

//==============================================
//==============================================

void pow_ui_field_elt( struct field_elt_t *out, struct field_elt_t *in, unsigned int exp ) {
  out->field = in->field;
  mpz_pow_ui( out->value, in->value, exp);
  mod_field_elt( out );
  //if (fDebug==false) {
  //mpz_mod( out->value, out->value, out->field->p );
  //}
  //else {
  //mpz_t mod;
  //mpz_init_set_ui( mod, fOrder);
  //mpz_mod( out->value, out->value, mod );
  //mpz_clear( mod );
  //}
}

//==============================================
//==============================================

void insert_poly_ring_node( mpz_t coeff, unsigned int exp, struct field_t *f, struct poly_ring_t *p ) {

  struct poly_ring_node_t *node = malloc( sizeof ( struct poly_ring_node_t ) );
  node->exponent = exp;
  node->coefficient = malloc( sizeof ( struct field_elt_t ) );
  mpz_init( node->coefficient->value);
  mpz_set( node->coefficient->value, coeff );
  node->coefficient->field = f;
  node->next = NULL;
  if( p->length==0 ) {
    p->head = p->tail = node;
  }
  else {
    p->tail->next = node;
    p->tail = node;
  }
  p->length++; 
}

//==============================================
//==============================================

void print_poly( struct poly_ring_t *poly ) {
  if (poly->length!=0) {
    struct poly_ring_node_t *node = poly->head;

    while( node != NULL ) {  
      gmp_printf("%Zdx^%d", node->coefficient->value, node->exponent);
      
      node = node->next;
      
      if( node != NULL ) { printf(" + "); }
    }
    printf("\n");
    free( node );
  }
  else {
    printf("empty poly.\n");
  }
}

//==============================================
//==============================================

void mul_poly( struct poly_ring_t *in1, struct poly_ring_t *in2, struct poly_ring_t *out ) {

  if( in1 == NULL || in2 == NULL ) {
    printf("error empty poly(s)... exiting now\n");
    exit(0);
  }

  //printf(" in1:");
  //print_poly( in1 );

  //printf(" in2:");
  //print_poly( in2 );

  //printf(" out:");
  //print_poly( out );
  
  struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
  temp->length = 0;

  //temp->coefficient = malloc( sizeof( struct field_elt_t ));
  //mpz_init( temp->coefficient->value );
  
  struct poly_ring_node_t *node_in1 = in1->head;
  struct poly_ring_node_t *node_in2 = in2->head;

  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  
  unsigned int exp;
  
  while( node_in1 != NULL ) {
    while( node_in2 != NULL ) {
      //printf(" TINTIN\n");
      mul_field_elt( node_in1->coefficient, node_in2->coefficient, coeff);
      //printf(" MILOU\n");
      //gmp_printf(" printing coeffs\n %Zd\n %Zd\n %Zd\n", node_in1->coefficient, node_in2->coefficient, coeff);
      //gmp_printf(" base field %Zd\n", coeff->field->p );
      exp = node_in1->exponent + node_in2->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in2 = node_in2->next;
    }
    node_in1 = node_in1->next;
    node_in2 = in2->head;
  }
  //exit(0);
  
  //unsigned int largest_exp;
  //largest_exp = get_largest_exp( temp );
  
  //printf("largest exponent %d\n", largest_exp);

  out->length = temp->length;
  out->head = temp->head;
  out->tail = temp->tail;

  //printf("---- LEAVING MUL_POLY ----\n");
  //print_poly(temp);
  //print_poly(out);
  //exit(0);
  //free_poly_ring( temp->head );
  //free( node_in1 );
  //free( node_in2 );
  //free_field_elt( coeff );
}

//==============================================
//==============================================

unsigned int get_largest_exp( struct poly_ring_t *poly ) {

  if( poly == NULL ) {
    printf("error empty element(s)... exiting now\n");
    exit(0);
  }

  unsigned int exp = 0;
  struct poly_ring_node_t *node_poly = poly->head;

  while( node_poly != NULL ) {
    if( node_poly->exponent > exp ) {
      exp = node_poly->exponent;
    }
    node_poly = node_poly->next;
  }

  /* freeing memory */
  free( node_poly );

  return exp;
}

//==============================================
//==============================================

void poly_ordering( struct poly_ring_t *p ) {

  struct poly_ring_node_t *node = p->head;

  struct poly_ring_node_t *next_node = node->next;

  struct poly_ring_node_t *temp_node;

  int count = 0;

  //print_poly( p );
  
  while( node != NULL ) {
    //printf(" node exp %d\n",node->exponent);
    //print_poly( p );
    while( next_node != NULL ) {
      if( node->exponent == next_node->exponent ) {
	add_field_elt( node->coefficient, next_node->coefficient, node->coefficient);
	temp_node = next_node->next;
	delete_poly_ring_node( p, next_node );
	next_node = temp_node;
      }
      else {
	next_node = next_node->next;
      }
    }
    node = node->next;
    if (node!= NULL) {
      next_node = node->next;
    }
  }
  //print_poly( p );
}

//==============================================
//==============================================

void delete_poly_ring_node( struct poly_ring_t *p, struct poly_ring_node_t *node_delete ) {
  struct poly_ring_node_t *node = p->head;

  if (node==NULL) {
    printf("leaving delete_poly_ring_node");
  }
  
  if (node!=NULL && node==node_delete) {
    p->head=node->next;
    p->length-=1;
    free(node_delete);
    return;
  }

  while (node!=NULL && node->next!=node_delete)
    node = node->next;
  
  if (node_delete!=p->tail) {
    node->next=node_delete->next;
    p->length-=1;
    free(node_delete);
    return;
  }
  else {
    node->next=NULL;
    p->tail=node;
    p->length-=1;
    free(node_delete);
    return;
  }

}

//==============================================
//==============================================

void evaluate_poly( struct poly_ring_t *poly, struct field_elt_t *x, struct field_elt_t *y ) {
  
  if (poly->length!=0) {
    struct poly_ring_node_t *node = poly->head;
    
    struct field_elt_t *temp = malloc( sizeof ( struct field_elt_t ) );
    mpz_init( temp->value );
    
    while( node != NULL ) {
      pow_ui_field_elt( temp, x, node->exponent );    //temp <- x^exp      [mod p]
      mul_field_elt( node->coefficient, temp, temp ); //temp <- coeff*temp [mod p]
      add_field_elt( temp, y, y );                    //y <- y+temp        [mod p]
      
      //gmp_printf("%Zdx^%d", node->coefficient->value, node->exponent);
      node = node->next;
    }

    //free (node);
    //free_field_elt (temp);   
  }
  
  else {
    printf("empty polynomial... existing now!");
    exit(0);
  }
}

//==============================================
//==============================================

void add_poly( struct poly_ring_t *in1, struct poly_ring_t *in2, struct poly_ring_t *out ) {
  if( in1 == NULL || in2 == NULL ) {
    printf("error empty poly(s)... exiting now\n");
    exit(0);
  }
  
  struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
  temp->length = 0;
  
  struct poly_ring_node_t *node_in1 = in1->head;
  struct poly_ring_node_t *node_in2 = in2->head;
  
  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  
  unsigned int exp;
  
  while( node_in1 != NULL && node_in2 != NULL  ) {
    if( node_in1->exponent == node_in2->exponent ) {
      add_field_elt( node_in1->coefficient, node_in2->coefficient, coeff);
      exp = node_in1->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in1 = node_in1->next;
      node_in2 = node_in2->next;
      if (node_in2==NULL) {
	while( node_in1 != NULL ) {
	  mpz_set( coeff->value, node_in1->coefficient->value );
	  coeff->field = node_in1->coefficient->field;
	  exp = node_in1->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in1 = node_in1->next;
	}
      }
      if (node_in1==NULL) {
	while( node_in2 != NULL ) {
	  mpz_set( coeff->value, node_in2->coefficient->value );
	  coeff->field = node_in2->coefficient->field;
	  exp = node_in2->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in2 = node_in2->next;
	}
      }
    }
    
    else if( node_in1->exponent > node_in2->exponent ) {
      mpz_set( coeff->value, node_in1->coefficient->value );
      coeff->field = node_in1->coefficient->field;
      exp = node_in1->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in1 = node_in1->next;
      if (node_in1==NULL) {
	while( node_in2 != NULL ) {
	  mpz_set( coeff->value, node_in2->coefficient->value );
	  coeff->field = node_in2->coefficient->field;
	  exp = node_in2->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in2 = node_in2->next;
	}
      }
    }
    
    else if( node_in1->exponent < node_in2->exponent ) {
      mpz_set( coeff->value, node_in2->coefficient->value );
      coeff->field = node_in2->coefficient->field;
      exp = node_in2->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in2 = node_in2->next;
      if (node_in2==NULL) {
	while( node_in1 != NULL ) {
	  mpz_set( coeff->value, node_in1->coefficient->value );
	  coeff->field = node_in1->coefficient->field;
	  exp = node_in1->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in1 = node_in1->next;
	}
      }
    }

    else { printf("error in add_poly... exiting now!\n"); exit(0); }
    
  }
  
  out->length = temp->length;
  out->head = temp->head;
  out->tail = temp->tail;
  
  //free_poly_ring( temp->head );
  //free( node_in1 );
  //free( node_in2 );
  //free_field_elt( coeff );
}

//==============================================
//==============================================

void sub_poly( struct poly_ring_t *in1, struct poly_ring_t *in2, struct poly_ring_t *out ) {
  if( in1 == NULL || in2 == NULL ) {
    printf("error empty poly(s)... exiting now\n");
    exit(0);
  }
  
  struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
  temp->length = 0;
  
  struct poly_ring_node_t *node_in1 = in1->head;
  struct poly_ring_node_t *node_in2 = in2->head;
  
  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  
  unsigned int exp;
  
  while( node_in1 != NULL && node_in2 != NULL  ) {
    if( node_in1->exponent == node_in2->exponent ) {
      sub_field_elt( node_in1->coefficient, node_in2->coefficient, coeff);
      exp = node_in1->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in1 = node_in1->next;
      node_in2 = node_in2->next;
      if (node_in1==NULL) {
	while( node_in2 != NULL ) {
	  struct field_elt_t *foo = malloc( sizeof( struct field_elt_t ));
	  mpz_init( foo->value );
	  neg_field_elt( foo, node_in2->coefficient );
	  mpz_set( coeff->value, foo->value );
	  coeff->field = node_in2->coefficient->field;
	  exp = node_in2->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in2 = node_in2->next;
	}
      }
      if (node_in2==NULL) {
	while( node_in1 != NULL ) {
	  mpz_set( coeff->value, node_in1->coefficient->value );
	  coeff->field = node_in1->coefficient->field;
	  exp = node_in1->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in1 = node_in1->next;
	}
      }
    }
    
    else if( node_in1->exponent > node_in2->exponent ) {
      mpz_set( coeff->value, node_in1->coefficient->value );
      coeff->field = node_in1->coefficient->field;
      exp = node_in1->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in1 = node_in1->next;
      if (node_in1==NULL) {
	while( node_in2 != NULL ) {
	  struct field_elt_t *foo = malloc( sizeof( struct field_elt_t ));
	  mpz_init( foo->value );
	  neg_field_elt( foo, node_in2->coefficient );
	  mpz_set( coeff->value, foo->value );
	  coeff->field = node_in2->coefficient->field;
	  exp = node_in2->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in2 = node_in2->next;
	}
      }
    }
    
    else if( node_in1->exponent < node_in2->exponent ) {
      mpz_set( coeff->value, node_in2->coefficient->value );
      coeff->field = node_in2->coefficient->field;
      exp = node_in2->exponent;
      insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
      node_in2 = node_in2->next;
      if (node_in2==NULL) {
	while( node_in1 != NULL ) {
	  mpz_set( coeff->value, node_in1->coefficient->value );
	  coeff->field = node_in1->coefficient->field;
	  exp = node_in1->exponent;
	  insert_poly_ring_node( coeff->value, exp, coeff->field, temp );
	  node_in1 = node_in1->next;
	}
      }
    }
    
    else { printf("error in sub_poly... exiting now!\n"); exit(0); }
    
  }
  
  out->length = temp->length;
  out->head = temp->head;
  out->tail = temp->tail;

  //free_poly_ring( temp->head );
  //free( node_in1 );
  //free( node_in2 );
  //free_field_elt( coeff );
  
}

//==============================================
//==============================================

void div_poly( struct poly_ring_t *in1, struct poly_ring_t *in2, struct poly_ring_t *out ) {

  printf("--- enter div_poly -----\n");
  
  if( in1 == NULL || in2 == NULL ) {
    printf("error empty poly(s)... exiting now\n");
    exit(0);
  }

  struct poly_ring_t *dividend = malloc( sizeof( struct poly_ring_t ));
  dividend->length = 0;
  
  struct poly_ring_t *quotient = malloc( sizeof( struct poly_ring_t ));
  quotient->length = 0;  

  struct poly_ring_t *remainder = malloc( sizeof( struct poly_ring_t ));
  remainder->length = 0;

  struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
  temp->length = 0;
  
  dividend->length = in1->length;
  dividend->head = in1->head;
  dividend->tail = in1->tail;

  struct poly_ring_node_t *node_dividend = dividend->head;
  struct poly_ring_node_t *node_divisor  = in2->head;
  
  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  
  unsigned int hi_exp;
  hi_exp = get_largest_exp(dividend) - get_largest_exp(in2);

  unsigned int exp;
  
  if (hi_exp<0) {
    coeff->field = in1->head->coefficient->field;
    insert_poly_ring_node( coeff->value, 0, coeff->field, out );
  }
  else {
    while (node_dividend->exponent >= get_largest_exp(in2) && dividend->length > 0) {
      //printf("\n");
      //printf(" dividend:");
      //print_poly( dividend );
      //printf(" divisor:");
      //print_poly( in2 );
      div_field_elt( node_dividend->coefficient, node_divisor->coefficient, coeff);
      exp = node_dividend->exponent - node_divisor->exponent;
      //gmp_printf(" %d %Zd %Zd %Zd\n",exp, coeff->value,node_dividend->coefficient->value,node_divisor->coefficient->value);
      insert_poly_ring_node( coeff->value, exp, coeff->field, quotient );
      //printf(" quotient:");
      //print_poly( quotient );
      insert_poly_ring_node( coeff->value, exp, coeff->field, out );
      //printf(" out:");
      //print_poly( out );
      mul_poly( quotient, in2, temp );
      //printf(" temp:");
      //print_poly( temp );
      //sub_poly( in1, temp, out );
      sub_poly( dividend, temp, dividend );
      //printf(" dividend:");
      //print_poly( dividend );
      delete_poly_ring_node( dividend, dividend->head );
      delete_poly_ring_node( quotient, quotient->head );
      //printf(" dividend:");
      //print_poly( dividend );
      //printf(" quotient:");
      //print_poly( quotient );
      //printf (" dividend length %d:\n", dividend->length);
      //gmp_printf ("---> %Zd\n", node_dividend->coefficient->value);
      if (dividend->length>0) {node_dividend = dividend->head;}//node_dividend = node_dividend->next;}
      //gmp_printf (" ---> %Zd\n", node_dividend->coefficient->value);
    }
    insert_poly_ring_node( coeff->value, exp, coeff->field, remainder );
    remainder->length = dividend->length;
    remainder->head = dividend->head;
    remainder->tail = remainder->tail;
    //printf(" out:");
    //print_poly( out );
    printf(" remainder:");
    print_poly( remainder );
    //exit(0);
  }

  //out->length = temp->length;
  //out->head = temp->head;
  //out->tail = temp->tail;

  free(temp);
}

//==============================================
//==============================================

void cst_mul_poly( struct poly_ring_t *in, struct field_elt_t *c, struct poly_ring_t *out ) {

  struct poly_ring_t *temp = malloc( sizeof( struct poly_ring_t ));
  temp->length = 0;
  
  struct poly_ring_node_t *node_in = in->head;

  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );

  while( node_in != NULL ) {
    mul_field_elt( node_in->coefficient, c, coeff);
    insert_poly_ring_node( coeff->value, node_in->exponent, coeff->field, temp );
    node_in = node_in->next;
  }

  out->length = temp->length;
  out->head = temp->head;
  out->tail = temp->tail;
  
  free(temp);
}

//==============================================
//==============================================

void random_monic_poly( struct field_t *f, struct poly_ring_t *poly, unsigned int hi_exp, gmp_randstate_t s ) {

  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  coeff->field = f;
  
  int i;
  for (i=hi_exp; i>=0; i--) {
    
    if (i==hi_exp) {
      mpz_set_ui (coeff->value, 1); 
      insert_poly_ring_node( coeff->value, hi_exp, coeff->field, poly );
    }
    else {
      nonzero_rnd_field_elt (coeff, s);
      insert_poly_ring_node( coeff->value, i, coeff->field, poly );
    }
  }
  free (coeff);
}

//==============================================
//==============================================

void random_poly( struct field_t *f, struct poly_ring_t *poly, unsigned int hi_exp, gmp_randstate_t s ) {

  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  coeff->field = f;
  
  int i;
  for (i=hi_exp; i>=0; i--) {
    nonzero_rnd_field_elt (coeff, s);
    insert_poly_ring_node( coeff->value, i, coeff->field, poly );
  }
  free (coeff);
}

//==============================================
//==============================================

bool is_poly_irreducible( struct field_t *f, struct poly_ring_t *poly ) {
  
  struct poly_ring_t *u = malloc( sizeof( struct poly_ring_t ));
  u->length=0;

  struct field_elt_t *coeff = malloc( sizeof( struct field_elt_t ));
  mpz_init( coeff->value );
  coeff->field = f;

  mpz_set_ui (coeff->value, 1); 
  insert_poly_ring_node( coeff->value, 1, coeff->field, u ); //u(x)<-x

  print_poly( u );

  unsigned int m;
  m = floor( 0.5 * (double) (get_largest_exp(poly)) );

  printf("%d %d\n",get_largest_exp(poly),m);
  
  //int i;
  //for (i=1; i<=m; i++) {
  //}
  
  free (u);
  return true;
}

//==============================================
//==============================================

void free_poly_ring( struct poly_ring_node_t *head ) {
  
  struct poly_ring_node_t *temp;
  
  while (head != NULL) {
    temp = head;
    head = head->next;
    free (temp);
  }

}

//==============================================
//==============================================

void free_field_elt( struct field_elt_t *elt ) {

  elt->field=NULL;
  free( elt->field );

  mpz_clear( elt->value );

  free( elt );
}

//==============================================
//==============================================

void mod_field_elt( struct field_elt_t *elt ) {
  mpz_t foo;
  mpz_init( foo );
  mpz_set_str( foo, "205523667896953300194895899082072403858390252929", 10 );

  //if (mpz_cmp(foo, elt->field->p)!=0 ) {
  // gmp_printf(" OOPS %Zd\n", elt->field->p);
  //exit(0);
  //}
  
  if (fDebug==false) {
    mpz_mod( elt->value, elt->value, elt->field->p );
  }
  
  else {
    mpz_t mod;
    mpz_init_set_ui( mod, fOrder);
    mpz_mod( elt->value, elt->value, mod );
    mpz_clear( mod );
  }

  mpz_clear( foo );
}
