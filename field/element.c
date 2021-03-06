#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <gmp.h>

#include "element.h"

//==================================================
//==================================================

void multiply_elements( struct element_t *e1, struct element_t *e2, struct element_t *e3 ) {
  /**
   * Polynomial multiplication function
   * Multiplies two elements to a given variable
   * e1: the first polynomial expression
   * e2: the second polynomial expression
   * e3: the result
   **/
  if( e1 == NULL || e2 == NULL ) {
    printf("error empty element(s)... exiting now\n");
    exit(0);
  }

  if( e1->field != e2->field ) {
    printf("error elements are not defined over the same ground field... exiting now\n");
    exit(0);
  }

  struct element_t *temp = malloc( sizeof( struct element_t ));
  temp->length = 0;
  temp->field = e1->field;
  
  struct element_node_t *node_e1 = e1->head;
  struct element_node_t *node_e2 = e2->head;
  
  mpz_t coeff;
  mpz_init( coeff );

  mpz_t mod;
  mpz_init( mod );
  mpz_set( mod, e1->field->n );

  int exp;

  while( node_e1 != NULL ) {
    while( node_e2 != NULL ) {
      mpz_mul( coeff, node_e1->coeff, node_e2->coeff);
      mpz_mod( coeff, coeff, mod );

      exp = node_e1->exp + node_e2->exp;

      //gmp_printf("%Zd %Zd %Zd %Zd %d", node_e1->coeff, node_e2->coeff, coeff, mod, exp);
      //printf("\n");

      set_element( temp, coeff, exp );
      node_e2 = node_e2->next;
    }
    node_e1 = node_e1->next;
    node_e2 = e2->head;
  }

  print_element( temp ); 
  printf("length( temp ) %d\n", get_element_length(temp));

  int hi_exp;
  hi_exp = get_element_highest_degree( temp );
  printf("highest_degree %d\n", hi_exp);

  e3->length = 0;
  e3->field  = e1->field;

  int i;
  struct element_node_t *node_temp = temp->head;
  for( i=0; i<=hi_exp; i++ ) {
    mpz_init( coeff );
    while( node_temp != NULL ) {
      if( node_temp->exp == i ) {
	mpz_add( coeff, coeff, node_temp->coeff);
      }
      node_temp = node_temp->next;
    }
    mpz_mod( coeff, coeff, mod );
    gmp_printf("exp %d coeff %Zd\n", i, coeff);
    set_element( e3, coeff, i );
    node_temp = temp->head;
  }
  
  /* freeing memory */
  mpz_clear( coeff );
  mpz_clear( mod );
  free( node_e1 );
  free( node_e2 );
  free( node_temp );
  free( temp );
}

//==================================================
//==================================================

int get_element_length( struct element_t *e ) {
  return e->length;
}

//==================================================
//==================================================

void set_element( struct element_t *e, mpz_t coeff, unsigned int exp ) {

  struct element_node_t *node = malloc( sizeof ( struct element_node_t ) );
  
  mpz_init( node->coeff );

  mpz_set( node->coeff, coeff );

  node->exp = exp;
  
  node->next = NULL;

  //gmp_printf("%Zdx^%d", node->coeff, node->exp);
  //printf("\n");
  
  if( e->length==0 ) {
    e->head = e->tail = node;
  }
  else {
    e->tail->next = node;
    e->tail = node;
  }

  e->length++;

  /********************************************************************************
   * source: http://www.geeksforgeeks.org/adding-two-polynomials-using-linked-list/
   *******************************************************************************/
  /*  
  struct element_t *r;
  struct element_t *z;

  z = *temp;

  if( z==NULL ) {
    r = (struct element_t*)malloc(sizeof(struct element_t));
    mpz_init( r->coeff );
    mpz_set( r->coeff, coeff );
    r->exp = exp;
    *temp = r;
    r->next = (struct element_t*)malloc(sizeof(struct element_t));
    r = r->next;
    r->next = NULL;
  }

  else {
    r = (struct element_t*)malloc(sizeof(struct element_t));
    mpz_init( r->coeff );
    mpz_set( r->coeff, coeff );
    r->exp = exp;
    r->next = (struct element_t*)malloc(sizeof(struct element_t));
    r = r->next;
    r->next = NULL;
  }
  */
}

//==================================================
//==================================================

void add_elements( struct element_t *e1, struct element_t *e2, struct element_t *e3 ) {

  if( e1->field != e2->field ) {
    printf("error elements are not defined over the same ground field... exiting now\n");
    exit(0);
  }

  struct element_t *temp = malloc( sizeof( struct element_t ));
  temp->length = 0;
  temp->field = e1->field;

  struct element_node_t *node_e1 = e1->head;
  struct element_node_t *node_e2 = e2->head;

  mpz_t coeff;
  mpz_init( coeff );

  mpz_t mod;
  mpz_init( mod );
  mpz_set( mod, e1->field->n );

  int exp;
  
  //while( node_e1 != NULL ) {
  //while( node_e2 != NULL ) {
  while( node_e1 != NULL && node_e2 != NULL  ) {     
    if( node_e1->exp == node_e2->exp ) {
      mpz_add( coeff, node_e1->coeff, node_e2->coeff);
      mpz_mod( coeff, coeff, mod );
      
      exp = node_e1->exp;
      
      set_element( temp, coeff, exp );
      node_e1 = node_e1->next;
      node_e2 = node_e2->next;
    }
    else if( node_e1->exp > node_e2->exp ) {
      set_element( temp, node_e1->coeff, node_e1->exp );
      node_e1 = node_e1->next;
    }
    else if( node_e1->exp < node_e2->exp ) {
      set_element( temp, node_e2->coeff, node_e2->exp );
      node_e2 = node_e2->next;
    }
    else { printf("error in add_elements... exiting now!\n"); exit(0); }
  }
  //}

  while( node_e1 != NULL || node_e2 != NULL  ) {
    if( node_e1 != NULL ) {
      
    }
  }
  
  print_element( temp ); 
  printf("length( temp ) %d\n", get_element_length(temp));
  
  /* freeing memory */
  mpz_clear( coeff );
  mpz_clear( mod );
  free( node_e1 );
  free( node_e2 );
  free( temp );
  
}

//==================================================
//==================================================

void print_element( struct element_t *e ) {

  struct element_node_t *node = e->head;

  //int i;
  //for( i=0; i<e->length; i++ ) {

  while( node != NULL ) {  
    gmp_printf("%Zdx^%d", node->coeff, node->exp);

    node = node->next;

    if( node != NULL ) { printf(" + "); }
  }

  printf("\n");
  
  gmp_printf("field order %Zd\n", e->field->n); 

  free( node );
}

//==================================================
//==================================================

int get_element_highest_degree( struct element_t *e ) {

  if( e == NULL ) {
    printf("error empty element(s)... exiting now\n");
    exit(0);
  }

  int exp=-1;
  struct element_node_t *node_e = e->head;

  while( node_e != NULL ) {
    if( node_e->exp > exp ) {exp=node_e->exp;}
    node_e = node_e->next;
  }

  /* freeing memory */
  free( node_e );

  return exp;
}
