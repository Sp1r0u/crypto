#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "wire.h"

void wire_contribution_to_mul_gate( struct wire_contribution_t *ctrb, struct gate_t *g ) {
  struct node_t *node = malloc( sizeof ( struct node_t ) );
  
  node->gate = g;
  node->next = NULL;
  
  if( ctrb->length==0 ) {
    ctrb->head = ctrb->tail = node;
  }
  else {
    ctrb->tail->next = node;
    ctrb->tail = node;
  }

  ctrb->length++;

}

//==============================================
//==============================================

void init_wire_contribution( struct wire_t *wire ) {
 
  //wire->left_ctrb = (struct wire_contribution_t*) calloc( 1, sizeof( struct wire_contribution_t ));
  wire->left_ctrb   = malloc( sizeof( struct wire_contribution_t ));
  wire->right_ctrb  = malloc( sizeof( struct wire_contribution_t ));
  wire->output_ctrb = malloc( sizeof( struct wire_contribution_t ));
 
  wire->left_ctrb->length   = 0;
  wire->right_ctrb->length  = 0;
  wire->output_ctrb->length = 0;
  /*
  wire->left_ctrb->head = NULL;
  wire->left_ctrb->tail = NULL;

  wire->right_ctrb->head = NULL;
  wire->right_ctrb->tail = NULL;

  wire->output_ctrb->head = NULL;
  wire->output_ctrb->tail = NULL;

  wire->V = malloc( sizeof( struct poly_ring_t ));
  wire->V->length = 0;
  wire->V->head = NULL;
  wire->V->tail = NULL;
  */

  wire->V = NULL;
  wire->W = NULL;
  wire->Y = NULL;
  
}

//==============================================
//==============================================

void print_wire_contribution_to_mul_gate( struct wire_t *w) {

  if (w->left_ctrb->length>0) {
    struct node_t *node = w->left_ctrb->head;
    //printf("\nwire ID %d left ctrb(s) to mult. gates:\n",w->id);
    while( node != NULL ) {
      gmp_printf(" %Zd",node->gate->root->value);
      node = node->next;
    } 
    free( node );
  }
  else {
    //printf("\nwire ID %d NO left ctrb to mult. gates\n",w->id);
    printf(" no left ctrb.");
  }
  
  if (w->right_ctrb->length>0) {
    struct node_t *node = w->right_ctrb->head;
    //printf("\nwire ID %d right ctrb(s) to mult. gates\n",w->id);
    while( node != NULL ) {
      gmp_printf(" %Zd",node->gate->root->value);
      node = node->next;
    } 
    free( node );
  }
  else {
    //printf("wire ID %d NO right ctrb to mult. gates\n",w->id);
    printf(" no right ctrb.");
  }

   if (w->output_ctrb->length>0) {
    struct node_t *node = w->output_ctrb->head;
    //printf("\nwire ID %d right ctrb(s) to mult. gates\n",w->id);
    while( node != NULL ) {
      gmp_printf(" %Zd",node->gate->root->value);
      node = node->next;
    } 
    free( node );
  }
  else {
    //printf("wire ID %d NO right ctrb to mult. gates\n",w->id);
    printf(" no output ctrb.");
  }

}

//==============================================
//==============================================

void free_wire( struct wire_t *wire ) {

  if (wire!=NULL) {

    /* type */
    if (wire->type!=NULL) {
      free (wire->type);
    }

    /* value */
    if (wire->value!=NULL) {
      free_field_elt (wire->value);
    }
 
    /* left ctrb */
    if (wire->left_ctrb->length==0) {
      free (wire->left_ctrb);
    }
    else {
      free_wire_contribution (wire->left_ctrb->head);
      free (wire->left_ctrb);
    }

    /* right ctrb */
    if (wire->right_ctrb->length==0) {
      free (wire->right_ctrb);
    }
    else {
      free_wire_contribution (wire->right_ctrb->head);
      free (wire->right_ctrb);
    }

    /* output ctrb */
    if (wire->output_ctrb->length==0) {
      free (wire->output_ctrb);
    }
    else {
      free_wire_contribution (wire->output_ctrb->head);
      free (wire->output_ctrb);
    }
    
    /* V */
    if (wire->V!=NULL) {
      if (wire->V->length==0) {
	free (wire->V);
      }
      else {
	free_poly_ring (wire->V->head);
	free (wire->V);
      }
    }
    
    /* W */
    if (wire->W!=NULL) {
      if (wire->W->length==0) {
	free (wire->W);
      }
      else {
	free_poly_ring (wire->W->head);
      }
    }
    
    /* Y */
    if (wire->Y!=NULL) {
      if (wire->Y->length==0) {
	free (wire->Y);
      }
      else {
	free_poly_ring (wire->Y->head);
	free (wire->Y);
      }
    }

    free (wire);
  }
}

//==============================================
//==============================================

void free_wire_contribution( struct node_t *head ) {
  struct node_t *temp;

  while (head != NULL) {
    temp = head;
    head = head->next;
    free (temp);
  }

}
