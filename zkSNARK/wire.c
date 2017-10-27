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

void init_wire_contribution( struct wire_t *w ) {
  w->left_ctrb  = malloc( sizeof( struct wire_contribution_t ));
  w->right_ctrb = malloc( sizeof( struct wire_contribution_t ));

  w->left_ctrb->length = 0;
  w->right_ctrb->length = 0;
}


