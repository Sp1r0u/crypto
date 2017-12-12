#include <stdlib.h>
#include <stdio.h>

#include "wire.h"

void initCtrb (struct wire_contribution_t *this) {
  this->free = freeCtrb;
  this->length = 0;
}

//================================
//================================

void freeCtrb (struct wire_contribution_t *this) {
  struct node_t *tmp;
  
  while (this->head!=NULL) {
    tmp = this->head;
    this->head = this->head->next;
    free (tmp);
  }

  free (this);
}

//================================
//================================

void initWire (struct wire_t *this) {
 
  this->view = viewWire;
  this->free = freeWire;
  this->ctrb = ctrb2gate;

  this->viewCtrb = viewWireCtrb;
  
  this->type = malloc (1024*sizeof (char));
  
  this->value = malloc (sizeof (struct field_elt_t)); 
  this->value->init = initFieldElt;

  this->left_ctrb = malloc (sizeof (struct wire_contribution_t));
  this->right_ctrb = malloc (sizeof (struct wire_contribution_t));
  this->output_ctrb = malloc (sizeof (struct wire_contribution_t));
 
  this->left_ctrb->init = initCtrb;
  this->right_ctrb->init = initCtrb;
  this->output_ctrb->init = initCtrb;

  this->left_ctrb->init (this->left_ctrb);
  this->right_ctrb->init (this->right_ctrb);
  this->output_ctrb->init (this->output_ctrb);
  
  this->V = NULL;
  this->W = NULL;
  this->Y = NULL;
  
}

//================================
//================================

void viewWire (struct wire_t *this) {
  gmp_printf (" %d %s %Zd %#Zx\n", this->id, this->type, this->value, this->value);
}

//================================
//================================

void viewWireCtrb (struct wire_t *this) {
  struct node_t *node;
  printf (" wire[%d]", this->id);
  if (this->left_ctrb->length>0) {
    node = this->left_ctrb->head;
    while (node!=NULL) {
      gmp_printf(" %Zd", node->gate->root->value);
      node = node->next;
    } 
  }
  else {
    printf(" no left ctrb.");
  }
  
  if (this->right_ctrb->length>0) {
    node = this->right_ctrb->head;
    while (node!=NULL) {
      gmp_printf(" %Zd",node->gate->root->value);
      node = node->next;
    } 
  }
  else {
    printf(" no right ctrb.");
  }

  if (this->output_ctrb->length>0) {
    node = this->output_ctrb->head;
    while (node!=NULL) {
      gmp_printf(" %Zd",node->gate->root->value);
      node = node->next;
    } 
  }
  else {
    printf(" no output ctrb.");
  }
  printf ("\n");
  
  printf (" V(x) = ");
  this->V->view (this->V);
  printf (" W(x) = ");
  this->W->view (this->W);
  printf (" Y(x) = ");
  this->Y->view (this->Y);
  
}

//================================
//================================

void freeWire (struct wire_t *this) {

  if (this!=NULL) {

    // type
    if (this->type!=NULL) {
      free (this->type);
    }

    // value
    if (this->value!=NULL) {
      this->value->free (this->value);
    }
    
    // left ctrb
    if (this->left_ctrb->length==0) {
      free (this->left_ctrb);
    }
    else {
      this->left_ctrb->free (this->left_ctrb);
    }

    // right ctrb
    if (this->right_ctrb->length==0) {
      free (this->right_ctrb);
    }
    else {
      this->right_ctrb->free (this->right_ctrb);
    }

    // output ctrb
    if (this->output_ctrb->length==0) {
      free (this->output_ctrb);
    }
    else {
      this->output_ctrb->free (this->output_ctrb);
    }

    // V
    if (this->V!=NULL) {
      if (this->V->length==0) {
	free (this->V);
      }
      else {
	this->V->free (this->V);
      }
    }
    
    // W
    if (this->W!=NULL) {
      if (this->W->length==0) {
	free (this->W);
      }
      else {
	this->W->free (this->W);
      }
    }
    
    // Y
    if (this->Y!=NULL) {
      if (this->Y->length==0) {
	free (this->Y);
      }
      else {
	this->Y->free (this->Y);
      }
    }
    
    free (this);
  }

}

//================================
//================================

void ctrb2gate (struct wire_contribution_t *this, struct gate_t *gate) {

  struct node_t *node = malloc (sizeof (struct node_t));
  
  node->gate = gate;
  node->next = NULL;
  
  if (this->length==0) {
    this->head = this->tail = node;
  }
  else {
    this->tail->next = node;
    this->tail = node;
  }

  this->length++;

}
