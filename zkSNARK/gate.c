#include <stdlib.h>

#include "gate.h"

#define strSize 20

void initGate (struct gate_t *this, char str[strSize]) {
  this->view = viewGate;
  this->free = freeGate;
  this->eval = evalGate;

  this->type = malloc (1024*sizeof (char));

  if (strcmp(str,"mul")==0) {
    this->root = malloc (sizeof (struct field_elt_t)); 
    this->root->init = initFieldElt;
  }
  else {
    this->root = NULL;
  }

  this->lagrange_poly = NULL;
}

//================================
//================================

void viewGate (struct gate_t *this) {
  if (strcmp(this->type,"mul")==0) {
    gmp_printf(" %s %Zd %#Zx %d %d %d %d %d\n",
	       this->type, this->root->value, this->root->value, this->ninput, this->linput->id, this->rinput->id, this->noutput, this->output->id);
  }
  else {
    gmp_printf(" %s NULL %d %d %d %d %d\n",
	       this->type, this->ninput, this->linput->id, this->rinput->id, this->noutput, this->output->id);  
  }
}

//================================
//================================

void freeGate (struct gate_t *this) {

  if (this!=NULL) {

    // type
    if (this->type!=NULL) {
      free (this->type);
    }

    // root
    if (this->root!=NULL) {
      this->root->free (this->root);
    }

    // Lagrange
    if (this->lagrange_poly!=NULL ) {
      this->lagrange_poly->free (this->lagrange_poly);
    }
    
    free (this);
  }
}

//================================
//================================

void evalGate (struct gate_t *this) {
  
  if (strcmp(this->type,"add")==0) {
    this->output->value->field->addElt
      (this->linput->value, this->rinput->value, this->output->value);
  }
  
  else if (strcmp(this->type,"mul")==0){
    this->output->value->field->mulElt
      (this->linput->value, this->rinput->value, this->output->value);
  }

  else {
    printf ("error in evalGate: gate type %s not defined... exiting now!\n", this->type);
    exit (1);
  }
  
}
