#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "gate.h"

void eval_gate( struct gate_t *g ) {

  if (strcmp(g->type,"add")==0) {
    gmp_printf(" evaluating %s. gate\n",g->type);
    add_field_elt( g->linput->value, g->rinput->value, g->output->value );
  }
  
  else if (strcmp(g->type,"mul")==0){
    gmp_printf(" evaluating %s. gate %Zd\n",g->type,g->root->value);
    mul_field_elt( g->linput->value, g->rinput->value, g->output->value );
  }

  else {
    printf("error: gate type %s not defined... exiting now!\n",g->type);
    exit(1);
  }
  
}

//==============================================
//==============================================

void free_gate( struct gate_t *gate ) {

  if (gate!=NULL) {

    if (gate->type!=NULL) {
      free (gate->type);
    }

    if (gate->root!=NULL) {
      free_field_elt (gate->root);
    }

    if (gate->lagrange_poly!=NULL ) {
      free_poly_ring (gate->lagrange_poly->head);
      free (gate->lagrange_poly);
    }
    /* 
    if (gate->linput!=NULL) {
      free_wire (gate->linput);
    }
   
    if (gate->rinput!=NULL) {
      free_wire (gate->rinput);
    }

    if (gate->output!=NULL) {
      free_wire (gate->output); 
    }
    */
    free (gate);
    
    return;
  }

}
