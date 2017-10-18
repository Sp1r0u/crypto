#include <stdio.h>
#include <string.h>
#include <gmp.h>

#include "gate.h"

void eval_gate( struct gate_t *g ) {
  //gmp_printf("\nevaluating %s gate %Zd:\n",g->type,g->label);
  if (strcmp(g->type,"add")==0) {
    mpz_add( g->output->value, g->linput->value, g->rinput->value );
    mpz_mod( g->output->value, g->output->value, g->field->p );
  }
  else if (strcmp(g->type,"mul")==0){
    mpz_mul( g->output->value, g->linput->value, g->rinput->value );
    mpz_mod( g->output->value, g->output->value, g->field->p );
  }
  else {
    printf("error: gate type %s not defined... exiting now!\n",g->type);
    exit(1);
  }
}
