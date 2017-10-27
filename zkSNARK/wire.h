#ifndef _WIRE_H_
#define _WIRE_H_

#include "gate.h"

#include <gmp.h>

struct node_t {
  struct gate_t *gate;
  struct node_t *next;
};

struct wire_contribution_t {
  int length;
  struct node_t *head;
  struct node_t *tail;
};

struct wire_t {
  unsigned int id;
  mpz_t value;
  char *type; //type: input, intermediate, output wires
  struct wire_contribution_t *left_ctrb;
  struct wire_contribution_t *right_ctrb;
};

void init_wire_contribution( struct wire_t* );
void wire_contribution_to_mul_gate( struct wire_contribution_t*, struct gate_t* );
void print_wire_contribution_to_mul_gate( struct wire_t*);

#endif /* _WIRE_H_ */
