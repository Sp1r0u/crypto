#ifndef _WIRE_H_
#define _WIRE_H_

#include <gmp.h>

#include "field.h"
#include "gate.h"

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
  char *type; //type: input, intermediate, output wires
  struct field_elt_t *value; //an element of Fp
  struct wire_contribution_t *left_ctrb;
  struct wire_contribution_t *right_ctrb;
  struct wire_contribution_t *output_ctrb;
  struct poly_ring_t *V; //univariate polynomial ring
  struct poly_ring_t *W; //univariate polynomial ring
  struct poly_ring_t *Y; //univariate polynomial ring
};

void init_wire_contribution( struct wire_t* );

void wire_contribution_to_mul_gate( struct wire_contribution_t*, struct gate_t* );

void print_wire_contribution_to_mul_gate( struct wire_t*);

void free_wire( struct wire_t* );

void free_wire_contribution( struct node_t* );

#endif /* _WIRE_H_ */
