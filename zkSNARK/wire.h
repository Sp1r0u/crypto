#ifndef _WIRE_H_
#define _WIRE_H_

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

  void (*init) (struct wire_contribution_t*);
  void (*free) (struct wire_contribution_t*);
};

struct wire_t {
  unsigned int id;
  char *type; //type: input, intermediate, output wires
  struct field_elt_t *value; //an element of Fp
  struct wire_contribution_t *left_ctrb;
  struct wire_contribution_t *right_ctrb;
  struct wire_contribution_t *output_ctrb;
  struct poly_t *V; //univariate polynomial
  struct poly_t *W; //univariate polynomial
  struct poly_t *Y; //univariate polynomial

  void (*init) (struct wire_t*); 
  void (*view) (struct wire_t*);
  void (*free) (struct wire_t*);

  void (*ctrb) (struct wire_contribution_t*, struct gate_t*);
  void (*viewCtrb) (struct wire_t*);
};

//================================
//================================

void initWire (struct wire_t*);
void viewWire (struct wire_t*);
void freeWire (struct wire_t*);

void initCtrb (struct wire_contribution_t*);
void freeCtrb (struct wire_contribution_t*);

void ctrb2gate (struct wire_contribution_t*, struct gate_t*);
void viewWireCtrb (struct wire_t*);

#endif /* _WIRE_H_ */
