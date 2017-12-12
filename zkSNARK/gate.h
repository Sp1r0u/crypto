#ifndef _GATE_H_
#define _GATE_H_

#include <string.h>

#include "wire.h"
#include "field.h"
#include "poly.h"

struct gate_t {
  char *type; //type: add, mul, etc.. gate
  unsigned int ninput; //number of inputs
  unsigned int noutput; //number of outputs
  struct field_elt_t *root; //random non-zero element of Fp
  struct poly_t *lagrange_poly; //univariate polynomial ring  
  struct wire_t *linput; //left input wire
  struct wire_t *rinput; //right input wire
  struct wire_t *output; //output wire

  void (*init) (struct gate_t*, char str[]);
  void (*view) (struct gate_t*);
  void (*free) (struct gate_t*);
  void (*eval) (struct gate_t*);
};

//================================
//================================

void initGate (struct gate_t*, char str[]);
void viewGate (struct gate_t*);
void freeGate (struct gate_t*);
void evalGate (struct gate_t*);

#endif /* _GATE_H_ */
