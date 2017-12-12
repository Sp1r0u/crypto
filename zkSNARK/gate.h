#ifndef _GATE_H_
#define _GATE_H_

#include <gmp.h>

#include "wire.h"
#include "field.h"

struct gate_t {
  char *type; //type: add, mul, etc.. gate
  unsigned int ninput; //number of inputs
  unsigned int noutput; //number of outputs
  struct field_elt_t *root; //random non-zero element of Fp
  struct poly_ring_t *lagrange_poly; //univariate polynomial ring  
  struct wire_t *linput; //left input wire
  struct wire_t *rinput; //right input wire
  struct wire_t *output; //output wire
};

void eval_gate( struct gate_t* );

void free_gate( struct gate_t* );

#endif /* _GATE_H_ */
