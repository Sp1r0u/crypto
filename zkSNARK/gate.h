#ifndef _GATE_H_
#define _GATE_H_

#include <gmp.h>

#include "wire.h"
#include "field.h"

struct gate_t {
  char *type; //type: add, mul, etc.. gate
  mpz_t label; //random non-zero element of Fp
  unsigned int ninput; //number of inputs
  unsigned int noutput; //number of outputs
  struct wire_t *linput; //left input wire
  struct wire_t *rinput; //right input wire
  struct wire_t *output; //output wire
  struct field_t *field;
};

#endif /* _GATE_H_ */
