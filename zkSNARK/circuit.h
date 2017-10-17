#ifndef _CIRCUIT_H_
#define _CIRCUIT_H_

#include <stdlib.h>

#include "wire.h"
#include "gate.h"
#include "field.h"

struct circuit_t {
  unsigned int nWire; //number of wires
  unsigned int nMul_gate; //number of mult. gates
  struct wire_t *wire; //all of the wires in the circuit
  struct gate_t *mul_gate; //all of the mult. gates in the circuit (a subset of all_gate)
  struct gate_t *all_gate; //all gates in the circuit 
};

void init_circuit( struct circuit_t*, struct field_t*, FILE*, gmp_randstate_t );
//void init_circuit( struct circuit_t*, FILE*, gmp_randstate_t );

#endif /* _CIRCUIT_H_ */
