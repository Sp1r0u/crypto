#ifndef _CIRCUIT_H_
#define _CIRCUIT_H_

#include <stdlib.h>

#include "wire.h"
#include "gate.h"
#include "field.h"

struct circuit_t {
  unsigned int nWire; //number of wires
  unsigned int nGate; //number of gates
  unsigned int nMul_gate; //number of mult. gates
  unsigned int nAdd_gate; //number of add. gates
  struct wire_t *wire; //all of the wires in the circuit
  //struct gate_t *mul_gate; //all of the mult. gates in the circuit
  //struct gate_t *add_gate; //all of the add. gates in the circuit
  struct gate_t **gate; //all of the gates in the circuit
  unsigned int *index; //index of all the mult. gates in the circuit
};

void init_circuit( struct circuit_t*, struct field_t*, FILE*, gmp_randstate_t );
void display_circuit( struct circuit_t* );
void display_mul_gate( struct circuit_t* ); //show a subset of the circuit, i.e. mult. gates only 
void set_random_input_values( struct circuit_t*, struct field_t*, gmp_randstate_t );
void eval_circuit( struct circuit_t* );

#endif /* _CIRCUIT_H_ */
