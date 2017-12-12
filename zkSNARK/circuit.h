#ifndef _CIRCUIT_H_
#define _CIRCUIT_H_

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "field.h"
#include "wire.h"
#include "gate.h"

struct circuit_t {

  struct field_t *field;
  
  unsigned int nWire;     // number of wires
  unsigned int nIOWire;   // number of IO wires
  unsigned int nGate;     // number of gates
  unsigned int nMul_gate; // number of mult. gates
  unsigned int nAdd_gate; // number of add. gates
  struct wire_t **wire;   // all of the wires in the circuit
  struct gate_t **gate;   // all of the gates in the circuit
  unsigned int *index;    // index of all the mult. gates in the circuit

  void (*init) (struct circuit_t*, struct field_t*, FILE*, gmp_randstate_t);
  void (*view) (struct circuit_t*);
  void (*free) (struct circuit_t*);
  void (*eval) (struct circuit_t*);
  
  void (*newWire) (struct circuit_t*, unsigned int*, char[], char[]);
  void (*newGate) (struct circuit_t*, unsigned int*, unsigned int*, unsigned int, char[], char[], gmp_randstate_t);

  void (*rndInputs) (struct circuit_t*, gmp_randstate_t);
  
};

//================================
//================================

void initCircuit (struct circuit_t*, struct field_t*, FILE*, gmp_randstate_t);
void viewCircuit (struct circuit_t*);
void freeCircuit (struct circuit_t*);
void evalCircuit (struct circuit_t*);

void newWire (struct circuit_t*, unsigned int*, char[], char[]);
void newGate (struct circuit_t*, unsigned int*, unsigned int*, unsigned int, char[], char[], gmp_randstate_t);

void setRNDInputs (struct circuit_t*, gmp_randstate_t);

#endif /* _CIRCUIT_H_ */
