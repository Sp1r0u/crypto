#ifndef _WIRE_H_
#define _WIRE_H_

#include <gmp.h>

struct wire_t {
  unsigned int id;
  mpz_t value;
  char *type; //type: input, intermediate, output wires
};

#endif /* _WIRE_H_ */
