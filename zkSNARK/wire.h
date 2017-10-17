#ifndef _WIRE_H_
#define _WIRE_H_

struct wire_t {
  unsigned int id;
  unsigned int value;
  char *type; //type: input, intermediate, output wires
};

#endif /* _WIRE_H_ */
