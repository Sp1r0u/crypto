#ifndef _ELEMENT_H_
#define _ELEMENT_H_

#include "field.h"

struct element_node_t {
  mpz_t coeff;
  int exp;
  struct element_node_t *next;
};

struct element_t {
  int length;
  struct element_node_t *head;
  struct element_node_t *tail;
  struct field_t *field;
};

void set_element( struct element_t*, mpz_t, unsigned int );
void print_element( struct element_t* ); 
void multiply_elements( struct element_t*, struct element_t*, struct element_t* );
void add_elements( struct element_t*, struct element_t*, struct element_t* );

int get_element_length( struct element_t* );
int get_element_highest_degree( struct element_t* );

#endif /* _ELEMENT_H_ */
