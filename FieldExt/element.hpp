#ifndef _ELEMENT_HPP_
#define _ELEMENT_HPP_

#include <gmp.h> 

#include "field.hpp"

class CElement {
  
public:
  void getFieldCharacteristic (mpz_t);

  void setValue (mpz_t);
  void getValue (mpz_t); 
  
  CElement ();
  CElement (CField*);
  ~CElement ();

private:
  CField* field; //declare a pointer variable called field pointing to an object of type CField

  mpz_t value;
};    

#endif
