#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include <gmp.h> 

#include "field.hpp"

class CElement {
  
public:
  void getFieldCharacteristic (mpz_t);
  void setRndElement (class CElement*, gmp_randstate_t);
  
  void setValue (mpz_t);
  void getValue (mpz_t); 
  
  CElement ();
  CElement (class CField*);
  ~CElement ();

private:
  CField* field; //declare a pointer variable called field pointing to an object of type CField

  mpz_t value;
};    

#endif
