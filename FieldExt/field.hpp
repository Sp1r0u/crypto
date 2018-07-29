#ifndef FIELD_HPP
#define FIELD_HPP

#include <iostream>
#include <gmp.h>

#include "config.hpp"
#include "element.hpp"

class CField {

public:
  void setCharacteristic (std::string);
  void getCharacteristic (mpz_t);

  void setIrreduciblePoly     ();
  void getIrreduciblePoly     ();

  void getRndElement (class CElement*, gmp_randstate_t);
  
  CField ();
  CField (struct config_t);
  ~CField ();

private:
  mpz_t p; 

};

#endif
