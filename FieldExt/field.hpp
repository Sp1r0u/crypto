#ifndef FIELD_HPP
#define FIELD_HPP

#include <iostream>
#include <list>
#include <utility>
#include <iterator>
#include <inttypes.h>

#include <gmp.h>

#include "config.hpp"
#include "element.hpp"

class CField {

public:
  void setCharacteristic (std::string);
  void getCharacteristic (mpz_t);

  void setRndMonicPoly      (uint16_t, gmp_randstate_t);
  void buildIrreduciblePoly (uint16_t, gmp_randstate_t); // generate an irreducible polynomial 
  void showPoly             (std::list <std::pair<uint16_t, mpz_t>>);
  
  void getRndElt (mpz_t, gmp_randstate_t);
  
  void getRndElement (class CElement*, gmp_randstate_t);
  
  CField ();
  CField (struct config_t);
  ~CField ();

private:
  mpz_t p;
  
  std::list <std::pair<uint16_t, mpz_t>> P; // irreducible polynomial
  
};

#endif
