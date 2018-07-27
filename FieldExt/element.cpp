#include "element.hpp"

void CElement::getFieldCharacteristic (mpz_t out) {
  field->getCharacteristic (out);
}

//========================================================

void CElement::getValue (mpz_t out) {
  mpz_init_set (out, value); 
}

//========================================================

void CElement::setValue (mpz_t in) {
  mpz_set (value, in); 
}

//========================================================

CElement::CElement (CField* ptr) {
  if (ptr)
    field = ptr;

  mpz_init (value);
}

//========================================================

CElement::CElement () {
  field = NULL;
  mpz_init (value);
}

//========================================================

CElement::~CElement () {
  mpz_clear (value);
}

