#include "field.hpp"

CField::CField (struct config_t cfg) {
  if (!(cfg.p).empty())
    setCharacteristic (cfg.p);
}

//========================================================

CField::~CField () {
  mpz_clear (p);
}

//========================================================

void CField::setCharacteristic (std::string str) {
  mpz_init_set_str (p, str.c_str(), 10);
}

//========================================================

void CField::getCharacteristic (mpz_t out) {
  mpz_init_set (out, p);
}

//========================================================

void CField::getRndElement (CElement* ptr, gmp_randstate_t state) {
  mpz_t p, value;
  ptr->getFieldCharacteristic (p);
  ptr->getValue (value);
  while (mpz_cmp_ui (value, 0)==0) {
    mpz_urandomm (value, state, p);
  }

  gmp_printf ("getRndElemt %Zd\n", value);
  ptr->setValue (value);
  
  mpz_clear (p);
  mpz_clear (value);
}
