#include "field.hpp"

CField::CField (struct config_t cfg) {
  if (!(cfg.p).empty())
    setCharacteristic (cfg.p);
}

//========================================================
//========================================================

CField::~CField () {
  mpz_clear (p);

  if (!P.empty()) {
    std::list <std::pair<uint16_t, mpz_t>> :: iterator it;
    for (it=P.begin (); it!=P.end(); ++it) {
      mpz_clear (it->second);
    }
  }
}

//========================================================
//========================================================

void CField::buildIrreduciblePoly (uint16_t deg, gmp_randstate_t state) {
  setRndMonicPoly (deg, state);
  
  //std::list <std::pair<uint16_t, mpz_t>> u;
  //std::pair<uint16_t, mpz_t> mypair;
  
  showPoly (P);
}

//========================================================
//========================================================

void CField::showPoly (std::list <std::pair<uint16_t, mpz_t>> irredP) {

  std::list <std::pair<uint16_t, mpz_t>> :: iterator it;

  for (it=irredP.begin (); it!=irredP.end(); ++it) {
    gmp_printf ("%Zdx^%" PRIu16, it->second, it->first);
    if (std::next(it)!=irredP.end())
	printf (" + ");
  }
}

//========================================================
//========================================================

void CField::setRndMonicPoly (uint16_t deg, gmp_randstate_t state) {
  uint16_t i;
  std::pair<uint16_t, mpz_t> mypair[deg+1];
  
  for (i=0; i<deg; i++) {
      
    mypair[i].first = i;
    
    mpz_init (mypair[i].second);
    getRndElt (mypair[i].second, state);
    
    P.push_front (mypair[i]);
  }

  mypair[deg].first = deg;

  mpz_set_ui (mypair[deg].second, 1);
  P.push_front (mypair[deg]);

}

//========================================================
//========================================================

void CField::getRndElt (mpz_t elt, gmp_randstate_t state) {
  while (mpz_cmp_ui (elt, 0)==0) {
    mpz_urandomm (elt, state, p);
  }
}

//========================================================
//========================================================

void CField::setCharacteristic (std::string str) {
  mpz_init_set_str (p, str.c_str(), 10);
}

//========================================================
//========================================================

void CField::getCharacteristic (mpz_t out) {
  mpz_init_set (out, p);
}

//========================================================
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
