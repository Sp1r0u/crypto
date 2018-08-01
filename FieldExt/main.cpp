#include <iostream>
#include <cstdlib>
#include <fstream>
#include <inttypes.h>
#include <vector>
#include <string>
#include <time.h>

#include <boost/program_options.hpp>

#include "config.hpp"
#include "field.hpp"
#include "element.hpp"

//https://download.cosine.nl/gvacanti/parsing_configuration_files_c++_CVu265.pdf
//https://stackoverflow.com/questions/30376601/valgrind-memory-still-reachable-with-trivial-program-using-iostream


int main (int argc, char* argv[]) {

  config_t cfg;

  namespace po = boost::program_options;

  po::options_description desc;

  desc.add_options ()
    ("*FIELD_CHARACTERISTIC", po::value<std::string>(&cfg.p));
    
  po::variables_map vm;

  const bool allow_unregistered = true;

  po::store (po::parse_config_file (std::cin, desc, allow_unregistered), vm);

  po::notify(vm);

  // initialize the random state with default algorithm
  gmp_randstate_t RND_STATE;
  gmp_randinit_default (RND_STATE);
  
  // seed the state
  time_t seed;
  time (&seed);
  gmp_randseed_ui (RND_STATE, seed);
  
  CField Fp (cfg); 

  /**************** test here *******************/
  uint16_t i;
  uint16_t deg = 2;
  
  std::pair<uint16_t, mpz_t>             mypair[deg+1];
  std::list <std::pair<uint16_t, mpz_t>> foo;
  
  for (i=0; i<=deg; i++) {

    mypair[i].first  = i;

    mpz_init (mypair[i].second);
    mpz_set_ui (mypair[i].second, 999-i);

    foo.push_front (mypair[i]);
  }

  for (i=0; i<=deg; i++) {
    mpz_clear (mypair[i].second);
  }
  
  /**************** test here *******************/

  
  //Fp.buildIrreduciblePoly (2, RND_STATE);
  
  //std::cout << " address Fp " << &Fp << std::endl;
  
  /*
  CElement elt (&Fp);
  elt.setRndElement (&elt, RND_STATE);
  mpz_t echo;
  elt.getValue (echo);
  gmp_printf ("echo %Zd\n", echo);
  mpz_clear (echo);
  mpz_t dummy;
  elt.getFieldCharacteristic (dummy);
  gmp_printf ("dummy %Zd\n", dummy);
  mpz_clear (dummy);
    
  mpz_t foo;
  Fp.getCharacteristic (foo);
  gmp_printf(" foo %Zd\n", foo);
  mpz_clear (foo);
  */

  gmp_randclear (RND_STATE);
  
  return (EXIT_SUCCESS);
  
}
