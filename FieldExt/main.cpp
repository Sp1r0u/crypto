#include <iostream>
#include <cstdlib>
#include <fstream>
#include <inttypes.h>
#include <vector>
#include <string>
//#include <gmp.h>
#include <time.h>

#include <boost/program_options.hpp>

#include "config.hpp"
#include "field.hpp"
#include "element.hpp"

//https://download.cosine.nl/gvacanti/parsing_configuration_files_c++_CVu265.pdf
//https://stackoverflow.com/questions/30376601/valgrind-memory-still-reachable-with-trivial-program-using-iostream


int main (int argc, char* argv[]) {

  //struct config_t* cfg = (struct config_t*) malloc (sizeof (struct config_t) );
  config_t cfg;

  namespace po = boost::program_options;

  po::options_description desc;

  desc.add_options ()
    ("*FIELD_CHARACTERISTIC", po::value<std::string>(&cfg.p));
    //("*FIELD_CHARACTERISTIC", po::value<std::string>(&cfg->p));
    
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

  std::cout << " address Fp " << &Fp << std::endl;
  
  CElement elt (&Fp);
  mpz_t dummy;
  elt.getFieldCharacteristic (dummy);
  gmp_printf ("dummy %Zd\n", dummy);
  mpz_clear (dummy);
    
  mpz_t foo;
  Fp.getCharacteristic (foo);
  gmp_printf(" foo %Zd\n", foo);
  mpz_clear (foo);
  
  //free (cfg);
  
  return (EXIT_SUCCESS);
  
}
