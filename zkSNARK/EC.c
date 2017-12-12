#include <stdlib.h>
#include <string.h>

#include "EC.h"

void initEC (struct curve_t *this) {

  mpz_init (this->x);
  mpz_init (this->tx);
  mpz_init (this->nx);
  mpz_init (this->px);
  mpz_init (this->r);

  this->view = viewEC;
  this->free = freeEC;
  
  mpz_set_str (this->x, "1", 10);
  mpz_neg (this->x, this->x);

  mpz_set_str (this->tx, "1", 10);
  mpz_neg (this->tx, this->tx);
  
  mpz_set_str (this->px, "205523667896953300194896352429254920972540065223", 10);
  mpz_set_str (this->nx, "205523667896953300194895899082072403858390252929", 10);

  char str1[1024] = "";
  char str2[1024] = "";

  // type of the EC
  strcat (str1, "type f\n");

  // field characteristic
  strcat (str1, "q ");
  mpz_get_str (str2, 10, this->px); 
  strcat (str1, str2);
  strcat (str1, "\n");
  
  // order of the curve
  strcat (str1, "r ");
  mpz_get_str (str2, 10, this->nx);
  strcat (str1, str2);
  strcat (str1, "\n");

  // b
  strcat (str1, "b 40218105156867728698573668525883168222119515413\n");

  // beta
  strcat (str1, "beta 115334401956802802075595682801335644058796914268\n");
  
  // alpha0
  strcat (str1, "alpha0 191079354656274778837764015557338301375963168470\n");
  
  // alpha1
  strcat (str1, "alpha1 71445317903696340296199556072836940741717506375");

  pbc_param_init_set_str (this->pbc_params, str1);

  FILE *fp=fopen ("pbc_params.txt", "w");

  if (fp==NULL) {
    printf("error while opening file... exiting now\n");
    exit (1);
  }

  pbc_param_out_str (fp, this->pbc_params);

  pairing_init_pbc_param (this->pairing, this->pbc_params);

  element_init_G1 (this->P, this->pairing);
  element_init_G2 (this->Q, this->pairing);
 
  element_random (this->P);
  element_random (this->Q);

  // read pbc_pairing.h L18
  mpz_set (this->r, this->pairing->r);
  
  this->view (this);
  
  fclose (fp);
}

//================================
//================================

void viewEC (struct curve_t *this) {
  printf ("\n");
  printf (" **** EC Settings ****\n");
  gmp_printf(" x  = %Zd\n", this->x);
  gmp_printf(" tx = %Zd\n", this->tx);
  gmp_printf(" px = %Zd\n", this->px);
  gmp_printf(" nx = %Zd\n", this->nx);

  if (pairing_is_symmetric (this->pairing)==0) {
    printf ("\n");
    printf(" Sanety check: asymmetric pairing function\n");
  }
  else {
    printf ("\n");
    printf (" Sanety check: symmetric pairing function..., exiting now!\n");
    exit(1);
  }

  printf (" Length in bytes to represent an element of G1: %d\n", pairing_length_in_bytes_G1 (this->pairing));
  printf (" Length in bytes to represent an element of G2: %d\n", pairing_length_in_bytes_G2 (this->pairing));

  printf ("\n");
  element_printf (" P = %B\n", this->P);
  element_printf (" Q = %B\n", this->Q);
  
  gmp_printf(" r = %Zd\n", this->r);
  
}

//================================
//================================

void freeEC (struct curve_t *this) {
  mpz_clear (this->x);
  mpz_clear (this->tx);
  mpz_clear (this->nx);
  mpz_clear (this->px);
  mpz_clear (this->r);

  element_clear (this->P);
  element_clear (this->Q);

  pbc_param_clear (this->pbc_params);

  pairing_clear (this->pairing);  
    
  free (this);
  
}

//================================
//================================

