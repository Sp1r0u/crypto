#ifndef _CRS_H_
#define _CRS_H_

#include <gmp.h>

#include "circuit.h"
#include "BN.h"
#include "field.h"

#define NPARAMS 9

struct EKey_t {
};

struct VKey_t {
  element_t *P; //<P>=G1
  element_t *Q; //<Q>=G2
  element_t alpha_vQ; //alpha_v * Q
  element_t alpha_wQ; //alpha_w * Q
  element_t alpha_wP; //alpha_w * P
  element_t betaP; //beta * P
  element_t betaQ; //beta * Q
  element_t rytP;  //r_y * t(s) * P
  element_t *rvVP; //r_v * v_i(s) * P  (i=IO wire IDs)
  element_t *rwWQ; //r_w * w_i(s) * Q  (i=IO wire IDs)
  element_t *ryYP; //r_y * y_i(s) * P  (i=IO wire IDs)
};

struct crs_t {
  //param contains a list of random non-zero field element in that order
  //param[0]=rv
  //param[1]=rw
  //param[2]=ry
  //param[3]=s
  //param[4]=alpha_v
  //param[5]=alpha_w
  //param[6]=alpha_y
  //param[7]=beta
  //param[8]=gamma
  mpz_t param[NPARAMS];
  struct EKey_t *EK;
  struct VKey_t *VK;
};

void generate_crs( struct crs_t*, struct circuit_t*, struct curve_t*, struct field_t*, gmp_randstate_t );
void generate_EKey( struct EKey_t*, struct circuit_t* );
void generate_VKey( struct VKey_t*, struct circuit_t*, struct curve_t*, struct field_t*, mpz_t[] );
void evaluate_target_poly( struct circuit_t*, mpz_t, mpz_t, struct field_t* );

#endif /* _CRS_H_ */
