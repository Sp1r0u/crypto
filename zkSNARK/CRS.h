#ifndef _CRS_H_
#define _CRS_H_

#include <stdbool.h>
#include <gmp.h>

#include "circuit.h"
#include "BN.h"
#include "field.h"

extern bool fDebug;
extern unsigned int fOrder;

#define NPARAMS 9

struct extrakey_t {
  struct wire_t *wire;
  element_t elt;
};

struct EKey_t {
  struct extrakey_t **rvVP; //r_v * v_i(s) * P  (i=Int wire IDs)
  struct extrakey_t **rwWQ; //r_w * w_i(s) * Q  (i=Int wire IDs)
  struct extrakey_t **ryYP; //r_y * y_i(s) * P  (i=Int wire IDs)
  struct extrakey_t **alpha_vrvVP; //alpha_v * r_v * v_i(s) * P  (i=Int wire IDs)
  struct extrakey_t **alpha_wrwWQ; //alpha_w * r_w * w_i(s) * Q  (i=Int wire IDs)
  struct extrakey_t **alpha_yryYP; //alpha_y * r_y * y_i(s) * P  (i=Int wire IDs)
  struct extrakey_t **r_xbetaVP; //( r_v * beta * v_i(s) + r_w * beta * w_i(s) + r_y * beta * y_i(s) ) * P  (i=Int wire IDs)
};

struct VKey_t {
  element_t *P; //<P>=G1
  element_t *Q; //<Q>=G2
  element_t alpha_vQ; //alpha_v * Q
  element_t alpha_wQ; //alpha_w * Q
  element_t alpha_wP; //alpha_w * P
  element_t alpha_yQ; //alpha_y * P
  element_t betaP; //beta * P
  element_t betaQ; //beta * Q
  element_t rytP;  //r_y * t(s) * P
  //element_t **rvVP; //r_v * v_i(s) * P  (i=IO wire IDs)
  //element_t **rwWQ; //r_w * w_i(s) * Q  (i=IO wire IDs)
  //element_t **ryYP; //r_y * y_i(s) * P  (i=IO wire IDs)
  struct extrakey_t **rvVP; //r_v * v_i(s) * P  (i=IO wire IDs)
  struct extrakey_t **rwWQ; //r_w * w_i(s) * Q  (i=IO wire IDs)
  struct extrakey_t **ryYP; //r_y * y_i(s) * P  (i=IO wire IDs)
};

/* param contains a list of random non-zero field element in that order
 * param[0]=rv
 * param[1]=rw
 * param[2]=ry
 * param[3]=s
 * param[4]=alpha_v
 * param[5]=alpha_w
 * param[6]=alpha_y
 * param[7]=beta
 * param[8]=gamma */

struct crs_t {
  struct field_elt_t **param;
  struct EKey_t *EK;
  struct VKey_t *VK;
  unsigned int *nIO_V;
  unsigned int *nIO_W;
  unsigned int *nIO_Y;
  unsigned int *nInt_V;
  unsigned int *nInt_W;
  unsigned int *nInt_Y;
};

void generate_crs( struct crs_t*, struct circuit_t*, struct curve_t*, struct field_t*, gmp_randstate_t );

void generate_EKVK( struct crs_t*, struct circuit_t*, struct curve_t*, struct field_t* );

void generate_target_poly( struct circuit_t*, struct field_t*, struct poly_ring_t* );
void generate_TPoly( struct circuit_t*, struct field_t*, struct poly_ring_t* );

void generate_Lagrange_polynomial( struct circuit_t*, struct field_t* );

void generate_QAP_polynomial( struct circuit_t*, struct field_t*, struct crs_t* );

// returns subset of Vk and Ek elements (e.g. rxXY = {rvV1P, rwW3Q, ryY4P...})
void evaluate_rxXY( struct circuit_t*, struct crs_t* );

void print_EKVK( struct circuit_t*, struct crs_t* );

void free_crs( struct crs_t*, struct circuit_t* );

#endif /* _CRS_H_ */
