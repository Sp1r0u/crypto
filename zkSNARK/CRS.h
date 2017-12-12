#ifndef _CRS_H_
#define _CRS_H_

#include <gmp.h>

#include "field.h"
#include "wire.h"
#include "poly.h"
#include "circuit.h"

#define NPARAMS 9

struct extraKey_t {
  struct wire_t *wire;
  element_t elt;
  mpz_t value; // for testing purposes only
};

//================================
//================================

struct EKey_t {
  struct extraKey_t **RvVP; // r_v * v_i(s) * P  (i=Int wire IDs)
  struct extraKey_t **RwWQ; // r_w * w_i(s) * Q  (i=Int wire IDs)
  struct extraKey_t **RyYP; // r_y * y_i(s) * P  (i=Int wire IDs)
  struct extraKey_t **ALPHAvRvVP; // alpha_v * r_v * v_i(s) * P  (i=Int wire IDs)
  struct extraKey_t **ALPHAwRwWQ; // alpha_w * r_w * w_i(s) * Q  (i=Int wire IDs)
  struct extraKey_t **ALPHAwRwWP; // alpha_w * r_w * w_i(s) * P  (i=Int wire IDs)
  struct extraKey_t **ALPHAyRyYP; // alpha_y * r_y * y_i(s) * P  (i=Int wire IDs)
  struct extraKey_t **RxBETAXP;   // (r_v * beta * v_i(s) + r_w * beta * w_i(s) + r_y * beta * y_i(s)) * P  (i=Int wire IDs)
  struct extraKey_t **SQ; // s^j * Q (j=0,..., deg(t))
};

//================================
//================================

struct VKey_t { 
  element_t *P; // <P>=G1
  element_t *Q; // <Q>=G2
  element_t ALPHAvQ; // alpha_v * Q
  element_t ALPHAwQ; // alpha_w * Q
  element_t ALPHAwP; // alpha_w * P
  element_t ALPHAyQ; // alpha_y * Q
  element_t BETAP;   // beta * P
  element_t BETAQ;   // beta * Q
  element_t RyTP;    // r_y * t(s) * P

  struct extraKey_t **RvVP; // r_v * v_i(s) * P  (i=IO wire IDs)
  struct extraKey_t **RwWQ; // r_w * w_i(s) * Q  (i=IO wire IDs)
  struct extraKey_t **RyYP; // r_y * y_i(s) * P  (i=IO wire IDs)

};

//================================
//================================

/* param contains a list of random non-zero field element in that order
 * param[0]=r_v
 * param[1]=r_w
 * param[2]=r_y
 * param[3]=s
 * param[4]=alpha_v
 * param[5]=alpha_w
 * param[6]=alpha_y
 * param[7]=beta
 * param[8]=gamma */

struct crs_t {
  struct circuit_t *C;

  struct field_elt_t **param;
  
  struct EKey_t *EK;
  struct VKey_t *VK;

  unsigned int *nIO_V;
  unsigned int *nIO_W;
  unsigned int *nIO_Y;
  unsigned int *nInt_V;
  unsigned int *nInt_W;
  unsigned int *nInt_Y;

  struct poly_t *T; // target poly
  
  void (*init) (struct crs_t*, struct circuit_t*, gmp_randstate_t);
  void (*free) (struct crs_t*);
  void (*view) (struct crs_t*);

  void (*Lagrange) (struct crs_t*);

  void (*genQAP) (struct crs_t*);
  void (*genQAPP) (struct poly_t*, struct crs_t*, struct node_t*,
		   struct wire_t*, unsigned int*, unsigned int*);

  void (*genEKVK) (struct crs_t*);
  
  void (*genT) (struct poly_t*, struct circuit_t*);
  
  void (*genRxXY) (struct crs_t*);

};

//================================
//================================

void initCRS (struct crs_t*, struct circuit_t*, gmp_randstate_t);
void freeCRS (struct crs_t*);
void viewCRS (struct crs_t*);

void generateLagrangePoly (struct crs_t*);

void generateQAP (struct crs_t*);

void generateQAPPoly
(struct poly_t*, struct crs_t*, struct node_t*, struct wire_t*, unsigned int*, unsigned int*);

void generateEKeyVKey (struct crs_t*);

void generateTargetPoly (struct poly_t*, struct circuit_t*);

void generateRxXY (struct crs_t*);

#endif // _CRS_H_
