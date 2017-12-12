#ifndef _PROOF_H_
#define _PROOF_H_

#include "CRS.h"

struct proof_t {

  struct crs_t *CRS;
  
  void (*init) (struct proof_t*, struct crs_t*);
  void (*view) (struct proof_t*);
  void (*free) (struct proof_t*);

  void (*verif) (struct proof_t*);
  
  void (*genP) (struct proof_t*);
  void (*genH) (struct proof_t*);
  
  element_t aRvVP; // sum_i ( a_i * r_v * v_i(s) * P ) (i=Int wire IDs)
  element_t aRwWQ; // sum_i ( a_i * r_w * w_i(s) * Q ) (i=Int wire IDs)
  element_t aRwWP; // sum_i ( a_i * r_w * w_i(s) * P ) (i=Int wire IDs)
  element_t aRyYP; // sum_i ( a_i * r_y * y_i(s) * P ) (i=Int wire IDs)
  
  element_t aALPHAvRvVP; // sum_i ( a_i * alpha_v * r_v * v_i(s) * P ) (i=Int wire IDs)
  element_t aALPHAwRwWQ; // sum_i ( a_i * alpha_w * r_w * w_i(s) * Q ) (i=Int wire IDs)
  element_t aALPHAwRwWP; // sum_i ( a_i * alpha_w * r_w * w_i(s) * P ) (i=Int wire IDs)

  element_t aALPHAyRyYP; // sum_i ( a_i * alpha_y * r_y * y_i(s) * P ) (i=Int wire IDs)

  // sum_i ( a_i * ( r_v * beta * v_i(s) + r_w * beta * w_i(s) + r_y * beta * y_i(s) ) * P ) (i=Int wire IDs)
  element_t aRxBETAXP;

  // sum_i ( h_i * s^i * Q )
  element_t hSQ;
  
  // P(x)= sum_i (a_i * v_i(x)) * (a_i * w_i(x)) - (a_i * y_i(x)) 
  struct poly_t *P; //univariate polynomial

  // H(x) = P(x) / T(x);
  struct poly_t *H; //univariate polynomial
  
};

//================================
//================================

void initProof (struct proof_t*, struct crs_t*);
void viewProof (struct proof_t*);
void freeProof (struct proof_t*);

void verifyProof (struct proof_t*);

void genPPoly (struct proof_t*);
void genHPoly (struct proof_t*);
  

#endif /* _PROOF_H_ */
