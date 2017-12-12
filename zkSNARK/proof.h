#ifndef _PROOF_H_
#define _PROOF_H_

#include <gmp.h>

#include "BN.h"
#include "CRS.h"

struct proof_t {

  element_t arvVP; //a_i * r_v * v_i(s) * P (i=IO wire IDs)
  element_t aalpha_vrvVP; //a_i * alpha_v * r_v * v_i(s) * P  (i=Int wire IDs)

  element_t arwWQ; //a_i * r_w * w_i(s) * Q  (i=Int wire IDs)
  element_t aalpha_wrwWQ; //a_i * alpha_w * r_w * w_i(s) * Q  (i=Int wire IDs)

  element_t aryYP; //a_i * r_y * y_i(s) * P  (i=Int wire IDs)
  element_t aalpha_yryYP; //a_i * alpha_y * r_y * y_i(s) * P  (i=Int wire IDs)

  element_t ar_xbetaVP; //a_i * ( r_v * beta * v_i(s) + r_w * beta * w_i(s) + r_y * beta * y_i(s) ) * P  (i=Int wire IDs)

  struct poly_ring_t *H; //univariate polynomial ring, H(x)*T(x)=P(x) 
  
  void (*init)     (struct proof_t*, struct curve_t*);
  void (*view)     (struct proof_t*);
  void (*free)     (struct proof_t*);
  
  //void (*generate) (struct proof_t*, struct curve_t*, struct crs_t*, struct circuit_t*);
  void (*genPI) (struct proof_t*, struct curve_t*, struct crs_t*,
		 struct circuit_t*, struct field_t*);

  void (*genHPoly) ();
  void (*genTPoly) ();
  void (*genPPoly) ();
};

void init_proof (struct proof_t*, struct curve_t*);

void free_proof (struct proof_t*);
void view_proof (struct proof_t*);

void generate_proof (struct proof_t*, struct curve_t*, struct crs_t*,
		     struct circuit_t*, struct field_t*);

void generate_PPoly (struct proof_t*, struct circuit_t*, struct field_t*);

#endif /* _PROOF_H_ */
