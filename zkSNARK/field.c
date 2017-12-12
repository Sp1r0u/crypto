#include <stdlib.h>

#include "field.h"

void initFieldElt (struct field_elt_t *this, struct field_t *Fp) {
  mpz_init (this->value);

  this->field = Fp;
    
  this->view = viewFieldElt;
  this->free = freeFieldElt;
}

//================================
//================================

void initField (struct field_t *this, struct curve_t *EC) {
  this->EC = EC;
  this->free = freeField;
  
  this->rndElt = rndFieldElt;
  this->modElt = modFieldElt;
  this->invElt = invFieldElt;
  this->mulElt = mulFieldElt;
  this->divElt = divFieldElt;
  this->addElt = addFieldElt;
  this->subElt = subFieldElt;
  this->negElt = negFieldElt;

  this->powuiElt = powuiFieldElt;
}

//================================
//================================

void rndFieldElt (struct field_elt_t *this, gmp_randstate_t state) {

  if (_debug==false) {
    while (mpz_cmp_ui(this->value, 0)==0) {
      mpz_urandomm (this->value, state, this->field->EC->r);
    }
  }
  
  // bound for mpz_urandomm
  else {
    mpz_t upBnd;
    mpz_init_set_ui (upBnd, _order);
    while (mpz_cmp_ui (this->value, 0)==0) {
      mpz_urandomm (this->value, state, upBnd);
    }
    mpz_clear (upBnd);
  }
}

//================================
//================================

void freeField (struct field_t *this) {
  if (this!=NULL) {
    free (this);
  }
}

//================================
//================================

void freeFieldElt (struct field_elt_t *this) {

  if (this!=NULL) {
    mpz_clear (this->value);
    free (this);
  }

}

//================================
//================================

void viewFieldElt (struct field_elt_t *this) {
  if (this!=NULL) {
    gmp_printf (" elt %Zd\n",this->value);
  }
}

//================================
//================================

void modFieldElt (struct field_elt_t *this) {

  if (_debug==false) {
    mpz_mod (this->value, this->value, this->field->EC->r );
  }
  
  else {
    mpz_t mod;
    mpz_init_set_ui (mod, _order);
    mpz_mod (this->value, this->value, mod);
    mpz_clear (mod);
  }
}

//================================
//================================

void invFieldElt (struct field_elt_t *out, struct field_elt_t *in) {
  
  if (in->field!=out->field) {
    printf (" %p %p", in->field, out->field);
    printf ("error in invFieldElt: in and out are not elements of the same field... exiting now\n!");
    exit (1);
  }

  if (_debug==false) {
    if (mpz_invert (out->value, in->value, out->field->EC->r)==0) {
      printf("inverse does not exist... exiting now!");
      exit (1);
    }
  }
  else {
    mpz_t mod;
    mpz_init_set_ui (mod, _order);
    if (mpz_invert (out->value, in->value, mod)==0) {
      printf("inverse does not exist... exiting now!");
      exit (1);
    }
    mpz_clear( mod );
  }
}

//================================
//================================

void mulFieldElt (struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out) {

  if (in1->field!=in2->field) {
    printf (" %p %p", in1->field,in2->field);
    printf ("error in mulFieldElt: in1 and in2 are not elements of the same field... exiting now\n!");
    exit (1);
  }

  mpz_mul (out->value, in1->value, in2->value);
  out->field->modElt (out);

}

//================================
//================================

void divFieldElt (struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out) {

  if (in1->field!=in2->field) {
    printf (" %p %p", in1->field,in2->field);
    printf("error in divFieldElt: in1 and in2 are not elements of the same field... exiting now!\n");
    exit (1);
  }
  
  struct field_elt_t *tmp = malloc (sizeof (struct field_elt_t));
  tmp->init = initFieldElt;
  tmp->init (tmp, in1->field);

  tmp->field->invElt (tmp, in2);//temp<-in2^-1
  
  out->field->mulElt (in1, tmp, out);

  tmp->free (tmp);
}

//================================
//================================

void addFieldElt (struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out) {

  if (in1->field!=in2->field) {
    printf (" %p %p", in1->field,in2->field);
    printf ("error in addFieldElt: in1 and in2 are not elements of the same field... exiting now\n!");
    exit (1);
  }

  mpz_add (out->value, in1->value, in2->value);
  out->field->modElt (out);

}

//================================
//================================

void subFieldElt (struct field_elt_t *in1, struct field_elt_t *in2, struct field_elt_t *out) {

  if (in1->field!=in2->field) {
    printf (" %p %p", in1->field,in2->field);
    printf ("error in subFieldElt: in1 and in2 are not elements of the same field... exiting now\n!");
    exit (1);
  }

  mpz_sub (out->value, in1->value, in2->value);
  out->field->modElt (out);
  
}

//================================
//================================

void negFieldElt (struct field_elt_t *in, struct field_elt_t *out) {

  if (in->field!=out->field) {
    printf (" %p %p", in->field, out->field);
    printf ("error in negFieldElt: in and out are not element of the same field... exiting now\n!");
    exit (1);
  }

  mpz_neg (out->value, in->value);
  out->field->modElt (out);

}

//================================
//================================
void powuiFieldElt (struct field_elt_t *in, struct field_elt_t *out, unsigned int exp) {

  if (in->field!=out->field) {
    printf (" %p %p", in->field, out->field);
    printf ("error in powuiFieldElt: in and out are not element of the same field... exiting now\n!");
    exit (1);
  }

  mpz_pow_ui (out->value, in->value, exp);
  out->field->modElt (out);
  
}
