#include <stdio.h>
#include <stdlib.h>

#include "poly.h"

void viewPoly (struct poly_t *this) {

  if (this->length!=0) {
    struct poly_node_t *node = this->head;
    printf (" ");
    
    while (node!=NULL) {  
      gmp_printf("%Zdx^%d", node->coefficient->value, node->exponent);
      node = node->next;
      if (node!=NULL) {
	printf(" + ");
      }
    }
    printf ("\n");
  }
  
  else {
    printf(" empty poly.\n");
  }
}

//================================
//================================

void initPoly (struct poly_t *this) {
  this->length = 0;

  this->view = viewPoly;
  this->free = freePoly;

  this->reset = resetPoly;
  
  this->rndPoly = randomPoly;
  
  this->insertNode = insertPolyNode;
  this->removeNode = removePolyNode;

  this->mulPoly = multiplyPoly;
  this->divPoly = dividePoly;
  this->addPoly = addPoly;
  this->subPoly = subtractPoly;
  
  this->cstMulPoly = cstMultiplyPoly;
  
  this->simplifyPoly = simplifyPoly;

  this->sort = sortPoly;
  
  this->copyPoly = copyPoly;

  this->hiExp = getHighestExp;

  this->eval = evaluatePoly;
}

//================================
//================================

void evaluatePoly (struct poly_t *poly, struct field_elt_t *x, struct field_elt_t *y ) {

  if (x->field!=y->field) {
    printf (" %p %p", x->field, y->field);
    printf ("error in evaluatePoly: x and y are not elements of the same field... exiting now\n!");
    exit (1);
  }
  
  if (poly->length!=0) {
    struct poly_node_t *node = poly->head;
    
    struct field_elt_t *tmp_elt = malloc (sizeof (struct field_elt_t) );
    tmp_elt->init = initFieldElt;
    tmp_elt->init (tmp_elt, x->field);

    mpz_set_ui (y->value, 0);

    while (node!=NULL) {
      tmp_elt->field->powuiElt (x, tmp_elt, node->exponent); // tmp_elt=x^exp [mod r]

      tmp_elt->field->mulElt (node->coefficient, tmp_elt, tmp_elt); // tmp_elt=coeff*tmp_elt [mod r]
      
      y->field->addElt (tmp_elt, y, y); //y+=tmp_elt [mod r]

      node = node->next;
    }

    tmp_elt->free (tmp_elt);
  }
  
  else {
    printf("empty polynomial... existing now!");
    exit(0);
  }
}


//================================
//================================

void randomPoly (struct poly_t *this, struct field_t *Fp, unsigned int highExp, gmp_randstate_t state) {
  
  struct field_elt_t *coeff = malloc (sizeof (struct field_elt_t) );
  coeff->init = initFieldElt;

  coeff->init (coeff, Fp);

  int i;
  for (i=highExp; i>=0; i--) {
    coeff->field->rndElt (coeff, state);
    this->insertNode (this, coeff->value, i, coeff->field);
    mpz_set_ui (coeff->value, 0);
  }
  
  coeff->free (coeff);
}

//================================
//================================

void insertPolyNode (struct poly_t *this, mpz_t coeff, unsigned int exp, struct field_t *Fp) {
  struct poly_node_t *node = malloc (sizeof (struct poly_node_t) );
  node->exponent = exp;

  node->coefficient = malloc (sizeof (struct field_elt_t) );
  node->coefficient->init = initFieldElt;
  node->coefficient->init (node->coefficient, Fp);

  mpz_set (node->coefficient->value, coeff);

  node->next = NULL;
  
  if (this->length==0) {
    this->head = this->tail = node;
  }
  else {
    this->tail->next = node;
    this->tail = node;
  }
  this->length++;
  
}

//================================
//================================

void freePoly (struct poly_t *this) {
  
  if (this!=NULL) {
  
    if (this->length>0) {
      
      struct poly_node_t *tmp;
      
      while (this->head != NULL) {
	tmp = this->head;
	tmp->coefficient->free (tmp->coefficient);

	this->head = this->head->next;
	
	free (tmp);
      }
      
    }
    free (this);
  }
}

//================================
//================================

void removePolyNode (struct poly_t *this, struct poly_node_t *node2remove) {
  struct poly_node_t *node = this->head;
  
  if (node==NULL) {
    printf("leaving delete_poly_ring_node");
  }
  
  if (node!=NULL && node==node2remove) {
    this->head=node->next;
    this->length-=1;
    node2remove->coefficient->free (node2remove->coefficient);
    free (node2remove);
    return;
  }
  
  while (node!=NULL && node->next!=node2remove)
    node=node->next;
  
  if (node2remove!=this->tail) {
    node->next=node2remove->next;
    this->length-=1;
    node2remove->coefficient->free (node2remove->coefficient);
    free (node2remove);
    return;
  }
  else {
    node->next=NULL;
    this->tail=node;
    this->length-=1;
    node2remove->coefficient->free (node2remove->coefficient);
    free (node2remove);
    return;
  } 
}

//================================
//================================

void multiplyPoly (struct poly_t *in1, struct poly_t *in2, struct poly_t *out) {

  if (in1==NULL||in2==NULL) {
    printf("error empty poly(s)... exiting now\n");
    exit(0);
  }

  struct poly_node_t *node_in1 = in1->head;
  struct poly_node_t *node_in2 = in2->head;

  struct field_elt_t *coeff = malloc (sizeof (struct field_elt_t));
  coeff->init = initFieldElt;
  coeff->init (coeff, node_in1->coefficient->field);
  
  unsigned int exp;
  
  while (node_in1!=NULL) {
    while (node_in2!=NULL) {
      coeff->field->mulElt (node_in1->coefficient, node_in2->coefficient, coeff);
      exp = node_in1->exponent + node_in2->exponent;
      out->insertNode (out, coeff->value, exp, coeff->field);
      node_in2 = node_in2->next;
    }
    node_in1 = node_in1->next;
    node_in2 = in2->head;
  }
  
  out->simplifyPoly (out);

  coeff->free (coeff);
  
}

//================================
//================================

void dividePoly (struct poly_t *in1, struct poly_t *in2, struct poly_t *out) {
  
  if (in1==NULL||in2==NULL) {
    printf("error in divividePoly, empty poly(s)... exiting now\n");
    exit(0);
  }

  struct poly_t *dividend = malloc (sizeof (struct poly_t));
  dividend->init = initPoly;
  dividend->init (dividend);
  dividend->copyPoly (in1, dividend);
  
  struct poly_t *quotient = malloc (sizeof (struct poly_t));
  quotient->init = initPoly;
  quotient->init (quotient);

  struct poly_t *remainder = malloc (sizeof (struct poly_t));
  remainder->init = initPoly;
  remainder->init (remainder);

  struct poly_t *tmp = malloc (sizeof (struct poly_t));
  tmp->init = initPoly;
  tmp->init (tmp);

  struct poly_t *dividend_tmp = malloc (sizeof (struct poly_t));
  dividend_tmp->init = initPoly;
  dividend_tmp->init (dividend_tmp);
  
  struct poly_node_t *node_dividend = dividend->head;
  struct poly_node_t *node_divisor  = in2->head;

  struct field_elt_t *coeff = malloc (sizeof (struct field_elt_t));
  coeff->init = initFieldElt;
  coeff->init (coeff, dividend->head->coefficient->field);
   
  unsigned int exp;

  if ((int)(getHighestExp (dividend)) - (int)(getHighestExp(in2))<0) {
    out->insertNode (out, coeff->value, 0, coeff->field);
  }

  else {
    while (node_dividend->exponent >= getHighestExp (in2) && dividend->length > 0) {
      coeff->field->divElt (node_dividend->coefficient, node_divisor->coefficient, coeff);

      exp = node_dividend->exponent - node_divisor->exponent;

      quotient->insertNode (quotient, coeff->value, exp, coeff->field);

      out->insertNode (out, coeff->value, exp, coeff->field);

      tmp->mulPoly (quotient, in2, tmp);
 
      dividend->subPoly (dividend, tmp, dividend_tmp);
      
      dividend->reset (dividend);

      dividend_tmp->removeNode (dividend_tmp, dividend_tmp->head);

      dividend->copyPoly (dividend_tmp, dividend);

      quotient->removeNode (quotient, quotient->head);

      tmp->reset (tmp);

      dividend_tmp->reset (dividend_tmp);
      
      if (dividend->length>0) {
	node_dividend = dividend->head;
      }
    }
    remainder->copyPoly (dividend, remainder);
    //printf(" remainder:");
    //remainder->view (remainder);
  }

  dividend->free (dividend);
  dividend_tmp->free (dividend_tmp);
  quotient->free (quotient);
  remainder->free (remainder);
  tmp->free (tmp);
  coeff->free (coeff);
}

//================================
//================================

void addPoly (struct poly_t *in1, struct poly_t *in2, struct poly_t *out) {

  if (in1==NULL||in2==NULL) {
    printf("error in addPoly, null poly(s)... exiting now\n");
    exit (0);
  }

  struct poly_t *tmp = malloc (sizeof (struct poly_t));
  tmp->init = initPoly;
  tmp->init (tmp);
  
  struct poly_node_t *node_in1 = in1->head;
  struct poly_node_t *node_in2 = in2->head;
  
  struct field_elt_t *coeff = malloc (sizeof (struct field_elt_t));
  coeff->init = initFieldElt;
  coeff->init (coeff, in1->head->coefficient->field);

  while (node_in1!=NULL && node_in2!=NULL) {
    if (node_in1->exponent == node_in2->exponent) {
      coeff->field->addElt (node_in1->coefficient, node_in2->coefficient, coeff);
      tmp->insertNode (tmp, coeff->value, node_in1->exponent, coeff->field);
      node_in1 = node_in1->next;
      node_in2 = node_in2->next;
      if (node_in2==NULL) {
	while (node_in1!=NULL) {
	  tmp->insertNode (tmp, node_in1->coefficient->value, node_in1->exponent, node_in1->coefficient->field);
	  node_in1 = node_in1->next;
	}
      }
      if (node_in1==NULL) {
	while (node_in2!=NULL) {
	  tmp->insertNode (tmp, node_in2->coefficient->value, node_in2->exponent, node_in2->coefficient->field);
	  node_in2 = node_in2->next;
	}
      }
    } // end-if
    
    else if (node_in1->exponent > node_in2->exponent) {
      tmp->insertNode (tmp, node_in1->coefficient->value, node_in1->exponent, node_in1->coefficient->field);
      node_in1 = node_in1->next;
      if (node_in1==NULL) {
	while (node_in2!=NULL) {
	  tmp->insertNode (tmp, node_in2->coefficient->value, node_in2->exponent, node_in2->coefficient->field);
	  node_in2 = node_in2->next;
	}
      }
    } // end-else if 
    
    else if (node_in1->exponent < node_in2->exponent) {
      tmp->insertNode (tmp, node_in2->coefficient->value, node_in2->exponent, node_in2->coefficient->field);
      node_in2 = node_in2->next;
      if (node_in2==NULL) {
	while( node_in1 != NULL ) {
	  tmp->insertNode (tmp, node_in1->coefficient->value, node_in1->exponent, node_in1->coefficient->field);
	  node_in1 = node_in1->next;
	}
      }
    } // end-else if

    else { printf("error in add_poly... exiting now!\n"); exit(0); }
    
  }

  out->copyPoly (tmp, out);
  out->simplifyPoly (out);
  
  coeff->free (coeff);
  tmp->free (tmp);
  
}
  
//================================
//================================

void subtractPoly (struct poly_t *in1, struct poly_t *in2, struct poly_t *out) {

  if (in1==NULL||in2==NULL) {
    printf("error in subtractPoly, null poly(s)... exiting now\n");
    exit (0);
  }
  
  struct poly_t *tmp = malloc (sizeof (struct poly_t));
  tmp->init = initPoly;
  tmp->init (tmp);
  
  struct poly_node_t *node_in1 = in1->head;
  struct poly_node_t *node_in2 = in2->head;
  
  struct field_elt_t *coeff = malloc (sizeof (struct field_elt_t));
  coeff->init = initFieldElt;
  coeff->init (coeff, in1->head->coefficient->field);
  
  while (node_in1!=NULL && node_in2!=NULL) {
    if (node_in1->exponent == node_in2->exponent) {
      coeff->field->subElt (node_in1->coefficient, node_in2->coefficient, coeff);
      tmp->insertNode (tmp, coeff->value, node_in1->exponent, coeff->field);
      node_in1 = node_in1->next;
      node_in2 = node_in2->next;
      if (node_in1==NULL) {
	while (node_in2!=NULL) {
	  coeff->field->negElt (node_in2->coefficient, coeff); 
	  tmp->insertNode (tmp, coeff->value, node_in2->exponent, coeff->field);
	  node_in2 = node_in2->next;
	}
      }
      if (node_in2==NULL) {
	while (node_in1!=NULL) {
	  tmp->insertNode (tmp, node_in1->coefficient->value, node_in1->exponent, node_in1->coefficient->field);
	  node_in1=node_in1->next;
	}
      }
    } // end-if
    
    else if (node_in1->exponent > node_in2->exponent) {
      tmp->insertNode (tmp, node_in1->coefficient->value, node_in1->exponent, node_in1->coefficient->field);
      node_in1 = node_in1->next;
      if (node_in1==NULL) {
	while (node_in2!=NULL) {
	  coeff->field->negElt (coeff, node_in2->coefficient);
	  tmp->insertNode (tmp, coeff->value, node_in2->exponent, coeff->field);
	  node_in2 = node_in2->next;
	}
      }
    } // end-else if
    
    else if (node_in1->exponent < node_in2->exponent) {
      coeff->field->negElt (node_in2->coefficient, coeff);
      tmp->insertNode (tmp, coeff->value, node_in2->exponent, node_in2->coefficient->field);
      node_in2 = node_in2->next;
      if (node_in2==NULL) {
	while (node_in1!=NULL) {
	  tmp->insertNode (tmp, node_in1->coefficient->value, node_in1->exponent, node_in1->coefficient->field);
	  node_in1 = node_in1->next;
	}
      }
    } // end-else if 
    
    else {
      printf("error in sub_poly... exiting now!\n");
      exit(0);
    }
    
  }

  out->copyPoly (tmp, out);
  out->simplifyPoly (out);
 
  coeff->free (coeff);
  tmp->free (tmp);
  
}

//================================
//================================

void simplifyPoly (struct poly_t *this) {

  struct poly_node_t *node = this->head;
  
  struct poly_node_t *next_node = node->next;

  struct poly_node_t *tmp;

  while (node!=NULL) {
    while (next_node!=NULL) {
      if (node->exponent==next_node->exponent) {
	node->coefficient->field->addElt (node->coefficient, next_node->coefficient, node->coefficient);
	tmp = next_node->next;
	this->removeNode (this, next_node);
	next_node = tmp;
      }
      else {
	next_node = next_node->next;
      }
    }
    node = node->next;
    if (node!= NULL) {
      next_node = node->next;
    }
  }
  
  if (this->head->exponent < this->tail->exponent) {
    this->sort (this);
  }
  
}

//================================
//================================

void copyPoly (struct poly_t *in, struct poly_t *out) {
  
  if (in==NULL||out==NULL) {
    printf("error in copyPoly, empty poly(s)... exiting now\n");
    exit(0);
  }

  struct poly_node_t *node = in->head;
  
  while (node!=NULL) {
    out->insertNode
      (out, node->coefficient->value, node->exponent, node->coefficient->field);
    node = node->next;
  }
}

//================================
//================================

unsigned int getHighestExp (struct poly_t *this) {
  
  if (this==NULL) {
    printf("error in getHighestExp, empty poly(s)... exiting now\n");
    exit(0);
  }

  unsigned int exp = 0;
  struct poly_node_t *node = this->head;

  while (node!=NULL) {
    if (node->exponent > exp) {
      exp = node->exponent;
    }
    node = node->next;
  }
  return exp;
}

//================================
//================================

void resetPoly (struct poly_t *this) {
  
  if (this!=NULL) {
  
    if (this->length>0) {
      this->length=0;
      struct poly_node_t *tmp;
      
      while (this->head != NULL) {
	tmp = this->head;
	tmp->coefficient->free (tmp->coefficient);

	this->head = this->head->next;
	
	free (tmp);
      }
    }
  }
}

//================================
//================================

void cstMultiplyPoly (struct poly_t *in, struct field_elt_t *elt, struct poly_t *out) {

  if (in==NULL||elt==NULL) {
    printf ("error in cstMultiplypoly, empty poly or field element... exiting now\n");
    exit (0);
  }

  struct poly_node_t *node_in = in->head;

  struct field_elt_t *coeff = malloc (sizeof (struct field_elt_t));
  coeff->init = initFieldElt;
  coeff->init (coeff, node_in->coefficient->field);
  
  while (node_in!=NULL) {
    coeff->field->mulElt (node_in->coefficient, elt, coeff);
    out->insertNode (out, coeff->value, node_in->exponent, coeff->field);
    node_in = node_in->next;
  }

  coeff->free (coeff);
}

//================================
//================================

void sortPoly (struct poly_t *this) {

  struct poly_node_t *node;
  
  struct poly_t *tmp_poly = malloc (sizeof (struct poly_t) );
  tmp_poly->init = initPoly;
  tmp_poly->init (tmp_poly);
  tmp_poly->copyPoly (this, tmp_poly);
 
  this->reset (this);
  
  unsigned int hiExp;
  int flag;
  
  while (tmp_poly->length>0) {

    hiExp = tmp_poly->hiExp (tmp_poly);
    node = tmp_poly->head;

    flag = 0;
    
    while (node->next!=NULL || flag==0) {
      if (node->exponent==hiExp) {
	this->insertNode (this, node->coefficient->value, node->exponent, node->coefficient->field);
	tmp_poly->removeNode (tmp_poly, node);
	flag = 1;
      }
      else {
	node = node->next;
      }
    }
    
    if (tmp_poly->length==1) {
      this->insertNode (this, node->coefficient->value, node->exponent, node->coefficient->field);
      tmp_poly->removeNode (tmp_poly, node);
    }
  }

  tmp_poly->free (tmp_poly);
}

