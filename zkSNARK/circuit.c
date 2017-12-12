#include "circuit.h"

#define bufSize 1024
#define strSize 20

void initCircuit (struct circuit_t *this, struct field_t *Fp, FILE *fptr, gmp_randstate_t state) {

  this->field = Fp;
  
  this->view = viewCircuit;
  this->free = freeCircuit;
  this->eval = evalCircuit;
  
  this->newWire = newWire;
  this->newGate = newGate;

  this->rndInputs = setRNDInputs;
    
  char buf[bufSize];
  char str[strSize];
  unsigned int nInput_wire = 0;
  unsigned int nMul_gate   = 0;
  unsigned int nAdd_gate   = 0;
  
  while (fgets (buf, sizeof(buf), fptr)!=NULL) {
    sscanf(buf,"%s",str);
    if (strcmp(str,"input")==0) {
      nInput_wire+=1;
    }
    else if (strcmp(str,"mul")==0) {
      nMul_gate+=1;
    }
    else if (strcmp(str,"add")==0) {
      nAdd_gate+=1;
    }
  }

  this->nWire = nInput_wire + nMul_gate + nAdd_gate;
  this->nIOWire = nInput_wire;
  this->nMul_gate = nMul_gate;
  this->nAdd_gate = nAdd_gate;
  this->nGate = this->nMul_gate + this->nAdd_gate;
  this->index = malloc (this->nMul_gate * sizeof (unsigned int));

  // allocates an array of pointers to wire_t
  this->wire = malloc (this->nWire * sizeof (struct wire_t*)); 

  // allocates an array of pointers to gate_t
  this->gate = malloc (this->nGate * sizeof (struct gate_t*));

  rewind( fptr );

  int i, j, k;

  char foo [strSize];
  
  unsigned int *igate  = malloc (sizeof (unsigned int));
  unsigned int *iwire  = malloc (sizeof (unsigned int));
  unsigned int *iindex = malloc (sizeof (unsigned int));
  
  *(igate)  = 0;
  *(iwire)  = 0;
  *(iindex) = 0;

  while (fgets(buf, sizeof(buf), fptr)!=NULL ) {
    sscanf (buf,"%s",str);

    if (strcmp(str,"input")==0) {
      this->newWire (this, iwire, buf, "input");
    } // end if 'input'

    else if (strcmp(str,"mul")==0) {
      sscanf (buf,"%*s %*s %*d %*d %*d %s", foo);
      if (strcmp(foo,"out")==0) {
	this->newWire (this, iwire, buf, "intermediate");
      }
      this->newGate (this, igate, iindex, *(iwire), buf, "mul", state);
    } // end else-if 'mul'

    else if (strcmp(str,"add")==0) {
      sscanf(buf,"%*s %*s %*d %*d %*d %s", foo);
      if (strcmp(foo,"out")==0) {
	this->newWire (this, iwire, buf, "intermediate");
      }
      this->newGate (this, igate, iindex, *(iwire), buf, "add", state);
    } // end else-if 'add'

    else if( strcmp(str,"output")==0 ) {
      sscanf(buf,"%*s %d", &i);
      this->nIOWire += 1;
      k=-1;
      for (j=0; j<*(iwire); j++) {
	if (this->wire[j]->id==i) {
	  strcpy (this->wire[j]->type, "output");
	  k=1;
	}
      }
      if (k<0) {
	printf("error in initCircuit: couldn't find output wire %d,... exiting now!\n",i);
	exit (1);
      }
    } // end else-if 'output'
    
    else {
      printf ("error in initCircuit: wire type %s not defined,... exiting now!\n",str);
      exit (1);
    } // else
    
  } // end while

  this->rndInputs (this, state);
  
  this->view (this);

  free (igate);
  free (iwire);
  free (iindex);
  
  fclose (fptr);
}

//================================
//================================

void newGate (struct circuit_t *this, unsigned int *igate, unsigned int *iindex, unsigned int iwire, char buf[bufSize], char str[strSize], gmp_randstate_t state) {
  this->gate[*(igate)] = malloc (sizeof (struct gate_t));
  this->gate[*(igate)]->init = initGate;
  this->gate[*(igate)]->init (this->gate[*(igate)], str);
  
  int i, j, k, l, m;
  sscanf(buf,"%*s %*s %d %d %d %*s %d %d", &i, &j, &k, &l, &m);

  strcpy (this->gate[*(igate)]->type, str);

  if (strcmp(str,"mul")==0) {
    this->gate[*(igate)]->root->init (this->gate[*(igate)]->root, this->field);
    this->gate[*(igate)]->root->field->rndElt (this->gate[*(igate)]->root, state);
    this->index[*(iindex)]=*(igate);
    *(iindex)+=1;
  }

  this->gate[*(igate)]->ninput  = i;
  this->gate[*(igate)]->noutput = l;

  int dummy;
  for (dummy=0; dummy<iwire; dummy++) {
    
    if (this->wire[dummy]->id==j) {
      this->gate[*(igate)]->linput=this->wire[dummy];
      if (strcmp(str,"mul")==0) {
	this->wire[dummy]->ctrb (this->wire[dummy]->left_ctrb, this->gate[*(igate)]);
      }
    }
    
    if (this->wire[dummy]->id==k) {
      this->gate[*(igate)]->rinput=this->wire[dummy];
      if (strcmp(str,"mul")==0) {
	this->wire[dummy]->ctrb (this->wire[dummy]->right_ctrb, this->gate[*(igate)]);
      }
    }
    
    if (this->wire[dummy]->id==m) {
      this->gate[*(igate)]->output=this->wire[dummy];
      if (strcmp(str,"mul")==0) {
	this->wire[dummy]->ctrb (this->wire[dummy]->output_ctrb, this->gate[*(igate)]);
      }
    }
  }
  
  *(igate)+=1;
  
}

//================================
//================================

void newWire (struct circuit_t *this, unsigned int *iwire, char buf[bufSize], char str[strSize]) {
  this->wire[*(iwire)] = malloc (sizeof (struct wire_t));
  this->wire[*(iwire)]->init = initWire;
  this->wire[*(iwire)]->init (this->wire[*(iwire)]);

  int i;
  if (strcmp(str,"input")==0) {
    sscanf(buf,"%*s %d", &i);
  }
  else if (strcmp(str,"intermediate")==0) {
    sscanf (buf,"%*s %*s %*d %*d %*d %*s %*d %d", &i);
  }
  else {
    printf ("error: str %s not defined... exiting now!\n",str);
    exit (1);
  }

  this->wire[*(iwire)]->id = i;

  strcpy (this->wire[*(iwire)]->type, str);
  
  this->wire[*(iwire)]->value->init (this->wire[*(iwire)]->value, this->field);

  *(iwire)+=1;

}

//================================
//================================

void viewCircuit (struct circuit_t *this) {
  int i;
  printf ("\n");
  printf (" **** Circuit ****\n");
  printf (" nb wires: %d\n", this->nWire);
  printf (" nb nIOWire: %d\n", this->nIOWire);
  printf (" nb gates: %d\n", this->nGate); 
  printf (" nb mult. gates: %d\n", this->nMul_gate);
  printf (" nb add.gates: %d\n", this->nAdd_gate);
 
  printf ("\n");
  //printf (" wire ID | type(input, intermed., output) | value(dec, hex) | gate mul. ID(left contr.) | gate mul. ID(right contr.) | gate mul. ID(output contr.)\n");
  printf (" wire ID | type(input, intermed., output) | value(dec, hex)\n");
  for (i=0; i<this->nWire; i++) {
    this->wire[i]->view (this->wire[i]);
  }
  
  printf ("\n");
  printf(" gate type(add, mul,...) | root(dec, hex) | nInput | left wire ID | right wire ID | nOutput | output wire ID\n");
  for (i=0; i<this->nGate; i++) {
    this->gate[i]->view (this->gate[i]);
  }
}

//================================
//================================

void freeCircuit (struct circuit_t *this) {
  int i;
  
  if (this->index!=NULL) {
    free (this->index);
  }

  for (i=0; i<this->nWire; i++) {
    this->wire[i]->free (this->wire[i]);
  }
  free (this->wire);
  
  for (i=0; i<this->nGate; i++) {
    this->gate[i]->free (this->gate[i]);
  }
  free (this->gate);
  
  free (this);
  
}

//================================
//================================

void setRNDInputs (struct circuit_t *this, gmp_randstate_t state) {
  int i;
  for (i=0; i<this->nWire; i++) {
    
    if (strcmp(this->wire[i]->type,"input")==0) {
      
      this->wire[i]->value->field->rndElt
	(this->wire[i]->value, state);
      
    }
  }
}

//================================
//================================

void evalCircuit (struct circuit_t *this) {
  int i;
   for (i=0; i<this->nGate; i++) {
    this->gate[i]->eval (this->gate[i]);
  }
}
