#include <stdio.h>
#include <string.h>
#include <gmp.h>

#include "circuit.h"

#define bufSize 1024
#define strSize 20

void init_circuit( struct circuit_t *c, struct field_t *Fp, FILE *fptr, gmp_randstate_t rnd_state ) {
  //void init_circuit( struct circuit_t *c, FILE *fptr, gmp_randstate_t state ) {
  char buf[bufSize];
  char str[strSize];
  int nInput_wire = 0;
  int nMul_gate = 0;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    //printf("%s\n", buf);
    sscanf(buf,"%s",str);
    if( strcmp(str,"input")==0 ) {nInput_wire+=1;}
    else if( strcmp(str,"mul")==0 ) {nMul_gate+=1;}
  }

  c->nWire = nInput_wire+nMul_gate;
  c->nMul_gate = nMul_gate;
    
  c->wire = malloc( (c->nWire) * sizeof( struct wire_t ));
  c->mul_gate = malloc( (c->nMul_gate) * sizeof( struct gate_t ));

  printf ("\nnb wires: %d\n", c->nWire);
  printf ("nb mult_gates: %d\n", c->nMul_gate);
  
  rewind( fptr );
  int i_wire  = 0;
  int i_mulgate = 0;
  int i, j, k, l, m;
  char foo[strSize];
  int dummy;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    sscanf(buf,"%s",str);
    if( strcmp(str,"input")==0 ) {
      sscanf(buf,"%*s %d", &i);
      c->wire[i_wire].id = i;
      c->wire[i_wire].type = "input";
      i_wire+=1;
    }
    
    else if( strcmp(str,"mul")==0 ) {
      sscanf(buf,"%*s %*s %d %d %d %s %d %d", &i, &j, &k, foo, &l, &m);
      
      if( strcmp(foo,"out")==0 ) {
	c->wire[i_wire].id = m;
	c->wire[i_wire].type = "intermediate";
	i_wire+=1;
      }
      c->mul_gate[i_mulgate].type = "mul";
      c->mul_gate[i_mulgate].field = Fp;
      rnd_field_elt( c->mul_gate[i_mulgate].label, c->mul_gate[i_mulgate].field, rnd_state);
      c->mul_gate[i_mulgate].ninput = i;
      c->mul_gate[i_mulgate].noutput = l;
      for (dummy=0; dummy<c->nWire; dummy++) {
      	if      (c->wire[dummy].id==j) {c->mul_gate[i_mulgate].linput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==k) {c->mul_gate[i_mulgate].rinput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==m) {c->mul_gate[i_mulgate].output=&(c->wire[dummy]);}
      }
      i_mulgate+=1;
    }

    else if( strcmp(str,"output")==0 ) {
      sscanf(buf,"%*s %d", &i);
      k=-1;
      for (j=0; j<c->nWire; j++) {
	if (c->wire[j].id==i) {
	  c->wire[j].type = "out";
	  k=1;
	}
      }
      if (k<0) {printf("error: couldn't find output wire %d... exiting now!\n",i); exit(1);}
    }

    else {printf("error: %s not defined... exiting now!\n",str); exit(1);}
  }
}
