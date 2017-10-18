#include <stdio.h>
#include <string.h>
#include <gmp.h>

#include "circuit.h"

#define bufSize 1024
#define strSize 20

void init_circuit( struct circuit_t *c, struct field_t *Fp, FILE *fptr, gmp_randstate_t rnd_state ) {
  char buf[bufSize];
  char str[strSize];
  int nInput_wire = 0;
  int nMul_gate   = 0;
  int nAdd_gate   = 0;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    //printf("%s\n", buf);
    sscanf(buf,"%s",str);
    if( strcmp(str,"input")==0 ) {nInput_wire+=1;}
    else if( strcmp(str,"mul")==0 ) {nMul_gate+=1;}
    else if( strcmp(str,"add")==0 ) {nAdd_gate+=1;}
  }

  c->nWire = nInput_wire + nMul_gate + nAdd_gate;
  c->nMul_gate = nMul_gate;
  c->nAdd_gate = nAdd_gate;
  c->nGate = c->nMul_gate + c->nAdd_gate;
  c->index = malloc( nMul_gate * sizeof( unsigned int ));
  
  c->wire = malloc( (c->nWire) * sizeof( struct wire_t ));
  //c->mul_gate = malloc( (c->nMul_gate) * sizeof( struct gate_t ));
  //c->add_gate = malloc( (c->nAdd_gate) * sizeof( struct gate_t ));
  c->gate = malloc( (c->nGate) * sizeof( struct gate_t ));
  
  printf ("\nnb wires: %d\n", c->nWire);
  printf ("nb gates: %d\n", c->nGate); 
  printf ("nb mult. gates: %d\n", c->nMul_gate);
  printf ("nb add.gates: %d\n", c->nAdd_gate);
  
  rewind( fptr );
  int i_wire = 0;
  int i_gate = 0;
  int i_index = 0;
  //int i_mulgate = 0;
  //int i_addgate = 0;
  int i, j, k, l, m;
  char foo[strSize];
  int dummy;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    sscanf(buf,"%s",str);
    if( strcmp(str,"input")==0 ) {
      sscanf(buf,"%*s %d", &i);
      c->wire[i_wire].id = i;
      c->wire[i_wire].type = "input";
      mpz_init( c->wire[i_wire].value );
      i_wire+=1;
    }
    
    else if( strcmp(str,"mul")==0 ) {
      sscanf(buf,"%*s %*s %d %d %d %s %d %d", &i, &j, &k, foo, &l, &m);
      
      if( strcmp(foo,"out")==0 ) {
	c->wire[i_wire].id = m;
	c->wire[i_wire].type = "intermediate";
	mpz_init( c->wire[i_wire].value );
	i_wire+=1;
      }
      /*
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
      */
      c->gate[i_gate].type = "mul";
      c->gate[i_gate].field = Fp;
      rnd_field_elt( c->gate[i_gate].label, c->gate[i_gate].field, rnd_state);
      c->gate[i_gate].ninput = i;
      c->gate[i_gate].noutput = l;
      for (dummy=0; dummy<c->nWire; dummy++) {
      	if      (c->wire[dummy].id==j) {c->gate[i_gate].linput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==k) {c->gate[i_gate].rinput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==m) {c->gate[i_gate].output=&(c->wire[dummy]);}
      }
      c->index[i_index]=i_gate;
      i_index+=1;
      i_gate+=1;
    }

    else if( strcmp(str,"add")==0 ) {
      sscanf(buf,"%*s %*s %d %d %d %s %d %d", &i, &j, &k, foo, &l, &m);
      
      if( strcmp(foo,"out")==0 ) {
	c->wire[i_wire].id = m;
	c->wire[i_wire].type = "intermediate";
	mpz_init( c->wire[i_wire].value );
	i_wire+=1;
      }
      /*
      c->add_gate[i_addgate].type = "add";
      c->add_gate[i_addgate].field = Fp;
      mpz_init( c->add_gate[i_addgate].label );
      mpz_set_si( c->add_gate[i_addgate].label, -99 );
      c->add_gate[i_addgate].ninput = i;
      c->add_gate[i_addgate].noutput = l;
      for (dummy=0; dummy<c->nWire; dummy++) {
      	if      (c->wire[dummy].id==j) {c->add_gate[i_addgate].linput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==k) {c->add_gate[i_addgate].rinput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==m) {c->add_gate[i_addgate].output=&(c->wire[dummy]);}
      }
      i_addgate+=1;
      */

      c->gate[i_gate].type = "add";
      c->gate[i_gate].field = Fp;
      mpz_init( c->gate[i_gate].label );
      mpz_set_si( c->gate[i_gate].label, -99 );
      c->gate[i_gate].ninput = i;
      c->gate[i_gate].noutput = l;
      for (dummy=0; dummy<c->nWire; dummy++) {
      	if      (c->wire[dummy].id==j) {c->gate[i_gate].linput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==k) {c->gate[i_gate].rinput=&(c->wire[dummy]);}
	else if (c->wire[dummy].id==m) {c->gate[i_gate].output=&(c->wire[dummy]);}
      }
      i_gate+=1;
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

//==============================================
//==============================================

void display_circuit( struct circuit_t *c ) {
  int i;

  printf("\nwire ID | type(input, intermed., output) | value\n");
  for (i=0; i<c->nWire; i++) {
    gmp_printf("wire %d %s %Zd\n", c->wire[i].id, c->wire[i].type, c->wire[i].value );
  }
  
  /*  
  for (i=0; i<c->nMul_gate; i++) {
    //gmp_printf("%s %Zd %Zd %d %d %d %d %d\n",
    gmp_printf("%s %Zd %d %d %d %d %d\n",
	       c->mul_gate[i].type,
	       //c->mul_gate[i].field->p,
	       c->mul_gate[i].label, 
	       c->mul_gate[i].ninput,
	       c->mul_gate[i].linput->id,
	       c->mul_gate[i].rinput->id,
	       c->mul_gate[i].noutput,
	       c->mul_gate[i].output->id);
  }
  
  for (i=0; i<c->nAdd_gate; i++) {
    //gmp_printf("%s %Zd %Zd %d %d %d %d %d\n",
    gmp_printf("%s %Zd %d %d %d %d %d\n",
	       c->add_gate[i].type,
	       //c->add_gate[i].field->p,
	       c->add_gate[i].label, 
	       c->add_gate[i].ninput,
	       c->add_gate[i].linput->id,
	       c->add_gate[i].rinput->id,
	       c->add_gate[i].noutput,
	       c->add_gate[i].output->id);  
  }
  */
  printf("\ngate type(add, mul,...) | label | nInput | left wire ID | right wire ID | nOutput | output wire ID\n");
  for (i=0; i<c->nGate; i++) {
    //gmp_printf("%s %Zd %Zd %d %d %d %d %d\n",
    gmp_printf("%s %Zd %d %d %d %d %d\n",
	       c->gate[i].type,
	       //c->add_gate[i].field->p,
	       c->gate[i].label, 
	       c->gate[i].ninput,
	       c->gate[i].linput->id,
	       c->gate[i].rinput->id,
	       c->gate[i].noutput,
	       c->gate[i].output->id);  
  }
}

//==============================================
//==============================================

void display_mul_gate( struct circuit_t *c ) {
  int i;
  for (i=0; i<c->nMul_gate; i++) {
    gmp_printf("%s %Zd %d %d %d %d %d\n",
	       c->gate[c->index[i]].type,
	       //c->add_gate[i].field->p,
	       c->gate[c->index[i]].label, 
	       c->gate[c->index[i]].ninput,
	       c->gate[c->index[i]].linput->id,
	       c->gate[c->index[i]].rinput->id,
	       c->gate[c->index[i]].noutput,
	       c->gate[c->index[i]].output->id);
  }
}

//==============================================
//==============================================

void set_random_input_values( struct circuit_t *c, struct field_t *Fp, gmp_randstate_t rnd_state ) {
  int i;
  for (i=0; i<c->nWire; i++) {
    if ( strcmp(c->wire[i].type,"input")==0) {
      //mpz_init( c->wire[i].value );
      rnd_field_elt( c->wire[i].value, Fp, rnd_state);
    }
  }
}
