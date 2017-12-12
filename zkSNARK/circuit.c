#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>

#include "circuit.h"

#define bufSize 1024
#define strSize 20

//==============================================
//==============================================

void eval_circuit( struct circuit_t *c ) {
  int i;
  printf ("\n**** Circuit Evaluation ****\n");
  for (i=0; i<c->nGate; i++) {
    eval_gate( c->gate[i] );
  }
}

//==============================================
//==============================================

void init_circuit( struct circuit_t *c, struct field_t *Fp, FILE *fptr, gmp_randstate_t rnd_state ) {
  char buf[bufSize];
  char str[strSize];
  unsigned int nInput_wire = 0;
  unsigned int nMul_gate   = 0;
  unsigned int nAdd_gate   = 0;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    //printf("%s\n", buf);
    sscanf(buf,"%s",str);
    if( strcmp(str,"input")==0 )    {nInput_wire+=1;}
    else if( strcmp(str,"mul")==0 ) {nMul_gate+=1;}
    else if( strcmp(str,"add")==0 ) {nAdd_gate+=1;}
  }
  
  c->nWire = nInput_wire + nMul_gate + nAdd_gate;
  c->nIOWire = nInput_wire;
  c->nMul_gate = nMul_gate;
  c->nAdd_gate = nAdd_gate;
  c->nGate = c->nMul_gate + c->nAdd_gate;
  c->index = malloc( c->nMul_gate * sizeof( unsigned int ));
 
  c->wire = malloc( (c->nWire) * sizeof( struct wire_t* )); //allocates an array of pointers to wire_t
  c->gate = malloc( (c->nGate) * sizeof( struct gate_t* )); //allocates an array of pointers to gate_t

  //c->wire = (struct wire_t*) calloc( (c->nWire), sizeof( struct wire_t* )); //allocates an array of pointers to wire_t
  //c->gate = (struct gate_t*) calloc( (c->nGate), sizeof( struct gate_t* )); //allocates an array of pointers to gate_t

  printf ("\n**** Circuit Initialization ****\n\n");
  printf ("nb wires: %d\n", c->nWire);
  printf ("nb gates: %d\n", c->nGate); 
  printf ("nb mult. gates: %d\n", c->nMul_gate);
  printf ("nb add.gates: %d\n", c->nAdd_gate);

  rewind( fptr );

  int i, j, k;

  char foo[strSize];
  
  unsigned int *igate = malloc( sizeof( unsigned int ));
  unsigned int *iwire = malloc( sizeof( unsigned int ));
  unsigned int *iindex = malloc( sizeof( unsigned int ));
  
  *(igate)  = 0;
  *(iwire)  = 0;
  *(iindex) = 0;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    sscanf(buf,"%s",str);
   
    if( strcmp(str,"input")==0 ) {
      new_wire(iwire, buf, "input", c, Fp);
      //print_wire(c->wire[*(iwire)-1]);
    }
   
    else if( strcmp(str,"mul")==0 ) {
      sscanf(buf,"%*s %*s %*d %*d %*d %s", foo);
      if( strcmp(foo,"out")==0 ) {
	new_wire(iwire, buf, "intermediate", c, Fp);
	//print_wire(c->wire[*(iwire)-1]);
      }
      new_gate(igate, iindex, *(iwire), buf, "mul", c, Fp, rnd_state);
      //print_gate(c->gate[*(igate)-1]);
      //exit(0);
    }//end "mul"
    
    else if( strcmp(str,"add")==0 ) {
      sscanf(buf,"%*s %*s %*d %*d %*d %s", foo);
      if( strcmp(foo,"out")==0 ) {
	new_wire(iwire, buf, "intermediate", c, Fp);
	//print_wire(c->wire[*(iwire)-1]);
      }
      new_gate(igate, iindex, *(iwire), buf, "add", c, Fp, rnd_state);
    }//end "add"
    
    else if( strcmp(str,"output")==0 ) {
      sscanf(buf,"%*s %d", &i);
      c->nIOWire += 1;
      k=-1;
      for (j=0; j<*(iwire); j++) {
	if (c->wire[j]->id==i) {
	  //c->wire[j]->type = "output";
	  strcpy (c->wire[j]->type, "output");
	  k=1;
	}
      }
      if (k<0) {printf("error: couldn't find output wire %d... exiting now!\n",i); exit(1);}
    }//end "output"
   
    else {printf("error: %s not defined... exiting now!\n",str); exit(1);}
  }
  
  printf ("nb IO wires: %d\n", c->nIOWire);

  free (igate);
  free (iwire);
  free (iindex);
  
  //exit(0);
  //for (i=0; i<c->nWire; i++) {print_wire_contribution_to_mul_gate(c->wire[i]);}
  //display_circuit(c);
  //exit(0);
}

//==============================================
//==============================================

void display_circuit( struct circuit_t *c ) {
  int i;

  printf("\nwire ID | type(input, intermed., output) | value(dec, hex) | gate mul. ID(left contr.) | gate mul. ID(right contr.) | gate mul. ID(output contr.)\n");

  for (i=0; i<c->nWire; i++) {
    gmp_printf("\n %d %s %Zd %#Zx", c->wire[i]->id, c->wire[i]->type, c->wire[i]->value, c->wire[i]->value );
    print_wire_contribution_to_mul_gate(c->wire[i]);
  }

  printf("\n\ngate type(add, mul,...) | root(dec, hex) | nInput | left wire ID | right wire ID | nOutput | output wire ID\n\n");
  for (i=0; i<c->nGate; i++) {
    //gmp_printf("%s %Zd %Zd %d %d %d %d %d\n",
    if ( strcmp(c->gate[i]->type,"mul")==0) {
      display_mul_gate( c->gate[i] );
    }
    else {
      gmp_printf(" %s NULL %d %d %d %d %d\n",
		 c->gate[i]->type,
		 c->gate[i]->ninput,
		 c->gate[i]->linput->id,
		 c->gate[i]->rinput->id,
		 c->gate[i]->noutput,
		 c->gate[i]->output->id);  
    }
  }
}

//==============================================
//==============================================

void display_mul_gate( struct gate_t *g ) {
  gmp_printf(" %s %Zd %#Zx %d %d %d %d %d\n",
	     g->type, g->root->value, g->root->value, g->ninput, g->linput->id, g->rinput->id, g->noutput, g->output->id);
}

//==============================================
//==============================================

void set_random_input_values( struct circuit_t *c, gmp_randstate_t rnd_state ) {
  int i;
  for (i=0; i<c->nWire; i++) {
    if ( strcmp(c->wire[i]->type,"input")==0) {
      nonzero_rnd_field_elt( c->wire[i]->value, rnd_state);
      //rnd_field_elt( c->wire[i]->value, Fp, rnd_state);
    }
  }
}

//==============================================
//==============================================

void new_gate( unsigned int *ig,
	       unsigned int *iind,
	       unsigned int wcount,
	       char buf[bufSize],
	       char str[strSize],
	       struct circuit_t *c,
	       struct field_t *f,
	       gmp_randstate_t rnd_state) {

  int i, j, k, l, m;
  
  sscanf(buf,"%*s %*s %d %d %d %*s %d %d", &i, &j, &k, &l, &m);

  c->gate[*(ig)] = malloc( sizeof( struct gate_t ));

  c->gate[*(ig)]->type = malloc( 1024 * sizeof( char ));
  strcpy (c->gate[*(ig)]->type, str);


  if( strcmp(str,"mul")==0 ) {

    c->gate[*(ig)]->root = malloc( sizeof( struct field_elt_t ));

    c->gate[*(ig)]->root->field = f;

    mpz_init( c->gate[*(ig)]->root->value );

    nonzero_rnd_field_elt( c->gate[*(ig)]->root, rnd_state );

    //rnd_field_elt( c->gate[*(ig)]->label, c->gate[*(ig)]->field, rnd_state);
  }
  else {

    c->gate[*(ig)]->root = NULL;

  }

  c->gate[*(ig)]->lagrange_poly = NULL;

  //c->gate[*(ig)]->linput = NULL;
  //c->gate[*(ig)]->rinput = NULL;
  //c->gate[*(ig)]->output = NULL;
  
  c->gate[*(ig)]->ninput  = i;

  c->gate[*(ig)]->noutput = l;
  
  int dummy;
  for (dummy=0; dummy<wcount; dummy++) {
    if (c->wire[dummy]->id==j) {
      c->gate[*(ig)]->linput=c->wire[dummy];
      if( strcmp(str,"mul")==0 ) {wire_contribution_to_mul_gate(c->wire[dummy]->left_ctrb, c->gate[*(ig)]);}
    }
    
    if (c->wire[dummy]->id==k) {
      c->gate[*(ig)]->rinput=c->wire[dummy];
      if( strcmp(str,"mul")==0 ) {wire_contribution_to_mul_gate(c->wire[dummy]->right_ctrb, c->gate[*(ig)]);}
    }
    
    if (c->wire[dummy]->id==m) {
      c->gate[*(ig)]->output=c->wire[dummy];
      if( strcmp(str,"mul")==0 ) {wire_contribution_to_mul_gate(c->wire[dummy]->output_ctrb, c->gate[*(ig)]);}
    }
  }

  if( strcmp(str,"mul")==0 ) {
    c->index[*(iind)]=*(ig);
    *(iind)+=1;
  }
  
  *(ig)+=1;

}

//==============================================
//==============================================

void new_wire( unsigned int *iw, char buf[bufSize], char str[strSize], struct circuit_t *c, struct field_t *f ) {

  int i;

  if (strcmp(str,"input")==0) {
    sscanf(buf,"%*s %d", &i);
  }
  else if (strcmp(str,"intermediate")==0) {
    sscanf(buf,"%*s %*s %*d %*d %*d %*s %*d %d", &i);
  }
  else {
    printf("error: str %s not defined... exiting now!\n",str);
    exit(1);
  }

  //c->wire[*(iw)] = malloc( sizeof( struct wire_t ));
  c->wire[*(iw)] = (struct wire_t*) calloc( 1, sizeof( struct wire_t ));
  
  c->wire[*(iw)]->id = i;

  c->wire[*(iw)]->type = malloc( 1024 * sizeof( char ));
  strcpy (c->wire[*(iw)]->type, str);
  
  c->wire[*(iw)]->value = malloc( sizeof( struct field_elt_t )); 
  c->wire[*(iw)]->value->field = f;
  mpz_init( c->wire[*(iw)]->value->value );

  init_wire_contribution( c->wire[*(iw)] );

  //c->wire[*(iw)]->V = NULL;
  //c->wire[*(iw)]->W = NULL;
  //c->wire[*(iw)]->Y = NULL;
  
  *(iw)+=1;
}

//==============================================
//==============================================

void print_wire( struct wire_t *w) {
  printf("\nwire ID | type(input, intermed., output) | value(dec, hex) | gate mul. ID(left contr.) | gate mul. ID(right contr.)\n");
  gmp_printf("wire %d | %s | %Zd %#Zx |", w->id, w->type, w->value, w->value );
}

//==============================================
//==============================================

void print_gate( struct gate_t *g) {
  printf("\ngate type(add, mul,...) | root | nInput | left wire ID | right wire ID | nOutput | output wire ID\n");
  gmp_printf("%s | %Zd | %d | %d | %d | %d | %d\n",
	     g->type, g->root->value, g->ninput, g->linput->id, g->rinput->id, g->noutput, g->output->id);
}

//==============================================
//==============================================

void free_circuit( struct circuit_t *c ) {
  int i;
  
  if (c!=NULL) {
      
    for (i=0; i<c->nWire; i++) {

      if (c->wire[i]!=NULL) {
	free_wire (c->wire[i]);
      }
    }
    free (c->wire);

    for (i=0; i<c->nGate; i++) {
      if (c->gate[i]!=NULL) {
	free_gate (c->gate[i]);
      }
    }
    free (c->gate);
    
    if (c->index!=NULL) {
      free (c->index);
    }

    free (c);

  }

}

