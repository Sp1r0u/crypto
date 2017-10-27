#include <stdio.h>
#include <string.h>
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
  c->nIOWire = nInput_wire;
  c->nMul_gate = nMul_gate;
  c->nAdd_gate = nAdd_gate;
  c->nGate = c->nMul_gate + c->nAdd_gate;
  c->index = malloc( nMul_gate * sizeof( unsigned int ));
  
  //c->wire = malloc( (c->nWire) * sizeof( struct wire_t ));
  //c->mul_gate = malloc( (c->nMul_gate) * sizeof( struct gate_t ));
  //c->add_gate = malloc( (c->nAdd_gate) * sizeof( struct gate_t ));
  c->wire = malloc( (c->nWire) * sizeof( struct wire_t* )); //allocates an array of pointers to wire_t
  c->gate = malloc( (c->nGate) * sizeof( struct gate_t* )); //allocates an array of pointers to gate_t

  printf ("\n**** Circuit Initialization ****\n\n");
  printf ("nb wires: %d\n", c->nWire);
  printf ("nb gates: %d\n", c->nGate); 
  printf ("nb mult. gates: %d\n", c->nMul_gate);
  printf ("nb add.gates: %d\n", c->nAdd_gate);

  rewind( fptr );
  //int i_wire = 0;
  //int i_gate = 0;
  //int i_index = 0;
  //int i, j, k, l, m;
  int i, j, k;
  char foo[strSize];
  //int dummy;

  unsigned int *igate  = malloc( sizeof( unsigned int ));
  unsigned int *iwire  = malloc( sizeof( unsigned int ));
  unsigned int *iindex = malloc( sizeof( unsigned int ));

  *(igate)  = 0;
  *(iwire)  = 0;
  *(iindex) = 0;
  
  while ( fgets(buf, sizeof(buf), fptr)!=NULL ) {
    sscanf(buf,"%s",str);

    if( strcmp(str,"input")==0 ) {
      //----------------------
      new_wire(iwire, buf, c, "input");
      //print_wire(c->wire[*(iwire)-1]);
      //----------------------
      //sscanf(buf,"%*s %d", &i);
      //printf("i_wire %d\n",i_wire);
      //c->wire[i_wire] = malloc( sizeof( struct wire_t ));
      //c->wire[i_wire]->id = i;
      //c->wire[i_wire]->type = "input";
      //mpz_init( c->wire[i_wire]->value );
      //init_wire_contribution( c->wire[i_wire] );
      //i_wire+=1;
    }
    
    else if( strcmp(str,"mul")==0 ) {
      
      //sscanf(buf,"%*s %*s %d %d %d %s %d %d", &i, &j, &k, foo, &l, &m);

      sscanf(buf,"%*s %*s %*d %*d %*d %s", foo);
      
      if( strcmp(foo,"out")==0 ) {
	//----------------------
	new_wire(iwire, buf, c, "intermediate");
	//print_wire(c->wire[*(iwire)-1]);
	//----------------------
	//printf("i_wire %d\n",i_wire);
	//c->wire[i_wire] = malloc( sizeof( struct wire_t ));
	//c->wire[i_wire]->id = m;
	//c->wire[i_wire]->type = "intermediate";
	//mpz_init( c->wire[i_wire]->value );
	//init_wire_contribution( c->wire[i_wire] );
	//i_wire+=1;
      }
      //exit(0);
      //printf("i_gate %d\n",i_gate);
      //----------------------
      new_gate(igate, iindex, *(iwire), buf, c, "mul", Fp, rnd_state);
      //print_gate(c->gate[*(igate)-1]);
      //----------------------
      //c->gate[i_gate] = malloc( sizeof( struct gate_t ));
      //c->gate[i_gate]->type = "mul";
      //c->gate[i_gate]->field = Fp;
      //rnd_field_elt( c->gate[i_gate]->label, c->gate[i_gate]->field, rnd_state);
      //c->gate[i_gate]->ninput = i;
      //c->gate[i_gate]->noutput = l;
      //for (dummy=0; dummy<c->nWire; dummy++) {
      //if (c->wire[dummy]->id==j) {
      //c->gate[i_gate]->linput=c->wire[dummy];
      //new_wire_contribution(c->wire[dummy]->leftC, c->gate[i_gate]);
      //}
      //else if (c->wire[dummy]->id==k) {c->gate[i_gate]->rinput=c->wire[dummy];}
      //else if (c->wire[dummy]->id==m) {c->gate[i_gate]->output=c->wire[dummy];}
      //else {printf("error: wire id %d not defined... exiting now!\n",c->wire[dummy]->id); exit(1);}
      //}
      //c->index[i_index]=i_gate;
      //i_index+=1;
      //i_gate+=1;
    }//end "mul"

    else if( strcmp(str,"add")==0 ) {

      //sscanf(buf,"%*s %*s %d %d %d %s %d %d", &i, &j, &k, foo, &l, &m);

      sscanf(buf,"%*s %*s %*d %*d %*d %s", foo);
      
      if( strcmp(foo,"out")==0 ) {
	//----------------------
	new_wire(iwire, buf, c, "intermediate");
	//print_wire(c->wire[*(iwire)-1]);
	//----------------------
	//printf("i_wire %d\n",i_wire);
	//c->wire[i_wire] = malloc( sizeof( struct wire_t ));
	//c->wire[i_wire]->id = m;
	//c->wire[i_wire]->type = "intermediate";
	//mpz_init( c->wire[i_wire]->value );
	//init_wire_contribution( c->wire[i_wire] );
	//i_wire+=1;
      }
      //----------------------
      new_gate(igate, iindex, *(iwire), buf, c, "add", Fp, rnd_state);
      //print_gate(c->gate[*(igate)-1]);
      //----------------------
      //printf("i_gate %d\n",i_gate);
      //c->gate[i_gate] = malloc( sizeof( struct gate_t ));
      //c->gate[i_gate]->type = "add";
      //c->gate[i_gate]->field = Fp;
      //mpz_init( c->gate[i_gate]->label );
      //mpz_set_si( c->gate[i_gate]->label, -99 );
      //c->gate[i_gate]->ninput = i;
      //c->gate[i_gate]->noutput = l;
      //for (dummy=0; dummy<c->nWire; dummy++) {
      //if      (c->wire[dummy]->id==j) {c->gate[i_gate]->linput=c->wire[dummy];}
      //else if (c->wire[dummy]->id==k) {c->gate[i_gate]->rinput=c->wire[dummy];}
      //else if (c->wire[dummy]->id==m) {c->gate[i_gate]->output=c->wire[dummy];}
      //else {printf("error: wire id %d not defined... exiting now!\n",c->wire[dummy]->id); exit(1);}
      //}
      //i_gate+=1;
    }//end "add"

    else if( strcmp(str,"output")==0 ) {
      sscanf(buf,"%*s %d", &i);
      c->nIOWire += 1;
      k=-1;
      for (j=0; j<*(iwire); j++) {
	//if (c->wire[j].id==i) {
	if (c->wire[j]->id==i) {
	  //c->wire[j].type = "output";
	  c->wire[j]->type = "output";
	  k=1;
	}
      }
      if (k<0) {printf("error: couldn't find output wire %d... exiting now!\n",i); exit(1);}
    }//end "output"

    else {printf("error: %s not defined... exiting now!\n",str); exit(1);}
  }
  printf ("nb IO wires: %d\n", c->nIOWire);
  
  //for (i=0; i<c->nWire; i++) {print_wire_contribution_to_mul_gate(c->wire[i]);}
  //display_circuit(c);
  //exit(0);
}

//==============================================
//==============================================

void display_circuit( struct circuit_t *c ) {
  int i;

  printf("\nwire ID | type(input, intermed., output) | value(dec, hex) | gate mul. ID(left contr.) | gate mul. ID(right contr.)\n");
  for (i=0; i<c->nWire; i++) {
    gmp_printf("\n %d\n %s\n %Zd %#Zx\n", c->wire[i]->id, c->wire[i]->type, c->wire[i]->value, c->wire[i]->value );
    print_wire_contribution_to_mul_gate(c->wire[i]);
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
  printf("\ngate type(add, mul,...) | label | nInput | left wire ID | right wire ID | nOutput | output wire ID\n\n");
  for (i=0; i<c->nGate; i++) {
    //gmp_printf("%s %Zd %Zd %d %d %d %d %d\n",
    gmp_printf(" %s | %Zd | %d | %d | %d | %d | %d\n\n",
	       c->gate[i]->type,
	       //c->add_gate[i].field->p,
	       c->gate[i]->label, 
	       c->gate[i]->ninput,
	       c->gate[i]->linput->id,
	       c->gate[i]->rinput->id,
	       c->gate[i]->noutput,
	       c->gate[i]->output->id);  
  }
}

//==============================================
//==============================================

void display_mul_gate( struct circuit_t *c ) {
  int i;
  for (i=0; i<c->nMul_gate; i++) {
    gmp_printf("%s %Zd %#Zx %d %d %d %d %d\n",
	       c->gate[c->index[i]]->type,
	       //c->add_gate[i].field->p,
	       c->gate[c->index[i]]->label,
	       c->gate[c->index[i]]->label,
	       c->gate[c->index[i]]->ninput,
	       c->gate[c->index[i]]->linput->id,
	       c->gate[c->index[i]]->rinput->id,
	       c->gate[c->index[i]]->noutput,
	       c->gate[c->index[i]]->output->id);
  }
}

//==============================================
//==============================================

void set_random_input_values( struct circuit_t *c, struct field_t *Fp, gmp_randstate_t rnd_state ) {
  int i;
  for (i=0; i<c->nWire; i++) {
    if ( strcmp(c->wire[i]->type,"input")==0) {
      rnd_field_elt( c->wire[i]->value, Fp, rnd_state);
    }
  }
}

//==============================================
//==============================================

void new_gate( unsigned int *ig, unsigned int *iind, unsigned int wcount, char buf[bufSize], struct circuit_t *c,
	       char str[strSize], struct field_t *f, gmp_randstate_t rnd_state) {

  int i, j, k, l, m;

  sscanf(buf,"%*s %*s %d %d %d %*s %d %d", &i, &j, &k, &l, &m);

  c->gate[*(ig)] = malloc( sizeof( struct gate_t ));

  c->gate[*(ig)]->type = str;
  c->gate[*(ig)]->field = f;

  if( strcmp(str,"mul")==0 ) {
    rnd_field_elt( c->gate[*(ig)]->label, c->gate[*(ig)]->field, rnd_state);
  }
  else {
    mpz_init( c->gate[*(ig)]->label );
    mpz_set_si( c->gate[*(ig)]->label, -99 );
  }

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

void new_wire( unsigned int *iw, char buf[bufSize], struct circuit_t *c, char str[strSize] ) {

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

  c->wire[*(iw)] = malloc( sizeof( struct wire_t ));
  c->wire[*(iw)]->id = i;
  c->wire[*(iw)]->type = str;

  mpz_init( c->wire[*(iw)]->value );

  init_wire_contribution( c->wire[*(iw)] );
  
  *(iw)+=1;
}

//==============================================
//==============================================

void print_wire( struct wire_t *w) {
  printf("\nwire ID | type(input, intermed., output) | value(dec, hex) | gate mul. ID(left contr.) | gate mul. ID(right contr.)\n");
  gmp_printf("wire %d | %s | %Zd %#Zx |\n", w->id, w->type, w->value, w->value );
}

//==============================================
//==============================================

void print_gate( struct gate_t *g) {
  printf("\ngate type(add, mul,...) | label | nInput | left wire ID | right wire ID | nOutput | output wire ID\n");
  gmp_printf("%s | %Zd | %d | %d | %d | %d | %d\n", g->type, g->label, g->ninput, g->linput->id, g->rinput->id, g->noutput, g->output->id);
}

//==============================================
//==============================================

void print_wire_contribution_to_mul_gate( struct wire_t *w) {

  if (w->left_ctrb->length>0) {
    struct node_t *node = w->left_ctrb->head;
    //printf("\nwire ID %d left ctrb(s) to mult. gates:\n",w->id);
    while( node != NULL ) {
      gmp_printf(" %Zd\n",node->gate->label);
      node = node->next;
    } 
    free( node );
  }
  else {
    //printf("\nwire ID %d NO left ctrb to mult. gates\n",w->id);
    printf(" no left ctrb.\n");
  }
  
  if (w->right_ctrb->length>0) {
    struct node_t *node = w->right_ctrb->head;
    //printf("\nwire ID %d right ctrb(s) to mult. gates\n",w->id);
    while( node != NULL ) {
      gmp_printf(" %Zd\n",node->gate->label);
      node = node->next;
    } 
    free( node );
  }
  else {
    //printf("wire ID %d NO right ctrb to mult. gates\n",w->id);
    printf(" no right ctrb.\n");
  }
  
}
