#ifndef _TOOLBOX_H_
#define _TOOLBOX_H_

void str2mpz( char*, mpz_t );

void mpz2str( mpz_t, char* );

void generate_random_monic_polynomial( struct field_t*, struct element_t*, unsigned int, gmp_randstate_t );

bool is_poly_irreducible( struct element_t* );

#endif /* _TOOLBOX_H_ */
