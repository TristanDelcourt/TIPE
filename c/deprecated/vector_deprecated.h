#pragma once
#include <gmp.h> 

typedef struct bignumvector {
    unsigned long int n;
    mpz_t* vi;
} bnvect_t;

bnvect_t* bnvect_init(unsigned long int n);

bnvect_t* bnvect_init_set(unsigned long int n, mpz_t v0);
bnvect_t* bnvect_init_set_si(unsigned long int n, signed long v0);

void bnvect_set_si(bnvect_t* v, signed long v0);

int bnvect_add(bnvect_t* sum, bnvect_t* op1, bnvect_t* op2);
void bnvect_add_ui(bnvect_t* sum, bnvect_t* op1, unsigned long op2);
int bnvect_add_index_ui(bnvect_t* op1, unsigned long op2, unsigned long int i);

int bnvect_sub(bnvect_t* diff, bnvect_t* op1, bnvect_t* op2);
void bnvect_sub_ui(bnvect_t* sum, bnvect_t* op1, unsigned long op2);    

void bnvect_mod(bnvect_t* op1, mpz_t mod);
void bnvect_mod_ui(bnvect_t* op1, unsigned long mod);

int bnvect_assign(bnvect_t* op1, mpz_t op2, unsigned long int i);
int bnvect_assign_si(bnvect_t* op1, signed long op2, unsigned long int i);

int bnvect_init_get(mpz_t out, bnvect_t* op1, unsigned long int i);

void bnvect_print(bnvect_t* op1);
