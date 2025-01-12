#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "vector_deprecated.h"

bnvect_t* bnvect_init(unsigned long int n){
    bnvect_t* v = malloc(sizeof(bnvect_t));
    v->n = n;
    v->vi = malloc(n*sizeof(mpz_t));

    for(unsigned long int i = 0; i<n; i++){
        mpz_init(v->vi[i]);
    }
    return v;
}

bnvect_t* bnvect_init_set(unsigned long int n, mpz_t v0){
    bnvect_t* v = malloc(sizeof(bnvect_t));
    v->n = n;
    v->vi = malloc(n*sizeof(mpz_t));

    for(unsigned long int i = 0; i<n; i++){
        mpz_init_set(v->vi[i], v0);
    }
    return v;
}

bnvect_t* bnvect_init_set_si(unsigned long int n, signed long v0){
    bnvect_t* v = malloc(sizeof(bnvect_t));
    v->n = n;
    v->vi = malloc(n*sizeof(mpz_t));

    for(unsigned long int i = 0; i<n; i++){
        mpz_init_set_si(v->vi[i], v0);
    }
    return v;
}

void bnvect_set_si(bnvect_t* v, signed long v0){
    for(unsigned long int i = 0; i<v->n; i++){
        mpz_set_si(v->vi[i], v0);
    }
}

int bnvect_add(bnvect_t* sum, bnvect_t* op1, bnvect_t* op2){
    if(op1->n != op2->n)
        return -1;
    
    for(unsigned long int i = 0; i<op1->n; i++){
        mpz_add(sum->vi[i], op1->vi[i], op2->vi[i]);
    }

    return 0;
    
}

void bnvect_add_ui(bnvect_t* sum, bnvect_t* op1, unsigned long op2){    
    for(unsigned long int i = 0; i<op1->n; i++){
        mpz_add_ui(sum->vi[i], op1->vi[i], op2);
    }
    
}

int bnvect_add_index_ui(bnvect_t* op1, unsigned long op2, unsigned long int i){
    if(i > op1->n -1)
        return -1;

    mpz_add_ui(op1->vi[i], op1->vi[i], op2);

    return 0;
}

int bnvect_sub(bnvect_t* diff, bnvect_t* op1, bnvect_t* op2){
    if(op1->n != op2->n)
        return -1;
    
    for(unsigned long int i = 0; i<op1->n; i++){
        mpz_sub(diff->vi[i], op1->vi[i], op2->vi[i]);
    }

    return 0;
}

void bnvect_sub_ui(bnvect_t* sum, bnvect_t* op1, unsigned long op2){    
    for(unsigned long int i = 0; i<op1->n; i++){
        mpz_sub_ui(sum->vi[i], op1->vi[i], op2);
    }
    
}

void bnvect_mod(bnvect_t* op1, mpz_t mod){
    for(unsigned long int i = 0; i<op1->n; i++){
        mpz_mod(op1->vi[i], op1->vi[i], mod);
    }
}

void bnvect_mod_ui(bnvect_t* op1, unsigned long mod){
    for(unsigned long int i = 0; i<op1->n; i++){
        mpz_mod_ui(op1->vi[i], op1->vi[i], mod);
    }
}

int bnvect_assign(bnvect_t* op1, mpz_t op2, unsigned long int i){
    if(i > op1->n -1)
        return -1;

    mpz_set(op1->vi[i], op2);

    return 0;
}

int bnvect_assign_si(bnvect_t* op1, signed long op2, unsigned long int i){
    if(i >= op1->n)
        return -1;

    mpz_set_si(op1->vi[i], op2);

    return 0;
}

int bnvect_init_get(mpz_t out, bnvect_t* op1, unsigned long int i){
    if(i >= op1->n)
        return -1;

    mpz_init_set(out, op1->vi[i]);

    return 0;

}

void bnvect_print(bnvect_t* op1){
    printf("(\n");
    for(unsigned long int i = 0; i < op1->n -1; i++){
        gmp_printf(" %Zd,\n", op1->vi[i]);
    }
    gmp_printf(" %Zd\n", op1->vi[op1->n-1]);
    printf(")\n");
}
