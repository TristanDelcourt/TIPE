#pragma once
#include <gmp.h>

void mod_vect(int* v, int mod, int n1);
void add_vect(int* sum, int* op, int n1);
void div_vect(int* v, int d, int n1);
void sub_vect(int** v, int i, int j, int n1);
void prod_vect(mpz_t prod, mpz_t* z, int n1, system_t s);