#pragma once
#include <gmp.h>
#include <stdbool.h>

struct poly_s {
    mpz_t d;
    mpz_t N;

    mpz_t a;
    mpz_t b;
    mpz_t c;

    mpz_t zi;
    mpz_t qx;

    // used to make operations without declaring and freeing everytime
    mpz_t op1, op2, op3;
};

typedef struct poly_s* poly_t;

void get_next_poly(poly_t p);
poly_t init_poly(mpz_t N, int M);
void calc_poly(poly_t p, mpz_t x);
poly_t copy_poly(poly_t p);
void free_poly(poly_t p);