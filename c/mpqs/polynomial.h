#pragma once
#include <gmp.h>
#include <stdbool.h>

struct poly_s {
    mpz_t a;
    mpz_t b;
    mpz_t c;
    mpz_t qx;

    // used to make operations without declaring and freeing everytime
    mpz_t op1, op2, op3, op4;
    bool done_iter;
};

typedef struct poly_s* poly_t;

poly_t init_poly(mpz_t N, int M);