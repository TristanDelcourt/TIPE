#pragma once
#include <gmp.h>
#include "polynomial.h"
#include <sys/time.h>

struct sieve_arg_s {
    int* pb;
    int pb_len;
    int extra;
    int* r;
    float* plogs;
    int s;
    int t;
    int* relations_found;
    int** v;
    bool quiet;
    mpz_t* z;
    mpz_t* d;
    poly_t Qinit;
    struct timeval begin;
};
typedef struct sieve_arg_s sieve_arg_t;

bool already_added(mpz_t zi, mpz_t* z, int relations_found);
void* sieve_100_polys (void* args);
int** parallel_mpqs(mpz_t* z, mpz_t* d, mpz_t N, int pb_len, int* pb, int extra, int s, bool quiet);