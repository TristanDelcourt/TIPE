#pragma once
#include <gmp.h>
#include <stdbool.h>

int calculate_threshhold_mpqs(mpz_t sqrt_N, int s, int* pb, int pb_len, int delta);
float* prime_logs_mpqs(int* pb, int pb_len);
bool vectorize_mpqs(mpz_t n, int* v, int pb_len, int* pb);
bool already_added(mpz_t zi, mpz_t* z, int relations_found);