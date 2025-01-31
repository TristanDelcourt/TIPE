#pragma once

#include <gmp.h>

void tonelli_shanks_ui(mpz_t n, int p, int* x1, int* x2);
void tonelli_shanks_mpz(mpz_t a, mpz_t p, mpz_t x1, mpz_t x2);