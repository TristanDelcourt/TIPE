#pragma once

#include <gmp.h>

void fermat_naive(mpz_t n, mpz_t a, mpz_t b);

bool brent_pollard_rho(mpz_t n, int c, int max, mpz_t d);