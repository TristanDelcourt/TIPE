#pragma once

#include <gmp.h>
#include <stdbool.h>

// returns false id n is a trivail compsoite, true otherwise
bool trial_division_primetest(mpz_t n, int max);

// returns a^b mod m
int powmod(mpz_t a, mpz_t b, mpz_t m);
