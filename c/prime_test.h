#pragma once

#include <gmp.h>
#include <stdbool.h>

// returns true if n is a probable prime, false otherwise
bool probable_prime_test(mpz_t n);