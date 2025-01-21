#pragma once
#include <gmp.h>

// bruh
bool is_prime(int n);

// calculates pi(n), the number of prime numbers <= n
int pi(int n);

// returns a list of piB first primes
int* primes(int piB, int B);

/** Reduces the factor base of the algorithm, refer to:
 * Quadratic sieve factorisation algorithm
 * Bc. Ondˇrej Vladyka
 * Section 2.3.1 (p.16)
*/
int* prime_base(mpz_t n, int* pb_len, int* primes, int piB);