#pragma once
#include <gmp.h>
#include <stdbool.h>

int** qsieve(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra, int s, bool tests);