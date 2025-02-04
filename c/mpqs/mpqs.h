#pragma once

#include <gmp.h>
#include <stdbool.h>

int** mpqs(mpz_t* z, mpz_t* d, mpz_t N, int pb_len, int* pb, int extra, int s, int delta, bool quiet);