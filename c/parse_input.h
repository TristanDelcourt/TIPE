#pragma once
#include <gmp.h>

enum Algorithms {
    DIXON,
    QSIEVE
};

typedef struct input_s {
    char* output_file;
    int bound, sieving_interval;
    mpz_t N;
} input_t;

input_t* parse_input(int argc, char** argv);