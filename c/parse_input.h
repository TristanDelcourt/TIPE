#pragma once
#include <gmp.h>
#include <stdbool.h>

typedef enum {DEFAULT, DIXON, QSIEVE} TYPE;

typedef struct input_s {
    char* output_file;
    int bound, sieving_interval;
    mpz_t N;
    bool quiet;
    TYPE algorithm;
    int extra;
} input_t;

input_t* parse_input(int argc, char** argv);
void free_input(input_t* input);