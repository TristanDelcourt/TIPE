#include <gmp.h>
#include <assert.h>

unsigned modulo(int value, unsigned m) {
    int mod = value % (int)m;
    if (mod < 0) {
        mod += m;
    }
    return mod;
}

void mod_vect(int* v, int mod, int n1){
    for(int i = 0; i<n1; i++){
        v[i] = modulo(v[i], mod);
    }
}

void add_vect(int* sum, int* op, int n1){
    for(int i = 0; i<n1; i++){
        sum[i] += op[i];
    }
}


void div_vect(int* v, int d, int n1){
    for(int i = 0; i<n1; i++){
        assert(v[i]%d == 0);
        v[i] /= d;
    }
}

void sub_vect(int** v, int i, int j, int n1){
    for(int k = 0; k<n1; k++){
        v[i][k] = v[i][k] - v[j][k];
    }
}

void prod_vect(mpz_t prod, mpz_t* v, int n1){
    mpz_set_ui(prod, 1);
    for(int i = 0; i<n1; i++){
        mpz_mul(prod, prod, v[i]);
    }
}