#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "vector.h"

void print_list_vect(bnvect_t** l, int n){
    for(int i = 0; i<n; i++){
        bnvect_print(l[i]);
        printf("\n");
    }
}


bnvect_t** init_list(int n){
    bnvect_t** l = malloc(n*sizeof(bnvect_t*));
    return l;
}

bool factorise(mpz_t n, bnvect_t* vi, int piB){
    bnvect_set_si(vi, 0);
    mpz_t n_prime;

    for(int d = 2; d < 4; d++) {
        while (mpz_congruent_ui_p(n, 0, d)){
            mpz_divexact_ui(n, n, d);
            bnvect_add_index_ui(vi, 1, d-2);
        }

    }
    int d = 5;
    int add = 2;
    int i = 2;
    while (mpz_cmp_ui(n, 1)) {
        if(i>=piB){
            return false;
        }

        while (mpz_congruent_ui_p(n, 0, d)){
            mpz_divexact_ui(n, n, d);
            bnvect_add_index_ui(vi, 1, i);
        }
        d += add;
        add = 6 - add;
        i++;
    }

    return true;
}

void get_zi(bnvect_t* vi, mpz_t zi, mpz_t N, int piB){
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    bool found = false;

    while(!found){
        mpz_urandomm(zi, rstate, N);
        mpz_t zi_2;
        mpz_mul(zi_2, zi, zi);
        mpz_mod(zi_2, zi_2, N);
        //gmp_printf("%Zd\n", zi_2);
        found = factorise(zi_2, vi, piB);
    }
}


void get_all_zi(){
    return;
}

void main(int agrc, char**argv){
    int B = 10000;
    int piB = 168;
    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    bnvect_t* v = bnvect_init(piB);
    
    mpz_t zi;
    mpz_init(zi);

    get_zi(v, zi, N, piB);

    /*
    mpz_init_set_ui(zi, 16853);
    mpz_mul(zi, zi, zi);
    mpz_mod_ui(zi, zi, 20382493);
    
    bool out = factorise(zi, v, piB);
    */
    bnvect_print(v);

}