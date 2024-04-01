#include <gmp.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include "big_int.h"

bool trial_division_test(mpz_t n, int max){
    
    for(int d = 2; d < 4; d++) {
        if (mpz_congruent_ui_p(n, 0, d))
            return false;
    }
    int d = 5;
    int add = 2;
    
    while (d<=max && mpz_cmp_d(n, d*d) >= 0) {
        if(mpz_congruent_ui_p(n, 0, d))
            return false;
        d += add;
        add = 6 - add;
    }

    return true;
}

bool fermat_primality_test(mpz_t n){

    for(int i = 3; i<7; i+=2){
        mpz_t a;
        mpz_init_set_ui(a, i);

        mpz_t n_minus_1;
        mpz_init(n_minus_1);
        mpz_sub_ui(n_minus_1, n, 1);
        
        mpz_t out;
        mpz_init(out);
        powmod(out, a, n_minus_1, n);

        if(mpz_cmp_ui(out, 1) != 0){
            return false;
        }
    }

    return true;
}

void factor_out_powers_of_2(mpz_t n, mpz_t s, mpz_t d){
    mpz_set_ui(s, 0);
    mpz_set(d, n);

    while(mpz_congruent_ui_p(d, 0, 2)){
        mpz_fdiv_q_ui(d, d, 2);
        mpz_add_ui(s, s, 1);
    }
}

bool miller_rabin_primality_test(mpz_t n, int k){

    mpz_t n_minus_1;
    mpz_init(n_minus_1);
    mpz_sub_ui(n_minus_1, n, 1);
    mpz_t s, d;
    mpz_init(s);
    mpz_init(d);
    factor_out_powers_of_2(n_minus_1, s, d);

    mpz_t a;
    mpz_init(a);
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));

    for(int i = 0; i < k; i++){

        // Generate random number between 2 and n-2
        mpz_set(a, n);
        mpz_sub_ui(a, a, 3);
        mpz_urandomm(a, state, a);
        mpz_add_ui(a, a, 2);

        mpz_t x;
        mpz_init(x);
        powmod(x, a, d, n);

        mpz_t y;
        mpz_init(y);

        for(int j = 0; mpz_cmp_ui(s, j) > 0; j++){

            mpz_powm_ui(y, x, 2, n);

            if(mpz_cmp_ui(y, 1) == 0 && mpz_cmp_ui(x, 1) != 0 && mpz_cmp(x, n_minus_1) != 0){
                return false;
            }

            mpz_set(x, y);
        }

        if(mpz_cmp_ui(y, 1) != 0){
            return false;
        }

    }
    
    return true;
}

bool probable_prime_test(mpz_t n){
    printf("-------- Primetesting n --------\n");

    // Trial division test
    if(!trial_division_test(n, 1000000)){
        printf("Trivial division : composite\n");
        printf("-------- n is composite --------\n");

        return false;
    }
    printf("Trivial division : prime\n");

    // Fermat primality test
    if(!fermat_primality_test(n)){
        printf("Fermat primality : composite\n");
        printf("-------- n is composite --------\n");
        return false;
    }
    printf("Fermat primality : prime\n");

    // Miller-Rabin primality test
    if(!miller_rabin_primality_test(n, 1000)){
        printf("Miller-Rabin     : composite\n");
        printf("-------- n is composite --------\n");
        return false;
    }
    printf("Miller-Rabin     : prime\n");


    printf("---------- n is prime ----------\n");
    return true;
}