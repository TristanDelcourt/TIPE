#include <gmp.h>
#include <stdbool.h>

bool trial_division_primetest(mpz_t n, int max){
    
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

int powmod(mpz_t a, mpz_t b, mpz_t m){
    mpz_t n;
    mpz_init_set_ui(n, 1);
    
    while (mpz_sgn(b) > 0){
        if (mpz_congruent_ui_p(b, 1, 2)){
            mpz_mul(n, n, a);
            mpz_mod(n, n, m);
        }
        mpz_fdiv_q_ui(b, b, 2);
    }
    return mpz_get_ui(n);
}