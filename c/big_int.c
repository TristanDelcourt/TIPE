#include <gmp.h>

void powmod(mpz_t n, mpz_t a, mpz_t b, mpz_t m){
    mpz_t a_copy, b_copy;
    mpz_init_set(a_copy, a);
    mpz_init_set(b_copy, b);

    mpz_set_ui(n, 1);
    while (mpz_sgn(b_copy) > 0){
        if (mpz_congruent_ui_p(b_copy, 1, 2)){
            mpz_mul(n, n, a_copy);
            mpz_mod(n, n, m);
        }
        mpz_fdiv_q_ui(b_copy, b_copy, 2);
        mpz_mul(a_copy, a_copy, a_copy);
        mpz_mod(a_copy, a_copy, m);
    }
}