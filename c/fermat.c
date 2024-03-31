#include <gmp.h>

void fermat_naive(mpz_t n, mpz_t a, mpz_t b) {
        
    mpz_t SQRT, r;

    mpz_inits(SQRT, r);
    mpz_sqrtrem(SQRT, r, n);

    mpz_t u, v;
    mpz_inits(u, v);

    mpz_mul_ui(u, SQRT, 2);
    mpz_add_ui(u, u, 1);
    
    mpz_set_ui(v, 1);
    
    while (mpz_sgn(r) != 0){
        
        if (mpz_sgn(r) > 0){
            while (mpz_sgn(r) > 0){
                mpz_sub(r, r, v);
                mpz_add_ui(v, v, 2);
            }
        }

        if (mpz_sgn(r) < 0){
            mpz_add(r, r, u);
            mpz_add_ui(u, u, 2);
        }
    
    }

    mpz_add(a, u, v);
    mpz_sub_ui(a, a, 2);
    mpz_divexact_ui(a, a, 2);

    mpz_sub(b, u, v);
    mpz_divexact_ui(b, b, 2);
}