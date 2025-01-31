#include "polynomial.h"
#include <gmp.h>
#include <stdlib.h>
#include <assert.h>

#include "../tonellishanks.h"

poly_t init_poly(mpz_t N, int M){
    poly_t p = malloc(sizeof(struct poly_s));

    mpz_inits(p->a, p->b, p->c, p->op1, p->op2, p->op3, p->op4, p->qx, NULL);
    p->done_iter = false;

    // choose value of d according to 2.4.2
    // sqrt(( sqrt(2N))/M )
    mpz_mul_ui(p->op1, N, 2);
    mpz_sqrt(p->op1, p->op1);
    mpz_div_ui(p->op1, p->op1, M);
    mpz_sqrt(p->op1, p->op1);
    mpz_prevprime(p->op1, p->op1);
    
    // get next prime such that (n/p) = 1
    while(mpz_legendre(N, p->op1) != 1){
        mpz_nextprime(p->op1, p->op1);
    }

    mpz_mul(p->a, p->op1, p->op1);

    mpz_t x1, x2;
    mpz_inits(x1, x2, NULL);
    tonelli_shanks_mpz(N, p->op1, x1, x2);

    // t
    mpz_set_ui(p->op2, 1);

    // getting ready to test congruences
    mpz_mul(p->op3, x1, x1);
    mpz_sub(p->op3, p->op3, N);
    mpz_divexact(p->op3, p->op3, p->op1);
    mpz_neg(p->op3, p->op3);
    mpz_mod(p->op3, p->op3, p->op1);

    mpz_mul_ui(p->op4, x1, 2);
    while(mpz_congruent_p(p->op4, p->op3, p->op1) == 0){
        mpz_add_ui(p->op2, p->op2, 1);
        mpz_add(p->op4, p->op4, x1);
        mpz_add(p->op4, p->op4, x1);
    }

    mpz_set(p->b, p->op1);
    mpz_mul(p->b, p->b, p->op2);
    mpz_add(p->b, p->b, x1);

    mpz_mul(p->op2, p->b, p->b);
    assert(mpz_congruent_p(p->op2, N, p->a) != 0);

    mpz_sub(p->c, p->op2, N);
    mpz_divexact(p->c, p->c, p->a);

    return p;
}

