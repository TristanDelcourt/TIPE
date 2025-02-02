#include "polynomial.h"
#include <gmp.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "../tonellishanks.h"

void calc_coefficients(poly_t p){
    mpz_mul(p->a, p->d, p->d);

    mpz_t x1, x2;
    mpz_inits(x1, x2, NULL);
    tonelli_shanks_mpz(p->N, p->d, x1, x2);    

    // getting ready for congruence solve for raising solution
    mpz_mul_ui(p->op1, x1, 2);

    mpz_mul(p->op2, x1, x1);
    mpz_sub(p->op2, p->op2, p->N);
    mpz_divexact(p->op2, p->op2, p->d);
    mpz_neg(p->op2, p->op2);
    mpz_mod(p->op2, p->op2, p->d);

    mpz_t g, n, m;
    mpz_inits(g, n, m, NULL);
    mpz_gcdext(g, n, m, p->d, p->op1);
    assert(mpz_cmp_ui(g, 1) == 0);
    mpz_mul(p->op1, p->op2, m); // t
    mpz_clears(g, n, m, NULL);

    mpz_set(p->b, p->d);
    mpz_mul(p->b, p->b, p->op1);
    mpz_add(p->b, p->b, x1);

    mpz_mul(p->op1, p->b, p->b);
    assert(mpz_congruent_p(p->op1, p->N, p->a) != 0);

    mpz_sub(p->c, p->op1, p->N);
    mpz_divexact(p->c, p->c, p->a);

    mpz_clears(x1, x2, NULL);
}

void get_next_poly(poly_t p){
    mpz_nextprime(p->d, p->d);
    while(mpz_legendre(p->N, p->d) != 1){
        mpz_nextprime(p->d, p->d);
    }
    calc_coefficients(p);
}

poly_t init_poly(mpz_t N, int M){
    poly_t p = malloc(sizeof(struct poly_s));

    mpz_inits(p->d, p->N, p->a, p->b, p->c, p->op1, p->op2, p->op3, p->zi, p->qx, NULL);
    mpz_set(p->N, N);

    // choose value of d according to 2.4.2
    // sqrt( (sqrt(2N))/M )
    mpz_mul_ui(p->op1, N, 2);
    mpz_sqrt(p->op1, p->op1);
    mpz_div_ui(p->op1, p->op1, M);
    mpz_sqrt(p->op1, p->op1);
    mpz_prevprime(p->d, p->op1);
    
    // get next prime such that (n/p) = 1
    while(mpz_legendre(N, p->d) != 1){
        mpz_nextprime(p->d, p->d);
    }

    calc_coefficients(p);
    return p;
}

void calc_poly(poly_t p, mpz_t x){
    mpz_mul(p->zi, p->a, x);
    mpz_add(p->zi, p->zi, p->b);
    
    mpz_mul(p->qx, x, x);
    mpz_mul(p->qx, p->qx, p->a);

    mpz_mul(p->op1, p->b, x);
    mpz_mul_ui(p->op1, p->op1, 2);
    mpz_add(p->qx, p->qx, p->op1);

    mpz_add(p->qx, p->qx, p->c);

}

void free_poly(poly_t p){
    mpz_clears(p->d, p->N, p->a, p->b, p->c, p->op1, p->op2, p->op3, p->zi, p->qx, NULL);
    free(p);
}

poly_t copy_poly(poly_t p){
    poly_t cpy = malloc(sizeof(struct poly_s));

    mpz_inits(cpy->d, cpy->N, cpy->a, cpy->b, cpy->c, cpy->op1, cpy->op2, cpy->op3, cpy->zi, cpy->qx, NULL);
    
    mpz_set(cpy->d, p->d);
    mpz_set(cpy->N, p->N);

    mpz_set(cpy->a, p->a);
    mpz_set(cpy->b, p->b);
    mpz_set(cpy->c, p->c);

    return cpy;
}