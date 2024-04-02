#include <gmp.h>
#include <stdbool.h>
#include "big_int.h"

void fermat_naive(mpz_t n, mpz_t a, mpz_t b) {
        
    mpz_t SQRT, r;
    /*
    if (mpz_perfect_power_p(n)){
        mpz_sqrt(a, n);
        mpz_set(b, a);
        return;
    }
    */

    //mpz_inits(SQRT, r);
    mpz_init(SQRT);
    mpz_init(r);

    mpz_sqrtrem(SQRT, r, n);
    mpz_add_ui(SQRT, SQRT, 1);
    mpz_neg(r, r);

    mpz_t u, v;

    //mpz_inits(u, v);
    mpz_init(u);
    mpz_init(v);

    mpz_mul_ui(u, SQRT, 2);
    mpz_add_ui(u, u, 1);
    mpz_set_ui(v, 1);
        
    while (mpz_sgn(r) != 0){

        if (mpz_sgn(r) > 0){
            while (mpz_sgn(r) > 0){
                //gmp_printf("u: %Zd, v: %Zd, r: %Zd\n", u, v, r);
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

    mpz_clears(SQRT, r, u, v, NULL);

}

bool checkgcd(mpz_t n, mpz_t a, mpz_t d){
    mpz_t gcd;
    mpz_init(gcd);
    mpz_gcd(gcd, n, a); // TODO
    if (mpz_cmp_ui(gcd, 1) > 0){
        mpz_set(d, gcd);

        mpz_clear(gcd);
        return true;
    }
    //mpz_clear(gcd);
    return false;
}

bool brent_pollard_rho(mpz_t n, int c, int max, mpz_t d) {
    // using f(x) = x**2 + c

    mpz_t x_1, x_2, product;
    mpz_init_set_ui(x_1, 2);
    mpz_init_set_ui(x_2, 4 + c);
    mpz_init_set_ui(product, 1);

    int range = 1;
    int terms = 0;

    while (terms <= max){
        for(int j = 1; j <= range; j++){
            mpz_mul(x_2, x_2, x_2);
            mpz_add_ui(x_2, x_2, c);
            mpz_mod(x_2, x_2, n);

            mpz_t temp;
            mpz_init(temp);
            mpz_sub(temp, x_1, x_2);
            mpz_mul(product, product, temp);
            mpz_mod(product, product, n);

            terms += 1;
            if (terms % 10 == 0)
                if (!checkgcd(n, product, d))
                    mpz_set_ui(product, 1);
                else
                    //mpz_clears(x_1, x_2, product, temp, NULL);
                    return true;

        }
        mpz_set(x_1, x_2);
        range *= 2;
        for(int j = 1; j <= range; j++){
            mpz_mul(x_2, x_2, x_2);
            mpz_add_ui(x_2, x_2, c);
            mpz_mod(x_2, x_2, n);
        }

    }
    //mpz_clears(x_1, x_2, product, NULL);
    return false;
}

bool pollard_p_minus_1(mpz_t n, int c, int max, mpz_t d){
    mpz_t m;
    mpz_init_set_ui(m, c);

    for(int i = 0; i>max; i++){
        mpz_t a;
        mpz_init_set_ui(a, i);
        powmod(m, m, a, n);
        if(i%10 == 0){
            if (!checkgcd(n, m, d))
                return false;
            else
                return true;
        }
    }
}