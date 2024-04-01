#include <gmp.h>
#include <stdio.h> 
#include "prime_test.h"
#include "fermat.h"

int main(){

    mpz_t n;
    mpz_init(n);
    gmp_fscanf(stdin, "%Zd", n);

    /*
    if(!trial_division_primetest(n, 1000000)){
        printf("n is trivially factorable\n");
        return 0;
    }
    printf("Passed trial division test\n");


    for(int i = 3; i<7; i+=2){
        mpz_t a;
        mpz_init_set_ui(a, i);

        mpz_t n_minus_1;
        mpz_init(n_minus_1);
        mpz_sub_ui(n_minus_1, n, 1);
        
        if(powmod(a, n_minus_1, n) == 1){
            printf("n is trivially factorable\n");
            return 0;
        }
    }

    printf("Passed probable prime test\n");
    */

    mpz_t a, b;
    mpz_init(a);
    mpz_init(b);


    //fermat_naive(n, a, b);
    int c = 1;
    while(1){
        printf("c: %d\n", c);
        if (brent_pollard_rho(n, c, 10000, a)){
            if (mpz_divisible_p(n, a)){
                mpz_divexact(b, n, a);
                gmp_printf("n = %Zd * %Zd\n", a, b);
                return 0;
            }
            c += 1;
        }
    }

    
}