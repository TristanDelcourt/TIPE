#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>

bool vectorize(mpz_t n, int* v, int pb_len, int* pb){
    /** Attemps naive factorisation to 'n' with the primes in
     * the prime base 'pb' and putting the result into 'v', vector of powers of
     * the primes in the prime base
     * If it succeeds, returns true, otherwise, returns false
    */
    for(int i = 0; i<pb_len; i++){
        v[i] = 0;
    }
    
    for(int i = 0; i<pb_len && (mpz_cmp_ui(n, 1) != 0); i++){
        while (mpz_divisible_ui_p(n, pb[i])){
            v[i]++;
            mpz_divexact_ui(n, n, pb[i]);
        }
    }

    if(mpz_cmp_ui(n, 1) == 0)
        return true;
    return false;
}

int** qsieve(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra, bool tests){
    /** Gets pb_len+extra zis such that their product will simplify our searach of
     * a B-smooth relation, definied at:
     * Quadratic sieve factorisation algorithm
     * Bc. Ondˇrej Vladyka
     * Definition 1.11 (p.5)
     */

    //ceil(sqrt(n))
    mpz_t sqrt_N;
    mpz_init(sqrt_N);
    mpz_sqrt(sqrt_N, N);
    mpz_add_ui(sqrt_N, sqrt_N, 1);

    mpz_t zi;
    mpz_init_set(zi, sqrt_N);
    mpz_t qx;
    mpz_init(qx);

    int** v = malloc((pb_len+extra)*sizeof(int*));

    for(int i = 0; i < pb_len+extra; i++){
        bool found = false;
        int* vi = malloc(pb_len*sizeof(int));

        while(!found){
            mpz_add_ui(zi, zi, 1);

            //Q(x)
            mpz_mul(qx, zi, zi);
            mpz_sub(qx, qx, N);

            found = vectorize(qx, vi, pb_len, pb);
        }
        if(!tests){
            printf("\r");
            printf("%f%%", (float)i/(pb_len+extra-1)*100);
            fflush(stdout);
        }
        
        v[i] = vi;
        mpz_set(z[i], zi);
    }
    if(!tests) printf("\n");

    mpz_clears(sqrt_N, zi, qx, NULL);


    return v;
}