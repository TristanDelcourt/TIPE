#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

bool vectorize_qsieve(mpz_t n, int* v, int pb_len, int* pb){
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

void solve_modular(mpz_t a, int p, mpz_t x1, mpz_t x2){
    /** Finds solution to:
     * x**2 = a mod p
     */

    int t = 0;
    mpz_t temp, p1, a_cpy;
    mpz_inits(temp, p1, a, NULL);
    mpz_set_ui(p1, p);
    mpz_mod_ui(a_cpy, a, p);
    while(t<p){
        mpz_set_ui(temp, t*t);
        mpz_sub(temp, temp, a_cpy);
        if(mpz_legendre(temp, p1) == -1){
            break;
        }
        t++;
    }
    int e = (p+1)/2;
    mpz_set(x1, temp);
    mpz_sqrt(x1, x1);
    mpz_add_ui(x1, x1, t);
    mpz_pow_ui(x1, x1, e);
    mpz_sub(x2, x1, p1);

    mpz_clears(temp, p1, NULL);
}

int** qsieve(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra, int s, bool tests){
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

    float* sinterval = malloc(s*sizeof(float));
    for(int i = 0; i<s; i++){
        sinterval[i] = 0;
    }
    float t = /*sievingT hreshold = log(intervalStart2 − n) − log(maxF actor)*/ 1;

    mpz_t x1, x2;
    mpz_inits(x1, x2, NULL);


    mpz_t temp, p1;
    mpz_inits(temp, p1, NULL);
    for(int i = 1; i < pb_len; i++){
        mpz_add_ui(zi, zi, 1);
        
        solve_modular(N, pb[i], x1, x2);
        
        mpz_set_ui(p1, pb[i]);
        gmp_printf("x1=%Zd, x2=%Zd, p=%Zd\n", x1, x2, p1);

        mpz_mul(temp, x1, x1);
        assert(mpz_congruent_p(temp, N, p1));

        mpz_mul(temp, x2, x2);
        assert(mpz_congruent_p(temp, N, p1));
        
        
        
    /*
        
        bool found = false;
        int* vi = malloc(pb_len*sizeof(int));

        while(!found){

            //Q(x)
            mpz_mul(qx, zi, zi);
            mpz_sub(qx, qx, N);

            found = vectorize_qsieve(qx, vi, pb_len, pb);
        }
        if(!tests){
            printf("\r");
            printf("%.1f%%", (float)i/(pb_len+extra-1)*100);
            fflush(stdout);
        }
        
        v[i] = vi;
        mpz_set(z[i], zi);
    */
    }
   exit(1);

    if(!tests) printf("\n");

    mpz_clears(sqrt_N, zi, qx, NULL);


    return v;
}