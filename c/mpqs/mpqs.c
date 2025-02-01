#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "polynomial.h"
#include "../system.h"
#include "../tonellishanks.h"

bool vectorize_mpqs(mpz_t n, int* v, int pb_len, int* pb){
    /** Attemps naive factorisation to 'n' with the primes in
     * the prime base 'pb' and putting the result into 'v', vector of powers of
     * the primes in the prime base
     * If it succeeds, returns true, otherwise, returns false
    */
    for(int i = 0; i<pb_len; i++){
        v[i] = 0;
    }
    if(mpz_sgn(n)<0){
        v[pb_len] = 1;
        mpz_neg(n, n);
    }
    else{
        v[pb_len] = 0;
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

float* prime_logs_mpqs(int* pb, int pb_len){
    float* plogs = malloc(pb_len*sizeof(float));
    
    for(int i = 0; i<pb_len; i++){
        plogs[i] = log2(pb[i]);
    }

    return plogs;
}

int calculate_threshhold_mpqs(mpz_t sqrt_N, int s, int* pb, int pb_len){
    
    mpz_t qstart;
    mpz_init_set_ui(qstart, s);
    mpz_mul(qstart, qstart, sqrt_N);

    int t = mpz_sizeinbase(qstart, 2) - (int) log2(pb[pb_len-1]);
    mpz_clear(qstart);
    return t;
}

int** mpqs(mpz_t* z, mpz_t* d, mpz_t N, int pb_len, int* pb, int extra, int s, bool quiet){
    /** Gets pb_len+extra zis that are b-smooth, definied at:
     * Quadratic sieve factorisation algorithm
     * Bc. Ondˇrej Vladyka
     * Definition 1.11 (p.5)
     */

    //ceil(sqrt(n))
    mpz_t sqrt_N;
    mpz_init(sqrt_N);
    mpz_sqrt(sqrt_N, N);
    mpz_add_ui(sqrt_N, sqrt_N, 1);

    mpz_t x;
    mpz_init(x);
    poly_t Q = init_poly(N, s);

    int** v = malloc((pb_len+extra)*sizeof(int*));
    for(int i = 0; i<pb_len+extra; i++){
        v[i] = malloc((pb_len+1)*sizeof(int*)); // +1 for -1
    }
    float* sinterval = malloc(2*s*sizeof(float));
    float* plogs = prime_logs_mpqs(pb, pb_len);

    
    // TESTS
    mpz_t temp;
    mpz_init(temp);
    // END TESTS
    

    int* r = malloc(pb_len*sizeof(int));
    int* x1 = malloc(pb_len*sizeof(int));
    int* x2 = malloc(pb_len*sizeof(int));

    // find solution for 2
    

    int sol1, sol2;
    for(int i = 1; i < pb_len; i++){
        tonelli_shanks_ui(N, pb[i], &sol1, &sol2);
        r[i] = sol1;
    }

    mpz_t g, m, n, pi;
    mpz_inits(g, m, n, pi, NULL);

    int relations_found = 0;
    clock_t start;
    start = clock();
    int tries = 0;
    while(relations_found < pb_len + extra){

        // for 2
        mpz_set_ui(temp, 0);
        calc_poly(Q, temp);
        x1[0] = 0;
        if(mpz_divisible_ui_p(Q->qx, 2) == 0) x1[0] = 1;

        //others
        for(int i = 1; i<pb_len; i++){
            mpz_set_ui(pi, pb[i]);
            mpz_gcdext(g, m, n, Q->a, pi);
            assert(mpz_cmp_ui(g, 1) == 0);

            mpz_set_ui(temp, r[i]);
            mpz_sub(temp, temp, Q->b);
            mpz_mul(temp, temp, m);
            mpz_mod(temp, temp, pi);

            x1[i] = mpz_get_ui(temp);

            calc_poly(Q, temp);
            assert(mpz_divisible_ui_p(Q->qx, pb[i]) != 0);

            mpz_set_ui(temp, pb[i]);
            mpz_sub_ui(temp, temp, r[i]);
            mpz_sub(temp, temp, Q->b);
            mpz_mul(temp, temp, m);
            mpz_mod(temp, temp, pi);

            x2[i] = mpz_get_ui(temp);

            calc_poly(Q, temp);
            assert(mpz_divisible_ui_p(Q->qx, pb[i]) != 0);

            
            //realign sieving interval to [-s, s]
            int k = (x1[i] + s)/pb[i];
            x1[i] -= k * pb[i];
            x1[i] += s;

            k = (x2[i] + s)/pb[i];
            x2[i] -= k * pb[i];
            x2[i] += s;

            mpz_set_si(temp, -s);
            mpz_add_ui(temp, temp, x1[i]);
            calc_poly(Q, temp);
            assert(mpz_divisible_ui_p(Q->qx, pb[i]) != 0);
        }

        for(int i = 0; i<2*s; i++){
            sinterval[i] = 0;
        }

        // sieve for 2
        while(x1[0]<2*s){
            sinterval[x1[0]] += plogs[0];
            x1[0] += pb[0];
        }

        // sieve other primes
        for(int i = 1; i < pb_len; i++){

            while(x1[i]<2*s){
                sinterval[x1[i]] += plogs[i];
                x1[i] += pb[i];
            }

            while(x2[i]<2*s){
                sinterval[x2[i]] += plogs[i];
                x2[i] += pb[i];
            }
        }

        int t = calculate_threshhold_mpqs(sqrt_N, s, pb, pb_len);
        //printf("t = %d\n", t);

        bool found;
        bool update_time = false;
        for(int i = 0; i<2*s && relations_found < pb_len + extra; i++){
            if(sinterval[i] > t){
                tries++;
                mpz_set_si(x, -s);
                mpz_add_ui(x, x, i);
                calc_poly(Q, x);
                
                found = vectorize_mpqs(Q->qx, v[relations_found], pb_len, pb);

                if(found){
                    mpz_set(z[relations_found], Q->zi);
                    mpz_set(d[relations_found], Q->d);
                    relations_found++;
                    update_time = true;
                    found = false;
                    if(!quiet){
                        printf("\r");
                        printf("%.1f%% | %.1f%%", (float)relations_found/(pb_len+extra)*100, (float)relations_found/tries*100);
                        fflush(stdout);
                    }
                }
            }
        }

        if(update_time && !quiet) printf(" (~%.0fs left)        " , (double)(clock() - start)/CLOCKS_PER_SEC/relations_found*((pb_len+extra - relations_found)));
        get_next_poly(Q);
    }

    if(!quiet) printf("\n");

    mpz_clears(sqrt_N, temp, g, m, n, pi, x, NULL);
    free(x1);
    free(x2);
    free(r);
    free(sinterval);
    free(plogs);
    free_poly(Q);

    return v;
}