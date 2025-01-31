#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "../system.h"
#include "../tonellishanks.h"

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

float* prime_logs(int* pb, int pb_len){
    float* plogs = malloc(pb_len*sizeof(float));
    
    for(int i = 0; i<pb_len; i++){
        plogs[i] = log2(pb[i]);
    }

    return plogs;
}

int calculate_threshhold(mpz_t N, mpz_t sqrt_N, int s, int loop_number, int* pb, int pb_len){
    
    mpz_t qstart;
    mpz_init_set_ui(qstart, s);
    mpz_mul_ui(qstart, qstart, loop_number);
    mpz_add(qstart, qstart, sqrt_N);
    mpz_mul(qstart, qstart, qstart);
    mpz_sub(qstart, qstart, N);

    int t = mpz_sizeinbase(qstart, 2) - (int) log2(pb[pb_len-1]);
    mpz_clear(qstart);
    return t;
}

int** qsieve(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra, int s, bool quiet){
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

    mpz_t zi;
    mpz_init_set(zi, sqrt_N);
    mpz_t qx;
    mpz_init(qx);

    int** v = malloc((pb_len+extra)*sizeof(int*));
    for(int i = 0; i<pb_len+extra; i++){
        v[i] = malloc(pb_len*sizeof(int*));
    }
    float* sinterval = malloc(s*sizeof(float));
    float* plogs = prime_logs(pb, pb_len);

    
    // TESTS
    mpz_t temp;
    mpz_init(temp);
    // END TESTS
    

    int* x1 = malloc(pb_len*sizeof(int));
    int* x2 = malloc(pb_len*sizeof(int));

    // find solution for 2
    mpz_set(temp, sqrt_N);
    mpz_mul(temp, temp, temp);
    mpz_sub(temp, temp, N);
    x1[0] = 0;
    if(mpz_divisible_ui_p(temp, 2) != 0) x1[0] = 1;

    int sol1, sol2;
    for(int i = 1; i < pb_len; i++){

            tonelli_shanks_ui(N, pb[i], &sol1, &sol2);
            x1[i] = sol1;
            x2[i] = sol2;

            // change solution from x² = n [p] to (sqrt(N) + x)² = n [p]
            mpz_set_ui(temp, x1[i]);
            mpz_sub(temp, temp, sqrt_N);
            mpz_mod_ui(temp, temp, pb[i]);

            x1[i] = mpz_get_ui(temp);

            mpz_set_ui(temp, x2[i]);
            mpz_sub(temp, temp, sqrt_N);
            mpz_mod_ui(temp, temp, pb[i]);

            x2[i] = mpz_get_ui(temp);
    }
    mpz_clear(temp);

    int loop_number = 0;
    int relations_found = 0;
    int tries = 0;
    while(relations_found < pb_len + extra){
        
        for(int i = 0; i<s; i++){
            sinterval[i] = 0;
        }

        // sieve for 2
        while(x1[0]<s){
            sinterval[x1[0]] += plogs[0];
            x1[0] += pb[0];
        }
        x1[0] = x1[0] - s;

        // sieve other primes
        for(int i = 1; i < pb_len; i++){

            while(x1[i]<s){
                sinterval[x1[i]] += plogs[i];
                x1[i] += pb[i];
            }

            while(x2[i]<s){

                sinterval[x2[i]] += plogs[i];
                x2[i] += pb[i];
            }

            //next interval
            x1[i] = x1[i] - s;
            x2[i] = x2[i] - s;
        }

        int t = calculate_threshhold(N, sqrt_N, s, loop_number, pb, pb_len);
        //printf("t = %d\n", t);

        bool found;
        for(int i = 0; i<s && relations_found < pb_len + extra; i++){
            if(sinterval[i] > t){
                tries++;

                // zi = sqrt(n) + x where x = s*loopnumber + i
                mpz_set_ui(zi, s);
                mpz_mul_ui(zi, zi, loop_number);
                mpz_add_ui(zi, zi, i);
                mpz_add(zi, zi, sqrt_N);

                // qx = zi**2 - N
                mpz_mul(qx, zi, zi);
                mpz_sub(qx, qx, N);
                
                found = vectorize_qsieve(qx, v[relations_found], pb_len, pb);
                                
                if(found){
                    mpz_set(z[relations_found], zi);
                    relations_found++;
                    found = false;
                    if(!quiet){
                        printf("\r");
                        printf("%.1f%% | %.1f%%", (float)relations_found/(pb_len+extra)*100, (float)relations_found/tries*100);
                        fflush(stdout);
                    }
                }
            }
        }
        loop_number++;
    }

    if(!quiet) printf("\n");

    mpz_clears(sqrt_N, zi, qx, NULL);
    free(x1);
    free(x2);
    free(sinterval);
    free(plogs);

    return v;
}