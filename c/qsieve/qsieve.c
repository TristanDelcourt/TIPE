#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

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

void TonelliShanks(mpz_t a, int* p, int* x1, int* x2, int j){
    //Solve the equation: x^2 ≡ a (mod p)
        
        
    //1. Factor `p - 1` into `2^S * Q` where Q is odd.
    int Q = p[j] - 1;
    int S = 0;
    while(Q % 2 == 0){
        S += 1;
        Q /= 2;
    }

    // 2. Find a NR(p).
    mpz_t y, pj;
    mpz_init_set_ui(y, 2);
    mpz_init_set_ui(pj, p[j]);

    while(mpz_legendre(y, pj) != -1){
        mpz_add_ui(y, y, 1);
    }

    // 3. Calculate the four quantities.
    mpz_t R, c, t, E;
    mpz_inits(R, c, t, E, NULL);

    mpz_powm_ui(R, a, (Q+1)/2, pj);
    mpz_powm_ui(c, y, Q, pj);
    mpz_powm_ui(t, a, Q, pj);
    mpz_set_ui(E, S);

    mpz_t temp, b;
    mpz_inits(temp, b, NULL);

    while(mpz_cmp_ui(t, 1) != 0){
        int i = 1;
        while(mpz_cmp_ui(E, i) != 0){
            mpz_set_ui(temp, 1);
            mpz_mul_2exp(temp, temp, i);
            mpz_powm(temp, t, temp, pj);
            if(mpz_cmp_ui(temp, 1) == 0){
                break;
            }
            i++;
        }

        // b = pow(c, 2 ** (E - i - 1), p)
        mpz_set(temp, E);
        mpz_sub_ui(temp, temp, i);
        mpz_sub_ui(temp, temp, 1);
        mpz_ui_pow_ui(temp, 2, mpz_get_ui(temp));
        mpz_powm(b, c, temp, pj);

        // R = R * b % p
        mpz_mul(R, R, b);
        mpz_mod_ui(R, R, p[j]);

        // c = pow(b, 2, p)
        mpz_pow_ui(c, b, 2);
        mpz_mod_ui(c, c, p[j]);

        // t = c * t % p
        mpz_mul(t, c, t);
        mpz_mod_ui(t, t, p[j]);

        mpz_set_ui(E, i);
    }

    x1[j] = mpz_get_ui(R);
    x2[j] = p[j] - mpz_get_ui(R);

    mpz_clears(y, pj, R, c, t, E, temp, b, NULL);
}

float* prime_logs(int* pb, int pb_len){
    float* plogs = malloc(pb_len*sizeof(float));
    
    for(int i = 0; i<pb_len; i++){
        plogs[i] = log2(pb[i]);
    }

    return plogs;
}

int calculate_threshhold(/*mpz_t N, mpz_t sqrt_N, int s, int loop_number, */int* pb, int pb_len){
    
    /*
    mpz_t qstart;
    mpz_init_set_ui(qstart, s);
    mpz_mul_ui(qstart, qstart, loop_number);
    mpz_add(qstart, qstart, sqrt_N);
    mpz_mul(qstart, qstart, qstart);
    mpz_sub(qstart, qstart, N);
    */

    int t = /*mpz_sizeinbase(qstart, 2) - */(int) log2(pb[pb_len-1]) - 1;
    //printf("%d, %d\n", t, (int) log2(pb[pb_len-1]));
    return t;
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
    for(int i = 0; i<pb_len+extra; i++){
        v[i] = malloc(pb_len*sizeof(int*));
    }
    float* sinterval = malloc(s*sizeof(float));
    float* plogs = prime_logs(pb, pb_len);

    /*
    // TESTS
    mpz_t temp, p1;
    mpz_inits(temp, p1, NULL);
    // END TESTS
    */

    int* x1 = malloc(pb_len*sizeof(int));
    int* x2 = malloc(pb_len*sizeof(int));
    for(int i = 1; i < pb_len; i++){
            mpz_add_ui(zi, zi, 1);
            TonelliShanks(N, pb, x1, x2, i);

            /*
            // TESTS
            mpz_set_ui(p1, pb[i]);
            //printf("x1=%d, x2=%d, p=%d\n", x1, x2, pb[i]);

            mpz_ui_pow_ui(temp, x1[i], 2);
            assert(mpz_congruent_p(temp, N, p1));

            mpz_ui_pow_ui(temp, x2[i], 2);
            assert(mpz_congruent_p(temp, N, p1));
            // END TESTS
            */
    }

    int loop_number = 0;
    int relations_found = 0;
    int tries = 0;
    while(relations_found < pb_len + extra){
        
        for(int i = 0; i<s; i++){
            sinterval[i] = 0;
        }

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
            x1[i] = x1[i] % s;
            x2[i] = x2[i] % s;
        }

        int t = calculate_threshhold(/*N, sqrt_N, s, loop_number, */pb, pb_len);
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

                    if(!tests){
                        printf("\r");
                        //printf("sinterval = %f\n", sinterval[i]);
                        printf("%.1f%% | %f%%", (float)relations_found/(pb_len+extra-1)*100, (float)relations_found/tries*100);
                        fflush(stdout);
                    }
                }
            }
        }

        //printf("found = %d\n", relations_found);
        loop_number++;
    }

    if(!tests) printf("\n");

    mpz_clears(sqrt_N, zi, qx, NULL);
    free(x1);
    free(x2);
    free(sinterval);
    free(plogs);

    return v;
}