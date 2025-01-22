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

void TonelliShanks(mpz_t a, int p, int* x1, int* x2){
    //Solve the equation: x^2 ≡ a (mod p)
        
        
    //1. Factor `p - 1` into `2^S * Q` where Q is odd.
    int Q = p - 1;
    int S = 0;
    while(Q % 2 == 0){
        S += 1;
        Q /= 2;
    }

    // 2. Find a NR(p).
    mpz_t y, p1;
    mpz_init_set_ui(y, 2);
    mpz_init_set_ui(p1, p);

    while(mpz_legendre(y, p1) != -1){
        mpz_add_ui(y, y, 1);
    }

    // 3. Calculate the four quantities.
    mpz_t R, c, t, E;
    mpz_inits(R, c, t, E, NULL);

    mpz_powm_ui(R, a, (Q+1)/2, p1);
    mpz_powm_ui(c, y, Q, p1);
    mpz_powm_ui(t, a, Q, p1);
    mpz_set_ui(E, S);

    mpz_t temp, b;
    mpz_inits(temp, b, NULL);

    while(mpz_cmp_ui(t, 1) != 0){
        int i = 1;
        while(mpz_cmp_ui(E, i) != 0){
            mpz_set_ui(temp, 1);
            mpz_mul_2exp(temp, temp, i);
            mpz_powm(temp, t, temp, p1);
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
        mpz_powm(b, c, temp, p1);

        // R = R * b % p
        mpz_mul(R, R, b);
        mpz_mod_ui(R, R, p);

        // c = pow(b, 2, p)
        mpz_pow_ui(c, b, 2);
        mpz_mod_ui(c, c, p);

        // t = c * t % p
        mpz_mul(t, c, t);
        mpz_mod_ui(t, t, p);

        mpz_set_ui(E, i);
    }

    *x1 = mpz_get_ui(R);
    *x2 = p - mpz_get_ui(R);

    mpz_clears(y, p1, R, c, t, E, temp, b, NULL);
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
    float* plogs = prime_logs(pb, pb_len);

    int  x1, x2;

    // TESTS
    mpz_t temp, p1;
    mpz_inits(temp, p1, NULL);
    // END TESTS

    for(int i = 1; i < pb_len; i++){
        mpz_add_ui(zi, zi, 1);
        
        TonelliShanks(N, pb[i], &x1, &x2);
        
        // TESTS
        mpz_set_ui(p1, pb[i]);
        printf("x1=%d, x2=%d, p=%d\n", x1, x2, pb[i]);

        mpz_ui_pow_ui(temp, x1, 2);
        assert(mpz_congruent_p(temp, N, p1));

        mpz_ui_pow_ui(temp, x2, 2);
        assert(mpz_congruent_p(temp, N, p1));
        // END TESTS
        
        while(x1<s){
            sinterval[x1] += plogs[i];
            x1 += pb[i];
        }

        while(x2<s){
            sinterval[x1] += plogs[i];
            x2 += pb[i];
        }

        for(int i = 0; i<s; i++){
            if(sinterval[i] > t){
                
            }
        }
        
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