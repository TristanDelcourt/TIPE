#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
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


uint64_t modpow(uint64_t a, uint64_t b, uint64_t n) {
    uint64_t x = 1, y = a;
    while (b > 0) {
        if (b % 2 == 1) {
            x = (x * y) % n; // multiplying with base
        }
        y = (y * y) % n; // squaring the base
        b /= 2;
    }
    return x % n;
}

void ts(mpz_t n, int* p, int* x1, int*x2, int j) {
    uint64_t q = p[j] - 1;
    uint64_t ss = 0;
    uint64_t z = 2;
    uint64_t c, r, t, m;

    while ((q & 1) == 0) {
        ss += 1;
        q >>= 1;
    }

    mpz_t temp, pj;
    mpz_init(temp);
    mpz_init_set_ui(pj, p[j]);

    if (ss == 1) {
        //uint64_t r1 = modpow(n, (p + 1) / 4, p);
        mpz_powm_ui(temp, n, (p[j]+1)/4, pj);
        uint64_t r1 = mpz_get_ui(temp);

        x1[j] = r1;
        x2[j] = p[j] - r1;
        return;
    }

    while (modpow(z, (p[j] - 1) / 2, p[j]) != (uint64_t) p[j] - 1) { // uint_64 only there for the compiler to stop complaining
        z++;
    }

    c = modpow(z, q, p[j]);
    
    //r = modpow(n, (q + 1) / 2, p);
    mpz_powm_ui(temp, n, (q+1)/2, pj);
    r = mpz_get_ui(temp);

    //t = modpow(n, q, p);
    mpz_powm_ui(temp, n, q, pj);
    t = mpz_get_ui(temp);
    
    m = ss;

    while (true) {
        uint64_t i = 0, zz = t;
        uint64_t b = c, e;
        if (t == 1) {
            x1[j] = r;
            x2[j] = p[j] - r;
            return;
        }
        while (zz != 1 && i < (m - 1)) {
            zz = zz * zz % p[j];
            i++;
        }
        e = m - i - 1;
        while (e > 0) {
            b = b * b % p[j];
            e--;
        }
        r = r * b % p[j];
        c = b * b % p[j];
        t = t * c % p[j];
        m = i;
    }
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

int** qsieve(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra, int s, bool tests){
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


    for(int i = 1; i < pb_len; i++){
            ts(N, pb, x1, x2, i);

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
                
                //printf("sinterval = %f\n", sinterval[i]);
                
                if(found){
                    mpz_set(z[relations_found], zi);
                    relations_found++;
                    found = false;

                    if(!tests){
                        printf("\r");
                        printf("%.1f%% | %.1f%%", (float)relations_found/(pb_len+extra-1)*100, (float)relations_found/tries*100);
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