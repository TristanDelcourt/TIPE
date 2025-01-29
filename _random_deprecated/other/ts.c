#include <gmp.h>

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
    x2[j] = p[j] - x1[j];

    mpz_clears(y, pj, R, c, t, E, temp, b, NULL);
}

float* prime_logs(int* pb, int pb_len){
    float* plogs = malloc(pb_len*sizeof(float));
    
    for(int i = 0; i<pb_len; i++){
        plogs[i] = log2(pb[i]);
    }

    return plogs;
}
