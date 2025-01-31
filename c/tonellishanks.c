#include <stdint.h>
#include <gmp.h>


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

void tonelli_shanks_ui(mpz_t n, int p, int* x1, int* x2) {
    uint64_t q = p - 1;
    uint64_t ss = 0;
    uint64_t z = 2;
    uint64_t c, r, t, m;

    while ((q & 1) == 0) {
        ss += 1;
        q >>= 1;
    }

    mpz_t temp, pj;
    mpz_init(temp);
    mpz_init_set_ui(pj, p);

    if (ss == 1) {
        //uint64_t r1 = modpow(n, (p + 1) / 4, p);
        mpz_powm_ui(temp, n, (p+1)/4, pj);
        uint64_t r1 = mpz_get_ui(temp);

        *x1 = r1;
        *x2 = p - r1;
        return;
    }

    while (modpow(z, (p - 1) / 2, p) != (uint64_t) p - 1) { // uint_64 only there for the compiler to stop complaining
        z++;
    }

    c = modpow(z, q, p);
    
    //r = modpow(n, (q + 1) / 2, p);
    mpz_powm_ui(temp, n, (q+1)/2, pj);
    r = mpz_get_ui(temp);

    //t = modpow(n, q, p);
    mpz_powm_ui(temp, n, q, pj);
    t = mpz_get_ui(temp);
    
    m = ss;

    while(1){
        uint64_t i = 0, zz = t;
        uint64_t b = c, e;
        if (t == 1) {
            *x1 = r;
            *x2 = p - r;
            return;
        }
        while (zz != 1 && i < (m - 1)) {
            zz = zz * zz % p;
            i++;
        }
        e = m - i - 1;
        while (e > 0) {
            b = b * b % p;
            e--;
        }
        r = r * b % p;
        c = b * b % p;
        t = t * c % p;
        m = i;
    }
}

#include <gmp.h>

void tonelli_shanks_mpz(mpz_t a, mpz_t p, mpz_t x1, mpz_t x2){
    //Solve the equation: x^2 ≡ a (mod p)
        
        
    //1. Factor `p - 1` into `2^S * Q` where Q is odd.
    mpz_t Q;
    mpz_init_set(Q, p);
    mpz_sub_ui(Q, Q, 1);
    int S = 0;
    while(mpz_divisible_ui_p(Q, 2) != 0){
        S += 1;
        mpz_divexact_ui(Q, Q, 2);
    }

    // 2. Find a NR(p).
    mpz_t y;
    mpz_init_set_ui(y, 2);
    while(mpz_legendre(y, p) != -1){
        mpz_add_ui(y, y, 1);
    }

    // 3. Calculate the four quantities.
    mpz_t R, c, t, E;
    mpz_inits(R, c, t, E, NULL);

    mpz_powm(c, y, Q, p);
    mpz_powm(t, a, Q, p);
    mpz_set_ui(E, S);

    mpz_add_ui(Q, Q, 1);
    mpz_divexact_ui(Q, Q, 2);
    mpz_powm(R, a, Q, p);

    mpz_t temp, b;
    mpz_inits(temp, b, NULL);

    while(mpz_cmp_ui(t, 1) != 0){
        int i = 1;
        while(mpz_cmp_ui(E, i) != 0){
            mpz_set_ui(temp, 1);
            mpz_mul_2exp(temp, temp, i);
            mpz_powm(temp, t, temp, p);
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
        mpz_powm(b, c, temp, p);

        // R = R * b % p
        mpz_mul(R, R, b);
        mpz_mod(R, R, p);

        // c = pow(b, 2, p)
        mpz_pow_ui(c, b, 2);
        mpz_mod(c, c, p);

        // t = c * t % p
        mpz_mul(t, c, t);
        mpz_mod(t, t, p);

        mpz_set_ui(E, i);
    }

    mpz_set(x1, R);
    mpz_sub(x2, p, x1);

    mpz_clears(y, R, c, t, E, Q, temp, b, NULL);
}