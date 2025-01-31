#include <stdint.h>
#include <gmp.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>

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

void tonelli_shanks_ui(mpz_t n, unsigned long int p, int* x1, int* x2) {
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

    while (modpow(z, (p - 1) / 2, p) != (unsigned long int) p - 1) { // uint_64 only there for the compiler to stop complaining
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

void tonelli_shanks_mpz(mpz_t n, mpz_t p, mpz_t x1, mpz_t x2){
    mpz_t q, z;
    mpz_init_set(q, p);
    mpz_sub_ui(q, q, 1);
    int ss = 0;
    mpz_init_set_ui(z, 2);

    while(mpz_divisible_ui_p(q, 2) != 0){
        ss += 1;
        mpz_divexact_ui(q, q, 2);
    }

    mpz_t op1;
    mpz_init(op1);

    if (ss == 1) {
        //uint64_t r1 = modpow(n, (p + 1) / 4, p);
        mpz_add_ui(op1, p, 1);
        mpz_divexact_ui(op1, op1, 4);
        mpz_powm(op1, n, op1, p);

        mpz_set(x1, op1);
        mpz_sub(x2, p, x1);

        mpz_clears(q, z, op1, NULL);
        return;
    }

    mpz_t op2, op3;
    mpz_inits(op2, op3, NULL);

    mpz_sub_ui(op1, p, 1);
    mpz_divexact_ui(op1, op1, 2);
    mpz_powm(op2, z, op1, p);

    mpz_sub_ui(op3, p, 1);
    while(mpz_cmp(op3, op2) != 0){
        mpz_add_ui(z, z, 1);
        mpz_powm(op3, z, op1, p);
    }

    mpz_t c, r, t, m, i, zz, b, e;
    mpz_inits(c, r, t, m, i, zz, b, e, NULL);
    mpz_powm(c, z, q, p);

    mpz_add_ui(op1, q, 1);
    mpz_divexact_ui(op1, op1, 2);
    mpz_powm(r, n, op1, p);

    mpz_powm(t, n, q, p);

    mpz_set_ui(m, ss);

    while(1){
        

        mpz_set_ui(i, 0);
        mpz_set(zz, t);
        mpz_set(b, c);

        if(mpz_cmp_ui(t, 1) == 0){
            mpz_set(x1, r);
            mpz_sub(x2, p, x1);

            mpz_clears(c, r, t, m, i, zz, b, e, op1, op2, op3, q, z, NULL);
            return;
        }

        mpz_sub_ui(op1, m, 1);
        while(mpz_cmp_ui(zz, 1) != 0 && mpz_cmp(i, op1)<0){
            mpz_mul(zz, zz, zz);
            mpz_mod(zz, zz, p);
            mpz_add_ui(i, i, 1);
        }

        mpz_sub(e, m, i);
        mpz_sub_ui(e, e, 1);
        while(mpz_cmp_ui(e, 0)>0){
            mpz_mul(b, b, b);
            mpz_mod(b, b, p);
            mpz_sub_ui(e, e, 1);
        }

        mpz_mul(r, r, b);
        mpz_mod(r, r, p);

        mpz_mul(c, b, b);
        mpz_mod(c, c, p);

        mpz_mul(t, t, c);
        mpz_mod(t, t, p);

        mpz_set(m, i);
    }

}