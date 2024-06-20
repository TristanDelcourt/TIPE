#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void print_list(int* l, int n){
    for(int i = 0; i<n; i++){
        printf("%d ", l[i]);
    }
}

bool factorise(mpz_t n, int* vi, int piB){
    for(int i = 0; i<piB; i++){
        vi[i] = 0;
    }

    for(int d = 2; d < 4; d++) {
        while (mpz_congruent_ui_p(n, 0, d)){
            mpz_divexact_ui(n, n, d);
            vi[d-2]++;
        }

    }
    int d = 5;
    int add = 2;
    int i = 2;
    while (mpz_cmp_ui(n, 1)) {
        if(i>=piB){
            return false;
        }

        while (mpz_congruent_ui_p(n, 0, d)){
            mpz_divexact_ui(n, n, d);
            vi[i]++;
        }
        d += add;
        add = 6 - add;
        i++;
    }

    return true;
}

int* get_zi(mpz_t zi, mpz_t N, int piB, gmp_randstate_t rstate){
    bool found = false;
    int* vi = malloc(piB*sizeof(int));

    while(!found){
        mpz_urandomm(zi, rstate, N);
        mpz_t zi_2;
        mpz_init(zi_2);
        mpz_mul(zi_2, zi, zi);
        mpz_mod(zi_2, zi_2, N);
        //gmp_printf("%Zd\n", zi_2);
        found = factorise(zi_2, vi, piB);
    }

    return vi;
}


int** get_all_zi(mpz_t* z, mpz_t N, int piB){
    gmp_randstate_t rstate;
    //unsigned long seed;
    gmp_randinit_default(rstate);
    //gmp_randseed_ui(rstate, seed);

    int** v = malloc((piB+1)*sizeof(int*));
    for(int i = 0; i < piB+1; i++){
        v[i] = get_zi(z[i], N, piB, rstate);  
    }

    return v;
}

void dixon(mpz_t N, int piB){

    mpz_t* z = malloc((piB+1)*sizeof(mpz_t));
    for(int i = 0; i < piB +1; i++){
        mpz_init(z[i]);
    }

    int** v = get_all_zi(z, N, piB);
    for(int i = 0; i < piB+1; i++){
        gmp_printf("%zd\n", z[i]);
        print_list(v[i], piB);
        printf("\n");
    }

}

void main(int agrc, char**argv){
    int B = 10000;
    int piB = 168;
    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    dixon(N, piB);
    

}