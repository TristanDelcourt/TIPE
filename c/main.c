#include <stdbool.h>
#include <gmp.h>
#include <time.h>
#include <stdio.h>
#include "system.h"
#include "parse_input.h"

// Include aglorithms
#include "./qsieve/qsieve.h"


void factor(mpz_t N, int B, int extra, bool tests, TYPE algorithm){
    int piB = pi(B);
    if(!tests) printf("pi(B) = %d\n", piB);
    int* p = primes(piB, B);

    /*
    int pb_len;
    int* pb = prime_base(N, &pb_len, p, piB);
    if(!tests) printf("base reduction %f%%\n", (float)pb_len/piB*100);
    free(p);
    */
    int* pb = p;
    int pb_len = piB;

    mpz_t* z = malloc((pb_len+extra)*sizeof(mpz_t));
    for(int i = 0; i < pb_len+extra; i++){
        mpz_init(z[i]);
    }
    
    //Getting zis
    int** v;
    clock_t t1 = clock();
    switch(algorithm){
        case DIXON:
            printf(stderr, "ERROR: not done\n");
            break;
        case QSIEVE:
            v = qsieve(z, N, pb_len, pb, extra, tests);
            break;
    }
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!tests) printf("Time to get zi: %fs\n", time_spent);
    
    mpz_t f1, f2, Z1, Z2;
    mpz_inits(f1, f2, Z1, Z2, NULL);
    
    //gaussian init
    system_t s = init_gauss(v, pb_len+extra, pb_len);
    if(!tests) printf("2^%d solutions to iterate\n", s->n2 - s->arb);
    int* sum = malloc(pb_len*sizeof(int));

    bool done = false;
    while(!done){

        gaussian_step(s);

        sum_lignes(sum, v, s);
        div_vect(sum, 2, pb_len);

        prod_vect(Z1, z, pb_len+extra);
        rebuild(Z2, sum, pb, pb_len);

        mpz_sub(f1, Z1, Z2);
        mpz_add(f2, Z1, Z2);

        mpz_gcd(f1, f1, N);
        mpz_gcd(f2, f2, N);

        if((!(mpz_cmp_ui(f1, 1) == 0) && !(mpz_cmp(f1, N) == 0))
            || (!(mpz_cmp_ui(f2, 1) == 0) && (!(mpz_cmp(f2, N) == 0)))){
            done = true;
        }

        if(s->done){
            if(!tests) fprintf(stderr, "ERROR: no solution for this set of zi\n");
            exit(1);
        }
    }
    free(sum);
    free(pb);
    free_system(s);
    free_ll(v, pb_len+extra, pb_len);
    for(int i = 0; i < pb_len+extra; i++){
        mpz_clear(z[i]);
    }
    free(z);

    if(!tests) gmp_printf("%Zd = 0 [%Zd]\n%Zd = 0 [%Zd]\n", N, f1, N, f2);

    assert(mpz_divisible_p(N, f1));
    assert(mpz_divisible_p(N, f2));

    mpz_clears(f1, f2, Z1, Z2, NULL);
}

int main(int argc, char** argv){
    input_t* input = parse_input(argc, argv);
    if(input==NULL){
        printf(stderr, "ERROR: Invalid input\n");
        return 1;
    }

    if(mpz_cmp_ui(input->N, 0) == 0){
        printf(stderr, "ERROR: No input number, use -n %%number%%\n");
        return 1;
    }

    if(input->bound == -1) input->bound = 500;
    if(input->sieving_interval == -1) input->sieving_interval = 1000;
    if(input->extra == -1) input->extra = 1;

    factor(
        input->N,
        input->bound,
        input->extra,
        input->quiet,
        input->algorithm
    );

    return 0;
}