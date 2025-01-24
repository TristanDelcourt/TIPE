#include <stdbool.h>
#include <gmp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "system.h"
#include "vector.h"
#include "parse_input.h"
#include "factorbase.h"
#include "list_matrix_utils.h"

// Include algorithms
#include "./dixon/dixon.h"
#include "./qsieve/qsieve.h"

/**
 * 
 * 
 * START OF ALGORITHM
 * 
 */


void rebuild(mpz_t prod, int* v, int* primes, int n1){
    /** Rebuilds the product of primes to the power of half
     * the solution found by the gaussian solve

     * EX:
     * v = (1, 2, 3, 1)
     * primes = [2, 3, 5, 7]
     * prod = 2**1 * 3** 2 * 5**3 * 7**1
     * returns prod
     * 
    */

    mpz_set_ui(prod, 1);
    mpz_t temp;
    mpz_init(temp);
    for(int i = 0; i<n1; i++){
        mpz_ui_pow_ui(temp, primes[i], v[i]);
        mpz_mul(prod, prod, temp);
    }
    mpz_clear(temp);
}

void sum_lignes(int* sum, int** v, system_t s){
    /** Sums the lines of vectors into 'sum' according the solution of the
     * output of the system 's', such that each power is even
     */
    for(int i = 0; i<s->n1; i++){
        sum[i] = 0;
    }

    for(int i = 0; i<s->n2; i++){
        if(s->sol[i]){
            add_vect(sum, v[s->perm[i]], s->n1);
        }
    }
}

void factor(mpz_t N, int B, int extra, bool quiet, int sinterval, TYPE algorithm){
    int piB = pi(B);
    if(!quiet) printf("pi(B) = %d\n", piB);
    int* p = primes(piB, B);

    int pb_len;
    int* pb;
    switch(algorithm){
        case DIXON:
            pb = p;
            pb_len = piB;
            break;
        case QSIEVE:
            pb = prime_base(N, &pb_len, p, piB);
            if(!quiet) printf("base reduction %f%%\n", (float)pb_len/piB*100);
            free(p);
            break;
        case DEFAULT:
            pb = prime_base(N, &pb_len, p, piB);
            if(!quiet) printf("base reduction %f%%\n", (float)pb_len/piB*100);
            free(p);
            break;
    }

    mpz_t* z = malloc((pb_len+extra)*sizeof(mpz_t));
    for(int i = 0; i < pb_len+extra; i++){
        mpz_init(z[i]);
    }
    
    //Getting zis
    int** v;
    clock_t t1 = clock();
    switch(algorithm){
        case DIXON:
            v = dixon(z, N, pb_len, pb, extra, quiet);
            break;
        case QSIEVE:
            v = qsieve(z, N, pb_len, pb, extra, sinterval, quiet);
            break;
        case DEFAULT:
            v = qsieve(z, N, pb_len, pb, extra, sinterval, quiet);
            break;
    }
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!quiet) printf("Time to get zi: %fs\n", time_spent);
    
    mpz_t f1, f2, Z1, Z2;
    mpz_inits(f1, f2, Z1, Z2, NULL);
    
    //gaussian init
    system_t s = init_gauss(v, pb_len+extra, pb_len);
    if(!quiet) printf("2^%d solutions to iterate\n", s->n2 - s->arb);
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
            if(!quiet) fprintf(stderr, "ERROR: no solution for this set of zi\n");
            exit(1);
        }
    }
    free(sum);
    free(pb);
    free_system(s);
    free_ll(v, pb_len+extra);
    for(int i = 0; i < pb_len+extra; i++){
        mpz_clear(z[i]);
    }
    free(z);

    if(!quiet) gmp_printf("%Zd = 0 [%Zd]\n%Zd = 0 [%Zd]\n", N, f1, N, f2);

    assert(mpz_divisible_p(N, f1));
    assert(mpz_divisible_p(N, f2));

    mpz_clears(f1, f2, Z1, Z2, NULL);
}

int main(int argc, char** argv){
    input_t* input = parse_input(argc, argv);
    if(input==NULL){
        fprintf(stderr, "ERROR: Invalid input\n");
        return 1;
    }

    if(mpz_cmp_ui(input->N, 0) == 0){
        fprintf(stderr, "ERROR: No input number, use -n %%number%%\n");
        return 1;
    }

    if(input->bound == -1) input->bound = 500;
    if(input->sieving_interval == -1) input->sieving_interval = 1000;
    if(input->extra == -1) input->extra = 1;

    clock_t t1 = clock();
    factor(
        input->N,
        input->bound,
        input->extra,
        input->quiet,
        input->sieving_interval,
        input->algorithm
    );
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!input->quiet) printf("Total time: %fs\n", time_spent);

    return 0;
}