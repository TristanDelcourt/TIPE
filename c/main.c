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

#include "./mpqs/polynomial.h"
#include "./mpqs/mpqs.h"


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

void factor(input_t* input){

    int piB = pi(input->bound);
    if(!input->quiet) printf("pi(B) = %d\n", piB);
    int* p = primes(piB, input->bound);

    int pb_len;
    int* pb;
    switch(input->algorithm){
        case DIXON:
            pb = p;
            pb_len = piB;
            break;
        case QSIEVE:
            pb = prime_base(input->N, &pb_len, p, piB);
            if(!input->quiet) printf("base reduction %f%%\n", (float)pb_len/piB*100);
            free(p);
            break;
        case MPQS:
            pb = prime_base(input->N, &pb_len, p, piB);
            pb[pb_len] = -1;
            if(!input->quiet) printf("base reduction %f%%\n", (float)pb_len/piB*100);
            free(p);
            break;
    }
    int target_nb = pb_len + input->extra;

    mpz_t* z = malloc((target_nb)*sizeof(mpz_t));
    for(int i = 0; i < target_nb; i++){
        mpz_init(z[i]);
    }
    
    //Getting zis
    int** v;
    clock_t t1 = clock();
    switch(input->algorithm){
        case DIXON:
            v = dixon(z, input->N, pb_len, pb, input->extra, input->quiet);
            break;
        case QSIEVE:
            v = qsieve(z, input->N, pb_len, pb, input->extra, input->sieving_interval, input->quiet);
            break;
        case MPQS:
            v = mpqs(z, input->N, pb_len, pb, input->extra, input->sieving_interval, input->quiet);
            break;
    }
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!input->quiet) printf("Time to get zi: %fs\n", time_spent);
    
    mpz_t f, Z1, Z2, test1, test2;
    mpz_inits(f, Z1, Z2, test1, test2, NULL);
    
    //gaussian init
    system_t s = init_gauss(v, target_nb, pb_len);
    if(!input->quiet) printf("2^%d solutions to iterate\n", s->n2 - s->arb);
    int* sum = malloc(pb_len*sizeof(int));

    bool done = false;
    while(!done){
        gaussian_step(s);

        prod_vect(Z1, z, target_nb, s);
        sum_lignes(sum, v, s);
        div_vect(sum, 2, pb_len);
        rebuild(Z2, sum, pb, pb_len);

        mpz_set(test1, Z1);
        mpz_mul(test1, test1, test1);
        mpz_set(test2, Z2);
        mpz_mul(test2, test2, test2);
        assert(mpz_congruent_p(test1, test2, input->N) != 0);

        mpz_sub(f, Z1, Z2);
        mpz_gcd(f, f, input->N);

        if(mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, input->N) != 0){
            assert(mpz_divisible_p(input->N, f));
            if(!input->quiet) gmp_printf("%Zd = 0 [%Zd]\n", input->N, f);
            done = true;
        }
        
        mpz_add(f, Z1, Z2);
        mpz_gcd(f, f, input->N);

        if(mpz_cmp_ui(f, 1) != 0 && mpz_cmp(f, input->N) != 0){
            assert(mpz_divisible_p(input->N, f));
            if(!input->quiet) gmp_printf("%Zd = 0 [%Zd]\n", input->N, f);
            done = true;
        }

        if(s->done){
            if(!input->quiet) fprintf(stderr, "ERROR: no solution for this set of zi\n");
            exit(1);
        }
    }

    free(sum);
    free(pb);
    free_system(s);
    free_ll(v, target_nb);
    for(int i = 0; i < target_nb; i++){
        mpz_clear(z[i]);
    }
    free(z);



    mpz_clears(f, Z1, Z2, NULL);
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

    if(input->bound == -1) input->bound = 1000;
    if(input->sieving_interval == -1) input->sieving_interval = 100000;
    if(input->extra == -1) input->extra = 1;

    poly_t q = init_poly(input->N, input->sieving_interval);

    clock_t t1 = clock();
    factor(input);
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!input->quiet) printf("Total time: %fs\n", time_spent);

    free_input(input);

    return 0;
}