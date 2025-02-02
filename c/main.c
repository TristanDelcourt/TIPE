#include <stdbool.h>
#include <gmp.h>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "system.h"
#include "vector.h"
#include "parse_input.h"
#include "factorbase.h"
#include "list_matrix_utils.h"

// Include algorithms
// Dixon's method
#include "./dixon/dixon.h"

// The Quadratic Sieve
#include "./qsieve/qsieve.h"

// Multipolynomial Quadratic Sieve
#include "./mpqs/polynomial.h"
#include "./mpqs/mpqs.h"
#include "./mpqs/parallel_mpqs.h"


/**
 * 
 * 
 * START OF ALGORITHM
 * 
 */


void rebuild_mpqs(mpz_t prod, mpz_t* d, int* v, int* primes, int n1, system_t s){
    mpz_set_ui(prod, 1);
    mpz_t temp;
    mpz_init(temp);
    for(int i = 0; i<n1; i++){
        if(s->sol[i]){
            mpz_mul(prod, prod, d[s->perm[i]]);
        }
        mpz_ui_pow_ui(temp, primes[i], v[i]);
        mpz_mul(prod, prod, temp);
    }
    mpz_clear(temp);
}

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
        case PMPQS:
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
    mpz_t* d;
    struct timeval t1, t2;
    gettimeofday(&t1, 0);
    switch(input->algorithm){
        case DIXON:
            v = dixon(z, input->N, pb_len, pb, input->extra, input->quiet);
            break;
        case QSIEVE:
            v = qsieve(z, input->N, pb_len, pb, input->extra, input->sieving_interval, input->quiet);
            break;
        case MPQS:
            d = malloc(target_nb*sizeof(mpz_t));
            for(int i = 0; i < target_nb; i++){
                mpz_init(d[i]);
            }
            v = mpqs(z, d, input->N, pb_len, pb, input->extra, input->sieving_interval, input->quiet);
            break;
        case PMPQS:
            d = malloc(target_nb*sizeof(mpz_t));
            for(int i = 0; i < target_nb; i++){
                mpz_init(d[i]);
            }
            v = parallel_mpqs(z, d, input->N, pb_len, pb, input->extra, input->sieving_interval, input->quiet);
            break;
    }

    gettimeofday(&t2, 0);
    long seconds = t2.tv_sec - t1.tv_sec;
    long microseconds = t2.tv_usec - t1.tv_usec;
    double time_spent = seconds + microseconds*1e-6;
    if(!input->quiet) printf("Time to get zi: %fs\n", time_spent);
    
    mpz_t f, Z1, Z2, test1, test2;
    mpz_inits(f, Z1, Z2, test1, test2, NULL);
    
    //gaussian init
    system_t s;
    int* sum;
    switch(input->algorithm){
        case DIXON:
            s = init_gauss(v, target_nb, pb_len);
            sum = malloc(pb_len*sizeof(int));
            break;
        case QSIEVE:
            s = init_gauss(v, target_nb, pb_len);
            sum = malloc(pb_len*sizeof(int));
            break;
        case MPQS:
            // for -1
            s = init_gauss(v, target_nb, pb_len+1);
            sum = malloc((pb_len+1)*sizeof(int));
            break;
        case PMPQS:
            // for -1
            s = init_gauss(v, target_nb, pb_len+1);
            sum = malloc((pb_len+1)*sizeof(int));
            break;
    }
    if(!input->quiet) printf("2^%d solutions to iterate\n", s->n2 - s->arb);

    bool done = false;
    while(!done){
        gaussian_step(s);

        prod_vect(Z1, z, target_nb, s);
        sum_lignes(sum, v, s);
        div_vect(sum, 2, pb_len);
        
        switch(input->algorithm){
            case DIXON:
                rebuild(Z2, sum, pb, pb_len);
                break;
            case QSIEVE:
                rebuild(Z2, sum, pb, pb_len);
                break;
            case MPQS:
                rebuild_mpqs(Z2, d, sum, pb, pb_len, s);
                break;
            case PMPQS:
                rebuild_mpqs(Z2, d, sum, pb, pb_len, s);
                break;
        }

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
    switch(input->algorithm){
        case DIXON:
            break;
        case QSIEVE:
            break;
        case MPQS:
            for(int i = 0; i < target_nb; i++) mpz_clear(d[i]);
            free(d);
            break;
        case PMPQS:
            for(int i = 0; i < target_nb; i++) mpz_clear(d[i]);
            free(d);
            break;
    }
    

    mpz_clears(f, Z1, Z2, test1, test2, NULL);
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

    struct timeval t1, t2;
    gettimeofday(&t1, 0);
    factor(input);
    gettimeofday(&t2, 0);
    long seconds = t2.tv_sec - t1.tv_sec;
    long microseconds = t2.tv_usec - t1.tv_usec;
    double time_spent = seconds + microseconds*1e-6;
    if(!input->quiet) printf("Total time: %fs\n", time_spent);

    free_input(input);

    return 0;
}