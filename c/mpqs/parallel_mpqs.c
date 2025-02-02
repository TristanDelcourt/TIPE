#include <gmp.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

#include "polynomial.h"
#include "common_mpqs.h"
#include "parallel_mpqs.h"
#include "../system.h"
#include "../tonellishanks.h"

pthread_mutex_t mutex;


void* sieve_100_polys (void* args){
    sieve_arg_t* arg = (sieve_arg_t*) args;

    poly_t Q = copy_poly(arg->Qinit);

    mpz_t temp, g, m, n, pi, x;
    mpz_inits(temp, g, m, n, pi, x, NULL);
    float* sinterval = malloc(2*arg->s*sizeof(float));
    int* x1 = malloc(arg->pb_len*sizeof(int));
    int* x2 = malloc(arg->pb_len*sizeof(int));

    for(int i = 0; i<100 && *(arg->relations_found) < arg->pb_len + arg->extra; i++){
        get_next_poly(Q);

        //get sol for 2
        mpz_set_ui(temp, 0);
        calc_poly(Q, temp);
        x1[0] = 0;
        if(mpz_divisible_ui_p(Q->qx, 2) == 0) x1[0] = 1;
        
        //get sol for others
        for(int i = 1; i<arg->pb_len; i++){
            mpz_set_ui(pi, arg->pb[i]);
            mpz_gcdext(g, m, n, Q->a, pi);
            assert(mpz_cmp_ui(g, 1) == 0);
            mpz_set_ui(temp, arg->r[i]);
            mpz_sub(temp, temp, Q->b);
            mpz_mul(temp, temp, m);
            mpz_mod(temp, temp, pi);
            x1[i] = mpz_get_ui(temp);
            calc_poly(Q, temp);
            assert(mpz_divisible_ui_p(Q->qx, arg->pb[i]) != 0);
            mpz_set_ui(temp, arg->pb[i]);
            mpz_sub_ui(temp, temp, arg->r[i]);
            mpz_sub(temp, temp, Q->b);
            mpz_mul(temp, temp, m);
            mpz_mod(temp, temp, pi);
            x2[i] = mpz_get_ui(temp);
            calc_poly(Q, temp);
            assert(mpz_divisible_ui_p(Q->qx, arg->pb[i]) != 0);

            //realign sieving interval to [-s, s]
            int k = (x1[i] + arg->s)/arg->pb[i];
            x1[i] -= k * arg->pb[i];
            x1[i] += arg->s;
            k = (x2[i] + arg->s)/arg->pb[i];
            x2[i] -= k * arg->pb[i];
            x2[i] += arg->s;
            mpz_set_si(temp, -arg->s);
            mpz_add_ui(temp, temp, x1[i]);
            calc_poly(Q, temp);
            assert(mpz_divisible_ui_p(Q->qx, arg->pb[i]) != 0);
        }

        //reset sieveing_interval
        for(int i = 0; i<2*arg->s; i++){
            sinterval[i] = 0;
        }

        // sieve for 2
        while(x1[0]<2*arg->s){
            sinterval[x1[0]] += arg->plogs[0];
            x1[0] += arg->pb[0];
        }
        // sieve other primes
        for(int i = 1; i < arg->pb_len; i++){
            while(x1[i]<2*arg->s){
                sinterval[x1[i]] += arg->plogs[i];
                x1[i] += arg->pb[i];
            }
            while(x2[i]<2*arg->s){
                sinterval[x2[i]] += arg->plogs[i];
                x2[i] += arg->pb[i];
            }
        }

        bool found;
        bool update_time = false;
        pthread_mutex_lock(&mutex);
        for(int i = 0; i<2*arg->s && *(arg->relations_found) < arg->pb_len + arg->extra; i++){
            if(sinterval[i] > arg->t){
                mpz_set_si(x, -arg->s);
                mpz_add_ui(x, x, i);
                calc_poly(Q, x);

                if(!already_added(Q->zi, arg->z, *(arg->relations_found))){
                    found = vectorize_mpqs(Q->qx, arg->v[*(arg->relations_found)], arg->pb_len, arg->pb);
                    if(found){
                        mpz_set(arg->z[*(arg->relations_found)], Q->zi);
                        mpz_set(arg->d[*(arg->relations_found)], Q->d);
                        *(arg->relations_found) += 1;
                        found = false;
                        update_time = true;
                        if(!arg->quiet){
                            printf("\r");
                            printf("%.1f%%", (float)(*(arg->relations_found))/(arg->pb_len+arg->extra)*100);
                            fflush(stdout);
                        }
                    }
                }
            }
        }
        
        struct timeval current;
        gettimeofday(&current, 0);
        long seconds = current.tv_sec - arg->begin.tv_sec;
        long microseconds = current.tv_usec - arg->begin.tv_usec;
        double elapsed = seconds + microseconds*1e-6;
        if(update_time && !arg->quiet) printf(" (~%.0fs left)        " , elapsed/(*arg->relations_found)*(arg->pb_len+arg->extra - (*arg->relations_found)));
        pthread_mutex_unlock(&mutex);
    }

    mpz_clears(temp, g, m, n, pi, x, NULL);
    free(x1);
    free(x2);
    free(sinterval);
    free_poly(Q);

    return NULL;
}

int** parallel_mpqs(mpz_t* z, mpz_t* d, mpz_t N, int pb_len, int* pb, int extra, int s, bool quiet){
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

    poly_t Q = init_poly(N, s);

    int** v = malloc((pb_len+extra)*sizeof(int*));
    for(int i = 0; i<pb_len+extra; i++){
        v[i] = malloc((pb_len+1)*sizeof(int*)); // +1 for -1
    }
    float* plogs = prime_logs_mpqs(pb, pb_len);
    

    int* r = malloc(pb_len*sizeof(int));
    int sol1, sol2;
    for(int i = 1; i < pb_len; i++){
        tonelli_shanks_ui(N, pb[i], &sol1, &sol2);
        r[i] = sol1;
    }
    int t = calculate_threshhold_mpqs(sqrt_N, s, pb, pb_len);

    sieve_arg_t* args = malloc(8*sizeof(sieve_arg_t));
    pthread_t* threads = malloc(8*sizeof(pthread_t));
    
    int relations_found = 0;
    struct timeval begin;
    gettimeofday(&begin, 0);
    while(relations_found < pb_len + extra){
        for(int i = 0; i<8; i++){
            args[i] = (sieve_arg_t) {
                pb,
                pb_len,
                extra,
                r,
                plogs,
                s,
                t,
                &relations_found,
                v,
                quiet,
                z,
                d,
                Q,
                begin
            };
            pthread_create(threads+i, NULL, sieve_100_polys, args+i);

            for(int i = 0; i<100; i++){
                get_next_poly(Q);
            }
        }

        for(int i = 0; i<8; i++){
            pthread_join(threads[i], NULL);
        }
    }
    if(!quiet) printf("\n");

    free(threads);
    free(args);
    free(r);
    free(plogs);
    free_poly(Q);
    mpz_clear(sqrt_N);

    return v;
}