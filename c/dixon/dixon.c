#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <time.h>
#include "system.h"
#include "vector.h"
#include "list_matrix_utils.h"

//#define MUL_KARATSUBA_THRESHHOLD 10

bool is_prime(int n) { 
    // Corner cases 
    if (n <= 1) 
        return false; 
    if (n <= 3) 
        return true; 
  
    // This is checked so that we can skip 
    // middle five numbers in below loop 
    if (n % 2 == 0 || n % 3 == 0) 
        return false; 
  
    for (int i = 5; i * i <= n; i = i + 6) 
        if (n % i == 0 || n % (i + 2) == 0) 
            return false; 
  
    return true; 
} 
  
int pi(int n) { 
    int k = 0;
    for (int i = 2; i <= n; i++) { 
        if (is_prime(i)) 
            k++;; 
    } 
    return k;
} 

int* primes(int piB, int B){
    int* p = malloc(piB*sizeof(int));
    int k = 0;
    for (int i = 2; i <= B; i++) { 
        if (is_prime(i)){
            p[k] = i;
            k++;
        }
    }
    return p;

}

/* Used for legendre symbol, exists in gmp already
bool euler_criterion(mpz_t n, int p){
    int e = (p-1)/2;
    mpz_t r, p1;
    mpz_init(r);
    mpz_init_set_ui(p1, p);
    mpz_powm_ui(r, n, e, p1);
    return(mpz_cmp_ui(r, 1) == 0);
}
*/

int* prime_base(mpz_t n, int* pb_len, int* primes, int piB){
    int* pb = malloc(piB*sizeof(int));
    pb[0] = 2;

    int j = 1;
    mpz_t p1;
    mpz_init(p1);
    for(int i = 1; i<piB; i++){
        mpz_set_ui(p1, primes[i]);
        if(mpz_legendre(n, p1) == 1){
            //printf("%d\n", primes[i]);
            pb[j] = primes[i];
            j++;
        }
    }
    *pb_len = j;
    pb = realloc(pb, j*sizeof(int));
    
    mpz_clear(p1);
    return pb;
}


void rebuild(mpz_t prod, int* v, int* primes, int n1){
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
    for(int i = 0; i<s->n1; i++){
        sum[i] = 0;
    }

    for(int i = 0; i<s->n2; i++){
        if(s->sol[i]){
            add_vect(sum, v[s->perm[i]], s->n1);
        }
    }
}


bool factorise(mpz_t n, int* v, int pb_len, int* pb){
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

int** get_all_zi(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra, bool tests){
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

    for(int i = 0; i < pb_len+extra; i++){
        bool found = false;
        int* vi = malloc(pb_len*sizeof(int));

        while(!found){
            mpz_add_ui(zi, zi, 1);

            //Q(x)
            mpz_mul(qx, zi, zi);
            mpz_sub(qx, qx, N);

            found = factorise(qx, vi, pb_len, pb);
        }
        if(!tests){
            printf("\r");
            printf("%f%%", (float)i/(pb_len+extra-1)*100);
            fflush(stdout);
        }
        
        v[i] = vi;
        mpz_set(z[i], zi);
    }
    if(!tests) printf("\n");

    mpz_clears(sqrt_N, zi, qx, NULL);


    return v;
}


void dixon(mpz_t N, int B, int extra, bool tests){
    int piB = pi(B);
    if(!tests) printf("pi(B) = %d\n", piB);
    int* p = primes(piB, B);

    int pb_len;
    int* pb = prime_base(N, &pb_len, p, piB);
    if(!tests) printf("base reduction %f%%\n", (float)pb_len/piB*100);
    free(p);
    /*
    int* pb = p;
    int pb_len = piB;
    */

    mpz_t* z = malloc((pb_len+extra)*sizeof(mpz_t));
    for(int i = 0; i < pb_len+extra; i++){
        mpz_init(z[i]);
    }
    
    //Getting zis
    clock_t t1 = clock();
    int** v = get_all_zi(z, N, pb_len, pb, extra, tests);
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!tests) printf("Time to get zi: %fs\n", time_spent);
    
    mpz_t f1, f2, Z1, Z2;
    mpz_inits(f1, f2, Z1, Z2, NULL);
    
    //gaussian init
    system_t s = init_gauss(v, pb_len+extra, pb_len);
    //printf("2^%d solutions to iterate\n", s->n2 - s->arb);
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


int main(int argc, char**argv){
    assert(argc == 4);
    srand(time(NULL));

    bool tests = atoi(argv[1]);
    int B = atoi(argv[3]);
    mpz_t N;
    mpz_init_set_str(N, argv[2], 10);

    clock_t t1 = clock();
    dixon(N, B, 1, tests);
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    if(!tests) printf("Total time: %fs\n", time_spent);

    mpz_clear(N);

    return 0;
}
