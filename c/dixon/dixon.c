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

bool euler_criterion(mpz_t n, int p){
    int e = (p-1)/2;
    mpz_t r, p1;
    mpz_init(r);
    mpz_init_set_ui(p1, p);
    mpz_powm_ui(r, n, e, p1);
    return(mpz_cmp_ui(r, 1) == 0);
}

int* prime_base(mpz_t n, int* pb_len, int* primes, int piB){
    int* pb = malloc(piB*sizeof(int));
    pb[0] = 2;

    int j = 1;
    mpz_t r, p1;
    mpz_inits(r, p1, NULL);
    for(int i = 1; i<piB; i++){
        int e = (primes[i]-1)/2;
        mpz_set_ui(p1, primes[i]);
        mpz_powm_ui(r, n, e, p1);
        if(mpz_cmp_ui(r, 1) == 0){
            //printf("%d %d\n", primes[i], j);
            pb[j] = primes[i];
            j++;
        }
    }

    printf("base reduction %f%%\n", (float)j/piB*100);
    *pb_len = j;
    pb = realloc(pb, j*sizeof(int));
    
    mpz_clears(r, p1, NULL);
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
        if(s->sol[s->perm[i]]){
            add_vect(sum, v[s->perm[i]], s->n1);
        }
    }
}


bool factorise(mpz_t n, int* v, int pb_len, int* pb){
    for(int i = 0; i<pb_len; i++){
        v[i] = 0;
    }
    
    for(int i = 0; i<pb_len && mpz_cmp_ui(n, 1); i++){
        while (mpz_divisible_ui_p(n, pb[i])){
            v[i]++;
            mpz_divexact_ui(n, n, pb[i]);
        }
    }

    if(!mpz_cmp_ui(n, 1))
        return true;
    return false;
}

int** get_all_zi(mpz_t* z, mpz_t N, int pb_len, int* pb, int extra){
    //ceil(sqrt(n))
    mpz_t sqrt_n;
    mpz_init(sqrt_n);
    mpz_root(sqrt_n, N, 2);
    mpz_add_ui(sqrt_n, sqrt_n, 1);

    mpz_t x;
    mpz_init_set_ui(x, 1);

    mpz_t zi;
    mpz_t qx;
    mpz_inits(zi, qx, NULL);

    int** v = malloc((pb_len+extra)*sizeof(int*));

    for(int i = 0; i < pb_len+extra; i++){
        bool found = false;
        int* vi = malloc(pb_len*sizeof(int));

        while(!found){
            mpz_add(zi, sqrt_n, x);            

            //Q(x)
            mpz_mul(qx, zi, zi);
            mpz_sub(qx, qx, N);

            found = factorise(qx, vi, pb_len, pb);

            mpz_add_ui(x, x, 1);
        }
        printf("\r");
        printf("%f%%", (float)i/(pb_len+extra-1)*100);
        fflush(stdout);
        
        v[i] = vi;
        mpz_set(z[i], zi);
        //gmp_printf("%Zd\n", z[i]);
        //print_list(vi, pb_len);
        //printf("\n\n");
    }
    printf("\n");

    mpz_clears(sqrt_n, zi, x, qx, NULL);


    return v;
}

/*
bool gen_solutions(int* sol, int* indices, int n){
    int i = 0;
    while(i<n && (sol[indices[i]] == 1)){
        sol[indices[i]] = 0;
        i++;   
    }
    if(i >= n){
        return true;
    }
    sol[indices[i]] = 1;
    return false;
}
*/

void dixon(mpz_t N, int B, int extra){
    int piB = pi(B);
    printf("pi(B) = %d\n", piB);
    int* p = primes(piB, B);

    int pb_len;
    int* pb = prime_base(N, &pb_len, p, piB);
    free(p);

    mpz_t* z = malloc((pb_len+extra)*sizeof(mpz_t));
    for(int i = 0; i < pb_len+extra; i++){
        mpz_init(z[i]);
    }
    
    //Getting zis
    clock_t t1 = clock();
    int** v = get_all_zi(z, N, pb_len, pb, extra);
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("Time to get zi: %fs\n", time_spent);

    //int* arbitrary_indices;
    //int len;
    
    mpz_t f1, f2, Z1, Z2;
    mpz_inits(f1, f2, Z1, Z2, NULL);
    
    //gaussian init
    system_t s = init_gauss(v, pb_len+extra, pb_len);
    print_list(s->perm, s->n2);
    printf("%d %d\n", s->n1, s->n2);

    for(int i = 0; i<s->arb; i++){
        printf("  ");
    }
    printf("|\n");

    int* sum = malloc(pb_len*sizeof(int));


    bool done = false;
    while(!done){

        gaussian_step(s);
        print_list(s->sol, pb_len+extra);

        sum_lignes(sum, v, s);
        print_list(sum, pb_len);
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
            fprintf(stderr, "ERROR: no solution for this set of zi\n");
            exit(1);
        }
    }
    free(sum);
    free(pb);
    free_ll(v, pb_len+extra, pb_len);
    for(int i = 0; i < pb_len+extra; i++){
        mpz_clear(z[i]);
    }
    free(z);

    gmp_printf("%Zd = 0 [%Zd]\n%Zd = 0 [%Zd]\n", N, f1, N, f2);

    assert(mpz_divisible_p(N, f1));
    assert(mpz_divisible_p(N, f2));

    mpz_clears(f1, f2, Z1, Z2, NULL);
}


int main(int argc, char**argv){
    assert(argc == 2);
    srand(time(NULL));

    int B = 20;
    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    clock_t t1 = clock();
    dixon(N, B, 5);
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("Total time: %fs\n", time_spent);

    mpz_clear(N);

    return 0;
}
