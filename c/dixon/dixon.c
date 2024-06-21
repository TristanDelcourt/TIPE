#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>

void print_list(int* l, int n){
    for(int i = 0; i<n; i++){
        printf("%d ", l[i]);
    }
    printf("\n");
}

void print_ll(int** ll, int n1, int n2){
    for(int i = 0; i<n1; i++){
        print_list(ll[i], n2);
    }
    printf("\n");
}

bool isPrime(int n) 
{ 
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
        if (isPrime(i)) 
            k++;; 
    } 
    return k;
} 

int* primes(int piB, int n){
    int* p = malloc(piB*sizeof(piB));
    int k = 0;
    for (int i = 2; i <= n; i++) { 
        if (isPrime(i)){
            p[k] = i;
            k++;
        }
    }
    return p;

}


/* DIXON */

bool factorise(mpz_t n, int* v, int piB, int* primes){
    for(int i = 0; i<piB; i++){
        v[i] = 0;
    }
    
    for(int i = 0; i<piB && mpz_cmp_ui(n, 1); i++){
        while (mpz_congruent_ui_p(n, 0, primes[i])){
            mpz_divexact_ui(n, n, primes[i]);
            v[i]++;
        }
    }

    if(!mpz_cmp_ui(n, 1))
        return true;
    return false;
      
    /*
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
    */
}

/*
int* get_vi(mpz_t zi, mpz_t N, int piB, mpz_t x){
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
*/


int** get_all_zi(mpz_t* z, mpz_t N, int piB, int* primes){
    mpz_t x;
    mpz_init(x);
    mpz_root(x, N, 2);
    mpz_t x_2_n;
    mpz_init(x_2_n);

    int** v = malloc((piB+1)*sizeof(int*));

    for(int i = 0; i < piB+1; i++){
        bool found = false;
        int* vi = malloc(piB*sizeof(int));

        while(!found){
            mpz_add_ui(x, x, 1);
            mpz_mul(x_2_n, x, x);
            mpz_sub(x_2_n, x_2_n, N);
            //gmp_printf("%Zd\n", x_2_n);
            found = factorise(x_2_n, vi, piB, primes);
        }

        v[i] = vi;
        mpz_set(z[i], x);
        //gmp_printf("%Zd\n", z[i]);
        //print_list(vi, piB);
        //printf("\n\n");
    }

    return v;
}

void swap_lines(int** v, int i, int j){
    int* temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}

int find_index(int** v, int from, int piB){
    for(int i = from; i < piB+1; i++){
        if(v[i][from]){
            return i;
        }
    }
    return -1;
}

unsigned modulo( int value, unsigned m) {
    int mod = value % (int)m;
    if (mod < 0) {
        mod += m;
    }
    return mod;
}

void mod_vect(int* v, int mod, int piB){
    for(int i = 0; i<piB; i++){
        v[i] = modulo(v[i], mod);
    }
}

void sub_vect(int** v, int i, int j, int piB){
    for(int k = 0; k<piB; k++){
        v[i][k] = v[i][k] - v[j][k];
    }
}

int* gaussian_solve(int** v, int piB){
    printf("Initial vectors\n");
    print_ll(v, piB+1, piB);
    
    for(int i = 0; i<piB+1; i++){
        mod_vect(v[i], 2, piB);
    }

    printf("Modded\n");
    print_ll(v, piB+1, piB);

    for(int i = 0; i<piB; i++){
        int k = find_index(v, i, piB);
        if(k != -1){
            swap_lines(v, i, k);

            for(int j = i + 1; j < piB+1; j++){
                if(v[j][i] == 1){
                    sub_vect(v, j, i, piB);
                    mod_vect(v[j], 2, piB);
                }
            }
        }
    }
    printf("Triangulate\n");
    print_ll(v, piB+1, piB);

    int* sol = malloc((piB+1)*sizeof(int));



}

void dixon(mpz_t N, int B){
    int piB = pi(B);
    int* p = primes(piB, B);
    //print_list(p, piB);

    mpz_t* z = malloc((piB+1)*sizeof(mpz_t));
    for(int i = 0; i < piB +1; i++){
        mpz_init(z[i]);
    }

    int** v = get_all_zi(z, N, piB, p);

    int* solution = gaussian_solve(v, piB);

}

void main(int argc, char**argv){
    assert(argc == 2);

    int B = 50;
    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    dixon(N, B);
}