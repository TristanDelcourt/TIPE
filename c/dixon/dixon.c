#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <assert.h>
#include <time.h>

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

int* primes(int piB, int B){
    int* p = malloc(piB*sizeof(piB));
    int k = 0;
    for (int i = 2; i <= B; i++) { 
        if (isPrime(i)){
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
    for(int i = 1; i<piB; i++){
        if(euler_criterion(n, primes[i])){
            pb[j] = primes[i];
            j++;
        }
    }

    //printf("base reduction %f percent\n", (float)j/piB*100);
    *pb_len = j;
    pb = realloc(pb, j);
    return pb;
}

void swap_lines(int** v, int i, int j){
    int* temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}

int find_index(int** v, int from, int n1){
    for(int i = from; i < n1; i++){
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

void mod_vect(int* v, int mod, int n1){
    for(int i = 0; i<n1; i++){
        v[i] = modulo(v[i], mod);
    }
}

void div_vect(int* v, int d, int n1){
    for(int i = 0; i<n1; i++){
        v[i] /= d;
    }
}

void sub_vect(int** v, int i, int j, int n1){
    for(int k = 0; k<n1; k++){
        v[i][k] = v[i][k] - v[j][k];
    }
}

void prod_vect(mpz_t prod, mpz_t* v, int n1){
    mpz_set_ui(prod, 1);
    for(int i = 0; i<n1; i++){
        mpz_mul(prod, prod, v[i]);
    }
}

void rebuild(mpz_t prod, int* v, int* primes, int n1){
    mpz_set_ui(prod, 1);
    mpz_t temp;
    mpz_init(temp);
    for(int i = 0; i<n1; i++){
        mpz_ui_pow_ui(temp, primes[i], v[i]);
        mpz_mul(prod, prod, temp);
    }
}

int* sum_lignes(int** v, int n1, int n2, int* sol){
    int* sum = malloc(n1*sizeof(int));

    for(int i = 0; i<n2; i++){
        sum[i] = 0;
        for(int j = 0; j<n1; j++){
            if(sol[j]){ 
                sum[i] += v[j][i];
            }
        }
    }

    return sum;
}

int** transpose(int** v, int n1, int n2){
    int** m = malloc(n2*sizeof(int*));

    for(int i = 0; i<n2; i++){
        m[i] = malloc(n1*sizeof(int));
        for(int j = 0; j<n1; j++){
            m[i][j] = v[j][i];
        }
    }

    return m;
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

int** get_all_zi(mpz_t* z, mpz_t N, int pb_len, int* pb){
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

    int** v = malloc((pb_len+1)*sizeof(int*));

    for(int i = 0; i < pb_len+1; i++){
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

        v[i] = vi;
        mpz_set(z[i], zi);
        //gmp_printf("%Zd\n", z[i]);
        //print_list(vi, pb_len);
        //printf("\n\n");
    }

    return v;
}


int* gaussian_solve(int** v, int n1, int n2){
    //printf("Initial vectors\n");
    //print_ll(v, n1, n2);
    
    int** m = transpose(v, n1, n2);

    //printf("Transposed\n");
    //print_ll(m, n2, n1);
    
    for(int i = 0; i<n2; i++){
        mod_vect(m[i], 2, n1);
    }

    //printf("Modded\n");
    //print_ll(m, n2, n1);

    for(int i = 0; i<n2; i++){
        int k = find_index(m, i, n2);
        if(k != -1){
            swap_lines(m, i, k);

            for(int j = i + 1; j < n2; j++){
                if(m[j][i] == 1){
                    sub_vect(m, j, i, n1);
                    mod_vect(m[j], 2, n1);
                }
            }
        }
    }
    
    //printf("Triangulate\n");
    //print_ll(m, n2, n1);

    int* sol = malloc(n1*sizeof(int));
    for(int i = 0; i<n1; i++){
        sol[i] =-1;
    }

    for(int i = n2-1; i>-1; i--){
        int j = 0;
        while(j < n1 && !m[i][j]){
            j++;
        }

        if(j<n1){
            sol[j] = 0;

            for(int k = n1-1; k>j; k--){
                if(sol[k] == -1){
                    sol[k] = rand() % 2;
                }
                sol[j] -= m[i][k] * sol[k];
                sol[j] = abs(sol[j]) % 2;
            }
        }
    }

    for(int i = 0; i<n1; i++){
        if(sol[i] == -1){
            sol[i] = rand() % 2;
        }
    }
    return sol;
}

void dixon(mpz_t N, int B){
    int piB = pi(B);
    int* p = primes(piB, B);

    int pb_len;
    int* pb = prime_base(N, &pb_len, p, piB);
    //print_list(pb, pb_len);

    mpz_t* z = malloc((pb_len+1)*sizeof(mpz_t));
    for(int i = 0; i < pb_len +1; i++){
        mpz_init(z[i]);
    }

    //Getting zis
    clock_t t1 = clock();
    int** v = get_all_zi(z, N, pb_len, pb);
    clock_t t2 = clock();
    double time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("Time to get zi: %fs\n", time_spent);

    //Solving matrix
    t1 = clock();
    int* solution = gaussian_solve(v, pb_len+1, pb_len);
    t2 = clock();
    time_spent = (double)(t2 - t1) / CLOCKS_PER_SEC;
    printf("Time to get solve matrix: %fs\n", time_spent);


    int* sum = sum_lignes(v, pb_len+1, pb_len, solution);
    div_vect(sum, 2, pb_len);

    mpz_t Z1, Z2;
    mpz_inits(Z1, Z2, NULL);
    prod_vect(Z1, z, pb_len+1);
    rebuild(Z2, sum, p, pb_len);

    mpz_t f1, f2;
    mpz_inits(f1, f2, NULL);
    mpz_sub(f1, Z1, Z2);
    mpz_add(f2, Z1, Z2);

    mpz_gcd(f1, f1, N);
    mpz_gcd(f2, f2, N);

    gmp_printf("%Zd = 0 [%Zd]\n%Zd = 0 [%Zd]\n", N, f1, N, f2);

    assert(mpz_divisible_p(N, f1));
    assert(mpz_divisible_p(N, f2));
}


void main(int argc, char**argv){
    assert(argc == 2);
    srand(time(NULL));

    int B = 1000;
    mpz_t N;
    mpz_init_set_str(N, argv[1], 10);

    dixon(N, B);
}
