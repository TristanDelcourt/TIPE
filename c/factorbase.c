#include <stdbool.h>
#include <gmp.h>
#include <stdlib.h>

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
        if (is_prime(i)) k++;
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