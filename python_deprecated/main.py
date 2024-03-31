from prime_test import *
from fermat import *

def definitly_not_prime(n):
    trial_division_factors(n)

def main():
    n = int(input())
    print("--")
    
    # division test for primes up to a maximum
    if not(trial_division_primetest(n)):
        definitly_not_prime(n)
        return
    print("Passed trial division test")
    
    if not(powmod(2, n-1, n) != 1 and powmod(3, n-1, n) != 1
           and powmod(5, n-1, n) != 1 and powmod(7, n-1, n) != 1):
        
        definitly_not_prime(n)
        return 
    print("Passed probable prime test")
    
    fermat_naive(n)
    
    print("--")
    

    

if __name__ == "__main__":
    
    main()