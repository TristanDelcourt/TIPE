from gmpy2 import mpz, isqrt_rem

def fermat_naive(n):
    n = mpz(n)
        
    SQRT, r = isqrt_rem(n)
    r *= -1
    
    u = 2 * SQRT + 1
    v = 1
    
    while r != 0:
        
        if r > 0:
            while r > 0:
                r -= v
                v += 2
        
        if r < 0:
            r += u
            u += 2
    
    a = (u + v -2)/2
    b = (u - v)/2
    
    print(f"n = {a} * {b}")