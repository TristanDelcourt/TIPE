from numpy import array, prod

def divide(f, d, i):
    i+=1
    p_i = d
    e_i = 1
    f = f//d
    while f%d == 0:
        e_i += 1
        f //= d
    return i, f, p_i, e_i

def print_factors(n, l):
    print("n = ", end='')
    for p, e in l:
        if e==1:
            print(f"{p} *", end=' ')
        else:
            print(f"{p}**{e} *", end=' ')
    

def trial_division_factors(n: int, *, max = 1_000_000, check = False):
    i = 0 # number of distinct prime factors
    f = n # still infactored portion
    
    factors = []
    
    for d in (2, 3):
        if f%d == 0:
            i, f, p_i, e_i = divide(f, d, i)
            factors.append((p_i, e_i))
            
    d = 5
    add = 2
    
    while d<=max and d*d <= f:
        if f%d == 0:
            i, f, p_i, e_i = divide(f, d, i)
            factors.append((p_i, e_i))
        d += add
        add = 6 - add
    
    if d*d>f :
        i += 1
        p_i = f
        e_i = 1
        factors.append((p_i, e_i))
    
    print_factors(n, factors)
    print(f)
    
    if check:
        if n != prod(array([pow(p, i) for (p, i) in factors])):
            print("Error in trial_division")


def trial_division_primetest(n: int, *, max = 1_000_000):
    i = 0 # number of distinct prime factors
    f = n # still infactored portion
        
    for d in (2, 3):
        if f%d == 0:
            return False

    d = 5
    add = 2
    
    while d<=max and d*d <= f:
        if f%d == 0:
            return False
        d += add
        add = 6 - add
    
    return True

def powmod(a, b, m):
    n = 1
    
    while b != 0:
        if b%2 == 1:
            n *= a%m
        b = b // 2
        a *= a%m
        
    return n

if __name__ == "__main__":
    trial_division_factors(360)