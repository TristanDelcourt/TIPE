import random
from math import prod, gcd, sqrt, ceil
from numpy import average
import sys
from copy import deepcopy
import timeit
from threading import Thread


B = 1000
MAX_TRIES = 20


# Renvoie la liste des nombres premiers inférieurs à B, ainsi que leur nombre
def prime_list():
    P = [2]
    i = 2
    while i <= B:
        if all(map(lambda x: i % x != 0, P)):
            P.append(i)
        i += 1
    return (P, len(P))


P, n = prime_list()


# Factorisation naïve
def naive(N):
    P = [2]
    d = {}
    i = 2
    while N > 1:
        if all(map(lambda x: i % x != 0, P)):
            P.append(i)
            while N % i == 0:
                if i not in d.keys():
                    d[i] = (True, 0)
                d[i] = (True, d[i][1] + 1)
                N //= i
        i += 1
    return d


# Renvoie la factorisation de n sous forme de vecteur de Z^n avec True ou une factorisation partielle avec False si n n'est pas B-friable
def basis_factorization(k: int) -> tuple[list[int], bool]:
    kk = k
    a = dict()
    for i in range(n):
        a[i] = 0
    i = 0
    while i < n and kk > 1:
        if kk % P[i] == 0:
            a[i] += 1
            kk //= P[i]
        else:
            i += 1
    # v = (list(a.values()), True) if kk == 1 else []
    # if kk == 1:
    #     assert prod(P[j] ** v[j] for j in range(n)) == k
    return list(a.values()), (True if kk == 1 else False)


# Renvoie S un ensemble de triplets (z_j, z_j^2 mod N, v_j) avec z_j un nombre tels que z_j^2 mod N est B-friable, v_j un vecteur de Z^n associé à la décomposition en facteurs premiers de z_j^2 mod N. (Recherche aléatoire)
def dixon_sieve(N):
    S = []
    z = random.randint(1, N)
    while len(S) < n + 1:
        if z not in S:
            zz = z**2 % N
            d = basis_factorization(zz)
            if d[1]:
                S.append((z, zz, d[0]))
        z = random.randint(1, N)
    return S


# Renvoie S un ensemble de triplets (z_j, z_j^2 mod N, v_j) avec z_j un nombre tels que z_j^2 mod N est B-friable, v_j un vecteur de Z^n associé à la décomposition en facteurs premiers de z_j^2 mod N. (Recherche sous la forme x^2 - n pour x >= floor(sqrt(n)))
def kraitchnik_sieve(N):
    S = []
    x = ceil(sqrt(N))
    while len(S) < n + 1:
        z = (x**2 - N) % N
        d = basis_factorization(z)
        if d[1]:
            S.append((x, z, d[0]))
        x += random.randint(1, 10)
    return S


def gauss_elim(M, Y):
    assert len(M[0]) > len(M)  # Pas de solution garantie sinon
    M = deepcopy(M)
    assert len(M) > 0
    assert len(Y) == len(M)

    # "Triangularisation" du système
    k, l = 0, 0  # Position du pivot
    while l < len(M[0]):
        # IdB: la matrice carrée issue de M avec les l premières colonnes est échelonnée en ligne
        # i0 est la première ligne >= k pour laquelle M[i_0][k] est non-nulle (si elle existe)
        i0 = k
        while i0 < len(M) and M[i0][l] == 0:
            i0 += 1
        if i0 < len(M):
            # M[k] <-> M[i0] et Y[k] <-> Y[i0]
            M[k], M[i0] = M[i0], M[k]
            Y[k], Y[i0] = Y[i0], Y[k]
            # M[i] <- M[i] - M[i][l] * M[k] et Y[i] <- Y[i] - M[i][l] * Y[k] pour i > k
            for i in range(k + 1, len(M)):
                a = M[i][l]
                Y[i] = abs(Y[i] - a * Y[k]) % 2
                for j in range(len(M[0])):
                    M[i][j] = abs(M[i][j] - a * M[k][j]) % 2
            k += 1
        l += 1

    # Solutions
    X = [-1 for _ in range(len(M[0]))]
    for i in range(len(M) - 1, -1, -1):  # On remonte les lignes depuis le bas
        j = 0
        while j < len(M[0]) and M[i][j] == 0:
            j += 1
        # j est la colonne de la variable principale
        if j < len(M[0]):
            X[j] = Y[i]
            # Variables libres
            for jj in range(len(M[0]) - 1, j, -1):
                if X[jj] == -1:
                    X[jj] = random.choice([0, 1])
                X[j] -= M[i][jj] * X[jj]
                X[j] = abs(X[j]) % 2
                # X[j] %= 2
    # Les colonnes entièrement nulles ont été ignorées
    for j in range(len(X)):
        if X[j] == -1:
            X[j] = random.choice([0, 1])

    # assert all(
    #     [sum([M[i][j] * X[j] for j in range(len(X))]) % 2 == 0 for i in range(n)]
    # )
    return X


# Renvoie une partie {(z_1, v_1),..., (z_r, v_r)} de S telle que v_1 + ... + v_r = 0
def nullify(S):
    for i in range(len(S)):
        if [S[i][2][j] % 2 for j in range(len(S[0][2]))] == [
            0 for j in range(len(S[0][2]))
        ]:
            return [i]
    M = [[S[i][2][j] % 2 for i in range(len(S))] for j in range(len(S[0][2]))]
    X = gauss_elim(M, [0 for _ in range(len(M))])
    A = [i for i in range(len(S)) if X[i] == 1]
    # assert all(
    #     [sum(S[i][2][j] for i in A for i in range(len(M))) % 2 == 0] for j in range(n)
    # )
    return A


# À partir d'une congruence de carrés, établit un facteur non-trivial de N. Si cela échoue, renvoie -1
def square_congruence(N, Z1, Z2):
    # assert Z1**2 % N == Z2**2 % N
    if Z1 % N == Z2 % N or Z1 % N == N - (Z2 % N):
        return -1
    return gcd(N, abs(Z1 - Z2))


# Applique l'algorithme de Dixon pour trouver un facteur non-trivial de N. Si cela échoue, renvoie -1
def dixon(N):
    S = dixon_sieve(N)
    A = nullify(S)
    Z1 = prod(S[i][0] for i in A) % N
    Z2 = prod(P[j] ** (sum(S[i][2][j] for i in A) // 2) for j in range(n)) % N
    return square_congruence(N, Z1, Z2)


# Applique l'algorithme de Dixon avec l'optimisation de Kraitchnik pour trouver un facteur non-trivial de N. Si cela échoue, renvoie -1
def kraitchnik(N):
    S = kraitchnik_sieve(N)
    A = nullify(S)
    Z1 = prod(S[i][0] for i in A) % N
    Z2 = prod(P[j] ** (sum(S[i][2][j] for i in A) // 2) for j in range(n)) % N
    # print(S, A, Z1, Z2)
    return square_congruence(N, Z1, Z2)


# Teste si k est probablement premier (True) ou composé (False)
def pseudo_prime_test(k):
    if k == 1:
        return False
    if k in P:
        return True
    a = random.randint(2, k - 1)
    return pow(a, k - 1, k) == 1


# Renvoie un dictionnaire contenant une factorisation quelconque de N dans un dictionnaire, de la forme f -> (b,p) où f est un facteur, b est si f n'est pas (a priori) composé, et p la multiplicité de f
def factor_once(N, factor_method):
    if N == 1:
        return {}
    if pseudo_prime_test(N):
        return {N: (True, 1)}
    d = {}
    F = basis_factorization(N)
    for i in range(n):
        if F[0][i] > 0:
            d[P[i]] = (True, F[0][i])
            N //= P[i] ** F[0][i]
    if d == {}:
        f = -1
        i = 0
        while f == -1 and i < MAX_TRIES:
            f = factor_method(N)
            i += 1
        if i >= MAX_TRIES:
            raise Exception(f"L'algorithme a échoué à factoriser {N}")
        # assert N % f == 0
        d[f] = (False, 1)
        d[N // f] = (False, 1)
    else:
        d[N] = (False, 1)
    return d


# Renvoie le nombre associé à la factorisation D
def number_of_dict(D):
    return prod([f**p for f, (b, p) in D.items()])


# Reprend chaque élement N du dictionnaire et le remplace par factor_once(N)
def factor_dict(D, factor_method):
    N1 = number_of_dict(D)
    I = [(f, (b, p)) for f, (b, p) in list(D.items()) if not b]
    for f, (b, p) in I:
        D1 = factor_once(f, factor_method)
        D.pop(f)
        for f1, (b1, p1) in D1.items():
            if f1 not in D.keys():
                D[f1] = (False, 0)
            D[f1] = (b1 or D[f1][0], D[f1][1] + p1)
        # assert number_of_dict(D) == N1


# Factorise N
def factor(N, factor_method, res1, res2, k):
    start = timeit.default_timer()

    D = {N: (False, 1)}
    while not all(map(lambda x: x[0], D.values())):
        # print(D)
        factor_dict(D, factor_method)
    
    end = timeit.default_timer()
    res1[k] = end-start
    res2[k] = D


def test(scale, size, factor_method):
    # TODO paralléliser
    L = [0 for _ in range(scale)]
    for i in range(scale):
        L[i] = random.randint(1, 9) * 10**size
        for j in range(size):
            L[i] += random.randint(0, 9) * 10**j
    
    
    threads = [None] * scale
    results_time = [None] * scale
    results_d = [None] * scale
    
    for i in range(scale):
        #print(L[i])
        threads[i] = Thread(target=factor, args=(L[i], factor_method, results_time, results_d, i))
        threads[i].start()
        #print(D)
    for i in range(len(threads)):
        threads[i].join()
    
    return results_time
    
    
    
    result[k] = round((end - start)/scale, 4)

if __name__ == "__main__":
    
    # scale = 20
    # threads = [None] * 10
    # result_time = [None] * 10
    # result_d = [None] * 10
    
    scale = 20
    random.seed()
    # assert len(sys.argv) > 0
    # N = int(sys.argv[1])
    for i in range(10, 20, 2):
        start = timeit.default_timer()
        
        print(
            f"Factorisation de nombres à {i} chiffres: {max(test(scale, i, kraitchnik))}s"
        )
        end = timeit.default_timer()
        print(end-start)
        # print(
        #    f"Factorisation de nombres à {i} chiffres: {round(test(scale, i, kraitchnik), 4)}s"
        #)