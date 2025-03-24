import math, sympy, itertools
import numpy as np

def build_factor_base(n, a=1 / 2 ** 0.5):
    limit = math.exp((math.log(n) * math.log(math.log(n))) ** 0.5)
    factor_base = [-1]
    for i in itertools.count(1, 1):
        p = sympy.prime(i)
        if p >= limit ** a:
            break
        if p == 2 and n % 2 == 1 or sympy.legendre_symbol(n, p) == 1:
            factor_base.append(p)
    return factor_base

def process_b_value(b_2, n, factor_base, b_values, relations):
    b2 = (b_2 ** 2) % n
    if b2 > n // 2:
        b2 -= n
    original_b2 = b2
    factorization = []
    
    for p in factor_base:
        count = 0
        if p == -1:
            if b2 < 0:
                b2 //= p
                count += 1
        else:
            while b2 % p == 0:
                b2 //= p
                count += 1
        factorization.append(count % 2)
    
    if b2 == 1 and b_2 not in b_values:
        relations.append((b_2, original_b2, factorization))
        b_values.add(b_2)