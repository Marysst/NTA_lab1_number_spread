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

def find_smooth_relations(n, factor_base):
    terms = []
    relations = []
    needed_relations = len(factor_base) + 1  
    b_values = set()
    b_2, b_1 = 1, 0
    v, u = 1, int(math.sqrt(n))
    a = u

    b_1, b_2 = b_2, (a * b_2 + b_1) % n
    print(b_1, b_2)
    process_b_value(b_2, n, factor_base, b_values, relations)
    
    while len(relations) < needed_relations:
        terms.append(a)
        v = (n - u ** 2) // v
        a = (int(math.sqrt(n)) + u) // v
        u = a * v - u
        b_1, b_2 = b_2, (a * b_2 + b_1) % n
        process_b_value(b_2, n, factor_base, b_values, relations)
    return relations