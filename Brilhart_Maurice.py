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