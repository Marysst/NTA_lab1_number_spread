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
    process_b_value(b_2, n, factor_base, b_values, relations)
    
    while len(relations) < needed_relations:
        terms.append(a)
        v = (n - u ** 2) // v
        a = (int(math.sqrt(n)) + u) // v
        u = a * v - u
        b_1, b_2 = b_2, (a * b_2 + b_1) % n
        process_b_value(b_2, n, factor_base, b_values, relations)
    return relations

def gauss_elimination_F2(A):
    A = A.T
    m, n = A.shape
    cur_row = 0
    free_variables = []
    
    for col in range(n):
        pivot_row = -1
        for row in range(cur_row, m):
            if A[row, col] == 1:
                pivot_row = row
                break
        
        if pivot_row == -1:
            free_variables.append(col)
            continue
        
        if pivot_row != cur_row:
            A[[cur_row, pivot_row]] = A[[pivot_row, cur_row]]
        
        for row in range(m):
            if row != cur_row and A[row, col] == 1:
                A[row] ^= A[cur_row]
        
        cur_row += 1
    
    solutions = []
    for free_var in free_variables:
        solution = np.zeros(n, dtype=int)
        solution[free_var] = 1
        for row in range(cur_row):
            pivot_col = np.argmax(A[row])
            if pivot_col < n:
                solution[pivot_col] = A[row, free_var]
        solutions.append(solution)
    return solutions

def compute_xy(zero_sums, relations, n):
    for i in zero_sums:
        x = y = 1
        count = 0
        for j in i:
            if j == 1:
                x *= relations[count][0]
                y *= relations[count][1]
            count += 1
        x = x % n
        root, exact = sympy.integer_nthroot(y, 2)
        if exact:
            return x, root
    print("No valid x, y found")
    return None, None

def drillhart_morrison(n):
    factor_base = build_factor_base(n)
    relations = find_smooth_relations(n, factor_base)
    A = np.array([r[2] for r in relations])
    zero_sums = gauss_elimination_F2(A)
    x, y = compute_xy(zero_sums, relations, n)
    
    if x is None or y is None:
        print("Factorization failed")
        return None, None
    
    r1 = sympy.gcd(int(x + y), n)
    r2 = sympy.gcd(int(x - y), n)
    return r1, r2

if __name__ == "__main__":
    n = 9073
    res = drillhart_morrison(n)
    print(n, "=", res[0], "*", res[1])