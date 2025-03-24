import random, math, sympy, itertools
import numpy as np

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def legendre_symbol(a, p):
    if p <= 0 or p % 2 == 0:
        raise ValueError("p має бути непарним простим числом")

    a = a % p
    if a == 0:
        return 0
    if a == 1:
        return 1
    if a % 2 == 0:
        return legendre_symbol(a // 2, p) * (-1 if (p % 8 in (3, 5)) else 1)
    return legendre_symbol(p % a, a) * (-1 if (a % 4 == 3 and p % 4 == 3) else 1)

def solovay_strassen(n, k=5):
    if n < 2:
        return False
    if n in (2, 3):
        return True
    if n % 2 == 0:
        return False

    for _ in range(k):
        a = random.randint(2, n - 1)
        g = gcd(a, n)
        if g > 1:
            return False  

        jacobian = legendre_symbol(a, n) % n
        mod_exp = pow(a, (n - 1) // 2, n)

        if jacobian != mod_exp:
            return False  

    return True  
    
def sieve(limit):
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    primes = []
    
    for num in range(2, limit + 1):
        if is_prime[num]:
            primes.append(num)
            for multiple in range(num * num, limit + 1, num):
                is_prime[multiple] = False
    return primes

def trial_division(n, limit=1000):
    if n < 2:
        return []
    
    primes = sieve(limit)
    factors = []

    for p in primes:
        if p * p > n:
            break
        while n % p == 0:
            factors.append(p)
            n //= p

    if n > 1:
        factors.append(n)  

    return factors

def f(x, n):
    return (x * x + 1) % n

def pollard_rho(n):
    if n % 2 == 0:
        return 2  

    for _ in range(5):  
        x = random.randint(2, n - 1)
        y = x
        d = 1

        while d == 1:
            x = f(x, n)
            y = f(f(y, n), n)
            d = gcd(abs(x - y), n)

        if 1 < d < n:
            return d  

    return None 

def full_pollard_factorization(n):
    factors = []

    def factorize(n):
        if n <= 1:
            return
        if solovay_strassen(n, 10):  
            factors.append(n)
            return
        
        divisor = None
        while not divisor:  
            divisor = pollard_rho(n)
            if divisor is None:
                print(f"Метод Полларда не зміг знайти дільник для {n}")
                return
        
        factors.append(divisor)
        factorize(n // divisor) 

    factorize(n)
    return sorted(factors)

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

def canonical_factorization(n):
    result = []
    
    def factorize(n):
        if solovay_strassen(n, 10):
            result.append(n)
            return True
        
        divisor = trial_division(n)
        if len(divisor) > 1:
            for i in range(len(divisor) - 1):
                result.append(divisor[i])
                n //= divisor[i]
            return factorize(n)
        
        divisor = pollard_rho(n)
        if divisor:
            result.append(divisor)
            n //= divisor
            if solovay_strassen(n, 10):
                result.append(n)
                return True
        
        divisor1, divisor2 = drillhart_morrison(n)
        if divisor1 and divisor2:
            if divisor1 == 1 or divisor2 == 1:
                result.extend([divisor1, divisor2])
                return False
            result.extend([divisor1, divisor2])
            return True
        
        print("Я не можу знайти канонічний розклад числа :(")
        return False
    
    if factorize(n):
        print(f"Канонічний розклад: {result}")
    else:
        print(f"Я знайшов такі дільники: {result}")
        