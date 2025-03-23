import random
import math

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
