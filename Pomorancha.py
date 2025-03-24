import math, random, itertools
import numpy as np
from sympy import gcd, isprime, prime, Matrix, legendre_symbol, log

def gauss_elimination_F2(A):
    A = A.T
    m, n = A.shape
    cur_row = 0
    free_variables = []
    
    # 1. Приведення матриці до ступінчастого вигляду
    for col in range(n):
        # Знаходимо перший рядок, де є одиниця в поточному стовпці
        pivot_row = -1
        for row in range(cur_row, m):
            if A[row, col] == 1:
                pivot_row = row
                break
        
        if pivot_row == -1:
            # Якщо такого рядка немає, додаємо стовпець до списку вільних змінних
            free_variables.append(col)
            continue
        
        # Якщо опорний рядок знайдений, міняємо його місцями з cur_row
        if pivot_row != cur_row:
            A[[cur_row, pivot_row]] = A[[pivot_row, cur_row]]
        
        # Оновлюємо всі інші рядки, які мають 1 в цьому стовпці
        for row in range(m):
            if row != cur_row and A[row, col] == 1:
                A[row] ^= A[cur_row]  # Операція XOR (додавання по модулю 2)
        
        # Переходимо до наступного рядка
        cur_row += 1
    
    # 2. Побудова базису розв'язків
    solutions = []
    for free_var in free_variables:
        # Створюємо вектор рішення з нульовими значеннями
        solution = np.zeros(n, dtype=int)
        
        # Встановлюємо 1 для вільної змінної
        solution[free_var] = 1
        
        # Для кожного опорного рядка визначаємо значення залежних змінних
        for row in range(cur_row):
            pivot_col = np.argmax(A[row])  # Пошук стовпця з одиницею
            if pivot_col < n:
                solution[pivot_col] = A[row, free_var]
        
        # Додаємо отриманий вектор рішення до множини розв'язків
        solutions.append(solution)
    
    return solutions

def quadratic_sieve(n, base_size, sieve_range, threshold):
    """ Реалізація квадратичного сита з обробкою помилок """
    
    if isprime(n):
        print("Помилка: число є простим, факторизація неможлива.")
        return None

    factor_base = [-1]
    for i in itertools.count(start=1, step=1):
        p = prime(i)
        if (p == 2 and n % 2 == 1) or (p != 2 and legendre_symbol(n, p) == 1):
            factor_base.append(p)
        if len(factor_base) == base_size:
            break
    
    if len(factor_base) < 2:
        print("Помилка: недостатньо чисел у факторній базі.")
        return None

    m = math.isqrt(n)
    factor_base_roots = {}
    for p in factor_base:
        if p > 0:
            for x in range(-sieve_range, sieve_range + 1):
                if (x + m) ** 2 % p == n % p:
                    if x in factor_base_roots.keys():
                        factor_base_roots[x].append(p)
                    else:
                        factor_base_roots[x] = [p]

    relations = []
    for x in range(-sieve_range, sieve_range + 1):
        a = x + m
        b = a ** 2 - n
        lg = log(abs(b), 10)
        that_we_subtract = sum(log(i, 10) for i in factor_base_roots.get(x, []))
        difference = lg - that_we_subtract

        if difference <= threshold:
            factorization = []
            original_b = b
            for p in factor_base:
                count = 0
                if p == -1:
                    if b < 0:
                        b //= p
                        count += 1
                else:
                    while b % p == 0:
                        b //= p
                        count += 1
                factorization.append(count % 2)
            if b == 1:
                relations.append((a, original_b, factorization))

    if not relations:
        print("Помилка: не знайдено жодного В-числа.")
        return None

    A = np.array([r[2] for r in relations])
    print(A)
    zero_sums = gauss_elimination_F2(A)

    if not zero_sums:
        print("Помилка: не знайдено жодної лінійної залежності.")
        return None

    for solution_system in zero_sums:
        x = y = 1
        system_number = -1
        for j in solution_system:
            system_number += 1
            if j == 1:
                x *= relations[system_number][0]
                y *= relations[system_number][1]

        if (x + y) % n == 0 or (x - y) % n == 0:
            continue

        factor_1 = gcd(x + y, n)
        factor_2 = gcd(x - y, n)

        print("Успіх! Знайдені фактори:", factor_1, factor_2)
        return factor_1, factor_2

    print("Не вдалося знайти нетривіальну факторизацію.")
    return None

def estimate_parameters(n):
    """Оцінює оптимальні параметри для квадратичного сита."""
    if n <= 0:
        raise ValueError("Число повинно бути додатним.")
    
    log_n = math.log(n)
    log_log_n = math.log(log_n)

    base_size = max(20, round(math.exp(0.3 * math.sqrt(log_n * log_log_n))))
    sieve_range = max(50, base_size * 2)
    threshold = max(5, round(log_n / 3, 2))
    
    return base_size, sieve_range, threshold