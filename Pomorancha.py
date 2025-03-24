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