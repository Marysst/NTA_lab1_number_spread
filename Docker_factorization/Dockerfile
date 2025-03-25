# Використовуємо офіційний Python-образ
FROM python:3.11

# Встановлюємо робочу директорію в контейнері
WORKDIR /app

# Копіюємо всі файли в контейнер
COPY . .

# Встановлюємо необхідні бібліотеки
RUN pip install --no-cache-dir numpy sympy

# Вказуємо команду, яка запускатиме скрипт факторизації
CMD ["python", "factorization.py"]