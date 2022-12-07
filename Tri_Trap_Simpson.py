# -*- coding: utf-8 -*-
from math import exp, log


def TriInt(h, n, ai):  # Метод средних прямоугольников
    Integral = 0
    for i in range(n):
        Integral = Integral + h*log(10*ai)
        ai = ai + h
    return Integral


def TrapInt(h, n, ai):  # Метод трапеций
    Integral = 0
    for i in range(n):
        Integral = Integral + h/2*(log(10*ai) + log(10*(ai+h)))
        ai = ai + h
    return Integral


def SimpsonInt(h, n, ai):  # Метод Симпсона
    Integral = 0
    for i in range(n):
        Integral = Integral + h/3*(log(10*ai) + 4*log(10*(ai + h)) + log(10*(ai + 2*h)))
        ai = ai + 2*h
    return Integral


a = 0.1  # Начало отрезка
b = 0.1*exp(1)  # Конец отрезка
n = 8  # Начальное число разбиений
eps = 1e-6  # Точность
k = 2  # Для прямоугольников и трапеций k=2, для Симпсона k=4
h = (b-a)/n  # Шаг

# Метод прямоугольников
print('МЕТОД СРЕДНИХ ПРЯМОУГОЛЬНИКОВ')
ai = a + h/2  # Середина отрезка
Integral1 = TriInt(h, n, ai)
print(f'h = {h: .10f}')
print(f'I_h = {Integral1: .8f}')
n = n*2  # Увеличение числа разбиений
h = h/2  # Уменьшение длины шага
ai = a + h/2
Integral2 = TriInt(h, n, ai)
print(f'h/2 = {h: .10f}')
print(f'I_h/2 = {Integral2: .8f}')
print(f'Уточнение по Ричардсону: {Integral2 + abs(Integral2 - Integral1)/(2**k - 1): .8f}')
while abs(Integral2 - Integral1)/(2**k - 1) > eps:  # Сравнение точности по Рунге
    Integral1 = Integral2
    print(f'h = {h: .10f}')
    print('I_h = ', f'{Integral1: .8f}')
    n = n*2
    h = h/2
    ai = a + h/2
    Integral2 = TriInt(h, n, ai)
    print(f'h/2 = {h: .10f}')
    print(f'I_h/2 = {Integral2: .8f}')
    print(f'Уточнение по Ричардсону: {Integral2 + abs(Integral2 - Integral1)/(2**k - 1): .8f}')

# Метод трапеций
print('\nМЕТОД ТРАПЕЦИЙ')
n = 8
h = (b-a)/n
Integral1 = TrapInt(h, n, a)
print(f'h = {h: .10f}')
print(f'I_h = {Integral1: .8f}')
n = n*2
h = h/2
Integral2 = TrapInt(h, n, a)
print(f'h/2 = {h: .10f}')
print(f'I_h/2 = {Integral2: .8f}')
print(f'Уточнение по Ричардсону: {Integral2 + abs(Integral2 - Integral1)/(2**k - 1): .8f}')
while abs(Integral2 - Integral1)/(2**k - 1) > eps:  # Сравнение точности по Рунге
    Integral1 = Integral2
    print(f'h = {h: .10f}')
    print('I_h = ', f'{Integral1: .8f}')
    n = n*2
    h = h/2
    Integral2 = TriInt(h, n, a)
    print(f'h/2 = {h: .10f}')
    print(f'I_h/2 = {Integral2: .8f}')
    print(f'Уточнение по Ричардсону: {Integral2 + abs(Integral2 - Integral1)/(2**k - 1): .8f}')

# Метод Симпсона
print('\nМЕТОД СИМПСОНА')
n = 8  # Начальное разбиение
k = 4
h = (b-a)/(2*n)  # Начальная длина шага
Integral1 = SimpsonInt(h, n, a)
print(f'h = {h: .10f}')
print(f'I_h = {Integral1: .8f}')
n = n*2
h = (b-a)/(2*n)
Integral2 = SimpsonInt(h, n, a)
print(f'h/2 = {h: .10f}')
print(f'I_h/2 = {Integral2: .8f}')
print(f'Уточнение по Ричардсону: {Integral2 + abs(Integral2 - Integral1)/(2**k - 1): .8f}')
while (abs(Integral2 - Integral1)/(2**k - 1) > eps):  # Сравнение точности по Рунге
    Integral1 = Integral2
    print(f'h = {h: .10f}')
    print('I_h = ', f'{Integral1: .8f}')
    n = n*2
    h = (b-a)/(2*n)
    Integral2 = SimpsonInt(h, n, a)
    print(f'h/2 = {h: .10f}')
    print(f'I_h/2 = {Integral2: .8f}')
    print(f'Уточнение по Ричардсону: {Integral2 + abs(Integral2 - Integral1)/(2**k - 1): .8f}')