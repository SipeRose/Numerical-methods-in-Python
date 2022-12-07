# -*- coding: utf-8 -*-
import numpy as np
import math
from prettytable import PrettyTable
from matplotlib import pyplot as plt


def progonka(A, a, b):   # Метод прогонки
    alpha1 = (-1) * A[0, 1] #альфа(1)
    beta1 = A[0, b - 1]       #бета(1)
    AlphaBeta = np.zeros((2, a - 1))
    AlphaBeta[0, 0] = alpha1
    AlphaBeta[1, 0] = beta1

    # Прямой ход метода прогонки
    for i in range(1, a - 1):
        AlphaBeta[0, i] = (-1)*(A[i, i+1]/(A[i, i-1]*AlphaBeta[0, i-1]+A[i, i]))  # Вычисление альфа(i)
        AlphaBeta[1, i] = (A[i, b - 1]-A[i, i-1]*AlphaBeta[1, i-1])/(A[i, i-1]*AlphaBeta[0, i-1]+A[i, i])  # Вычисление бета(i)
    # Обратный ход метода прогонки
    y = np.zeros(a)  # Вектор решений
    k2 = -A[a - 1, b - 3]  # Капа2
    mu2 = A[a - 1, b - 1]  # Мю2
    y[a - 1] = (k2*AlphaBeta[1, a - 2]+mu2)/(1-k2*AlphaBeta[0, a - 2])
    for i in range(a - 2, -1, -1):        # Вычисление y(i)
        y[i] = AlphaBeta[0,i]*y[i+1] + AlphaBeta[1, i]
    return y  # Вывод решения


def Graph(X, Y1, Y2, h1, k):  # Графическое представление
    global a
    for i in range(0, 11):
        X[i] = a
        Y1[i] = math.sqrt(X[i] + 2) + math.tan(X[i])
        Y2[i] = Sp[0, k + 1] + Sp[1, k + 1]*(a - Tabulation[0, k]) + Sp[2, k + 1]*(a - Tabulation[0, k])**2 + Sp[3, k + 1]*(a - Tabulation[0, k])**3
        if i != 10:
            a = a + h1
    plt.plot(Xlabel, Ylabel1)
    plt.plot(Xlabel, Ylabel2)


table = PrettyTable()  # Таблица
table.field_names = ['X', 'f(x)', 'Sp(x)', '|f(x) - Sp(x)|']

A = np.zeros((11, 12))   # Загатовка под трехдиагональную матрицу

a, b, n = -1, 1, 10  # Отрезок и количество шагов
h = (b - a)/n  # Шаг

Tabulation = np.zeros((2, 11))  # Табулируемая функция (точки и значения)
for i in range(0, 11):  # Табуляция функции
    if i == 5:
        a = 0
    Tabulation[0, i] = a
    Tabulation[1, i] = math.sqrt(a + 2) + math.tan(a)  # Данная функция
    a = a + h  # Следующая точка

A[0, 0] = 1  # Элементы трехдиагоналной матрицы2
A[10, 10] = 1

for i in range(1, 10):  # Элементы трехдиагоналной матрицы
    A[i, i-1] = h/3
    A[i, i] = 4/3 * h
    A[i, i + 1] = h/3
    A[i, 11] = (Tabulation[1, i + 1] - 2 * Tabulation[1, i] + Tabulation[1, i - 1])/h  # "Свободные члены"
U = progonka(A, 11, 12)  # Решение методом прогонки
Sp = np.zeros((4, 11))  # Массив коэффициаентов (a, b, c, d)
for i in range(1, 11):  # Вычисление коэффициентов Ai и Ci
    Sp[2, i] = U[i-1]
    Sp[0, i] = Tabulation[1, i-1]

for i in range(1, 11):  # Вычисление коэффициентов Bi и Di
    if i == 10:
        Sp[1, i] = (Tabulation[1, i] - Tabulation[1, i - 1])/h - h/3*(2*Sp[2, i])
        Sp[3, i] = (-1 * Sp[2, i])/(3*h)
    else:
        Sp[1, i] = (Tabulation[1, i] - Tabulation[1, i - 1])/h - h/3*(2*Sp[2, i] + Sp[2, i+1])
        Sp[3, i] = (Sp[2, i+1] - Sp[2, i])/(3*h)

a1 = -1 + h/2  # середина интервала
y1 = 0
y2 = 0
for i in range(1, 11):  # Вычисление значений функций и сплайна в серединах отрезков
    y1 = math.sqrt(a1 + 2) + math.tan(a1)
    y2 = Sp[0, i] + Sp[1, i]*(a1 - Tabulation[0, i - 1]) + Sp[2, i]*(a1 - Tabulation[0, i - 1])**2 + Sp[3, i]*(a1 - Tabulation[0, i - 1])**3
    table.add_row([f'{a1: .1f}', f'{y1: .5f}', f'{y2: .5f}', f'{abs(y1-y2): .4f}'])
    a1 = a1 + h
print(table)

# Графическое представление
h1, a = h/10, -1
Xlabel, Ylabel1, Ylabel2 = np.zeros(11), np.zeros(11), np.zeros(11)
for k in range(0, 10):
    Graph(Xlabel, Ylabel1, Ylabel2, h1, k)
plt.show()
