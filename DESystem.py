# -*- coding: utf-8 -*-
import numpy as np
from prettytable import PrettyTable


# f1(x1, x2, t)
def f1(u2, x):
    return -1*u2/x


# f2(x2, x2, t)
def f2(u1, x):
    return -1*u1/x


n = int(input('Введите количество разбиений '))  # Количесвто разбиений
K1 = np.zeros((2, 1))  # К1
K2 = np.zeros((2, 1))  # К2
K3 = np.zeros((2, 1))  # К3
K4 = np.zeros((2, 1))  # К4
a = 1  # Отрезок
b = 3
h = (b-a)/n  # Шаг
U = np.array([[2.0],  # Задача Коши
              [0.0]])
table = PrettyTable()  # Таблица
table.field_names = ['x', 'U1чис', 'U1ан', 'U2числ', 'U2ан']
#  Метод Рунге-Кутты 4го порядка
for i in range(n + 1):
    table.add_row([f'{a + i*h: .1f}', f'{U[0, 0]}', f'{a + i*h + 1/(a + i*h)}', f'{U[1, 0]}',
                   f'{-1*(a + i*h) + 1/(a + i*h)}'])
    K1[0, 0] = h * f1(U[1, 0], a + i*h)
    K1[1, 0] = h * f2(U[0, 0], a + i*h)
    K2[0, 0] = h * f1(U[1, 0] + K1[1, 0]/2, a + i*h + h/2)
    K2[1, 0] = h * f2(U[0, 0] + K1[0, 0]/2, a + i*h + h/2)
    K3[0, 0] = h * f1(U[1, 0] + K2[1, 0]/2, a + i*h + h/2)
    K3[1, 0] = h * f2(U[0, 0] + K2[0, 0]/2, a + i*h + h/2)
    K4[0, 0] = h * f1(U[1, 0] + K3[1, 0], a + i*h + h)
    K4[1, 0] = h * f2(U[0, 0] + K3[0, 0], a + i*h + h)
    U[0, 0] = U[0, 0] + 1/6 * (K1[0, 0] + 2*K2[0, 0] + 2*K3[0, 0] + K4[0, 0])
    U[1, 0] = U[1, 0] + 1/6 * (K1[1, 0] + 2*K2[1, 0] + 2*K3[1, 0] + K4[1, 0])

print(table)
