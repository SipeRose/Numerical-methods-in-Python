# -*- coding: utf-8 -*-
import math
import numpy as np
A = np.array([[9.0, -6.0, 3.0, -5.0],   # Исходная матрица
              [-6.0, -9.0, 5.0, -1.0],
              [3.0, 5.0, 7.0, 0.0],
              [-5.0, -1.0, 0.0, 1.0]])
T2 = np.array([[1.0, 0.0, 0.0, 0.0],   # Будущая матрица собственных векторов
               [0.0, 1.0, 0.0, 0.0],
               [0.0, 0.0, 1.0, 0.0],
               [0.0, 0.0, 0.0, 1.0]])
norma = 0  # Норма
eps = 1e-6  # Epsilon

for z in range(0, 4):     # Вычисление нормы Фробениуса
    for x in range(0, 4):
        if z != x:
            norma = norma + A[z, x]**2
norma = norma**(1/2)

c1 = 0  # Счетчик итераций
while norma > eps:      # Метод вращения Якоби
    c1 = c1 + 1
    T1 = np.array([[1.0, 0.0, 0.0, 0.0],    # Ортогональная матрица U
                   [0.0, 1.0, 0.0, 0.0],
                   [0.0, 0.0, 1.0, 0.0],
                   [0.0, 0.0, 0.0, 1.0]])
    k = 0  # Номер строки начального главного элемента
    q = 1  # Номер столбца начального главного элемента
    a = abs(A[0, 1])  # Начальный главный элемент
    for i in range(0, 4):    # Поиск главного элемента и его положения
        for j in range(0, 4):
            if abs(A[i, j]) > a and i != j:
                a = abs(A[i, j])
                k = i
                q = j

    if k < q:
        i = k    # Номер строки итогового главного элемента
        j = q    # Номер столбца итогового главного элемента
    else:
        i = q
        j = k

    # Поиск "угла поворота"
    q = A[i, i] - A[j, j]
    if q == 0:
        phi = (2**(1/2))/2
    else:
        if i > j:
            phi = 1/2 * math.atan((2 * A[i, j])/(A[j, j] - A[i, i]))
        else:
            phi = 1/2 * math.atan((2 * A[i, j])/(A[i, i] - A[j, j]))

    c = math.cos(phi)   # cos(phi)
    s = math.sin(phi)   # sin(phi)

    T1[i, i] = c
    T1[i, j] = -1*s
    T1[j, i] = s
    T1[j, j] = c

    B = np.linalg.inv(T1)   # Обратная к U
    C = B.dot(A)    # ОбратнаяU*A
    A = C.dot(T1)   # ОбратнаяU*A*U
    T2 = T2.dot(T1)     # U = U1*U2*...*Un

    # Поиск нормы новой матрицы А
    norma = 0
    for z in range(0, 4):
        for x in range(0, 4):
            if z != x:
                norma = norma + A[z, x]**2
    norma = norma**(1/2)    # Норма Фробениуса

# Деление каждого собственного вектора на его последний элемент
for i in range(0, 4):
    for j in range(0, 4):
        T2[i, j] = (T2[i, j])/(T2[3, j])

# Вывод результатов
print('')
print('Количество итераций', c1)
print('''
        Собственные значения и соответствующие им собственные вектора: ''')
for i in range(0, 4):
    print('Собственное число', i + 1, ':', round(A[i, i], 7))
    print('Собственный вектор', i + 1, ':', round(T2[0, i], 7), round(T2[1, i], 7), round(T2[2, i], 7),
          round(T2[3, i], 7))
    print('')

# Через встроенную функцию
print('Через встроенную функцию')
print(np.linalg.eig(A)[0])