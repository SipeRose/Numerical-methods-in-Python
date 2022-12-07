# -*- coding: utf-8 -*-
import math
import numpy as np
from matplotlib import pyplot as plt


def gauss(matrix):
    #  Прямой ход метода Гаусса
    for nrow in range(len(matrix)):
        pivot = nrow + np.argmax(abs(matrix[nrow:, nrow]))
        if pivot != nrow:
            matrix[[nrow, pivot]] = matrix[[pivot, nrow]]
        row = matrix[nrow]
        divider = row[nrow]
        if abs(divider) < 1e-10:
            raise ValueError(f"Матрица несовместна. Максимальный элемент в столбце {nrow}: {divider:.3g}")
        row /= divider
        for lower_row in matrix[nrow+1:]:
            factor = lower_row[nrow]
            lower_row -= factor*row
    return matrix


def make_identity(matrix):
    #  Обратный ход метода Гаусса
    for nrow in range(len(matrix)-1, 0, -1):
        row = matrix[nrow]
        for upper_row in matrix[:nrow]:
            factor = upper_row[nrow]
            upper_row -= factor*row
    return matrix


def min_square(p, q, tabulation, A):
    #  Метод наименьших квадратов
    for i1 in range(0, p):
        for j in range(0, p):
            A[i1, j] = 0
            for k in range(0, 11):
                A[i1, j] = A[i1, j] + tabulation[0, k]**(i1 + j)
    for i1 in range(0, p):
        A[i1, q - 1] = 0
        for k in range(0, 11):
            A[i1, q - 1] = A[i1, q - 1] + tabulation[0, k]**(i1) * tabulation[1, k]
    return A


a, b, n = -1, 1, 10  # Отрезок и количество шагов
h = (b - a)/n
Tabulation = np.zeros((2, 11))
for i in range(0, 11):
    if i == 5:
        a = 0
    Tabulation[0, i] = a  #  Табуляция
    Tabulation[1, i] = math.sqrt(a + 2) + math.tan(a)  #  Данная функция
    a = a + h

#  Линейная функция
B = np.zeros((2, 3))
min_square(2, 3, Tabulation, B)
gauss(B)
make_identity(B)
print(f'ф1(x)={B[0, 2]: .6f} +{B[1, 2]: .6f}x')

#  Квадратичная функция
B1 = np.zeros((3, 4))
min_square(3, 4, Tabulation, B1)
gauss(B1)
make_identity(B1)
print(f'ф2(x)={B1[0, 3]: .6f} +{B1[1, 3]: .6f}x {B1[2, 3]: .6f}x^2')

#  Кубическая функция
B2 = np.zeros((4, 5))
min_square(4, 5, Tabulation, B2)
gauss(B2)
make_identity(B2)
print(f'ф3(x)={B2[0, 4]: .6f} +{B2[1, 4]: .6f}x {B2[2, 4]: .6f}x^2 +{B2[3, 4]: .6f}x^3')

#  Массивы значений функций и невязок
F1, F2, F3, Eps1, Eps2, Eps3 = np.zeros(11), np.zeros(11), np.zeros(11), np.zeros(11), np.zeros(11), np.zeros(11)
sum1, sum2, sum3 = 0, 0, 0  # Суммы квадратов невязок
for i in range(0, 11):
    #  Табуляция новых функций
    F1[i] = B[0, 2] + B[1, 2] * Tabulation[0, i]
    F2[i] = B1[0, 3] + B1[1, 3] * Tabulation[0, i] + B1[2, 3] * Tabulation[0, i]**(2)
    F3[i] = B2[0, 4] + B2[1, 4] * Tabulation[0, i] + B2[2, 4] * Tabulation[0, i]**(2) + B2[3, 4] * Tabulation[0, i]**(3)
    #  Невязки
    Eps1[i], Eps2[i], Eps3[i] = F1[i] - Tabulation[1, i], F2[i] - Tabulation[1, i], F3[i] - Tabulation[1, i]
    # Суммы квадратов невязок
    sum1, sum2, sum3 = sum1 + Eps1[i]**(2), sum2 + Eps2[i]**(2), sum3 + Eps3[i]**(2)

print(' ______ ____________ __________ __________ __________ __________ __________ __________')
print('|  Xi  | Yi = F(Xi) |  Ф1(Xi)  |  Ф2(Xi)  |  Ф3(Xi)  |   eps1   |   eps2   |   eps3   |')
print('|------|------------|----------|----------|----------|----------|----------|----------|')

for i in range(0, 11):
    print(f'| {Tabulation[0, i]: .1f} |  {Tabulation[1, i]: .5f}  | {F1[i]: .5f} | {F2[i]: .5f} | {F3[i]: .5f}'
          f' | {Eps1[i]: .5f} | {Eps2[i]: .5f} | {Eps3[i]: .5f} |')

print('|______|____________|__________|__________|__________|__________|__________|__________|')
print(f'|           Сумма квадратов невязок                  |{sum1: .5f}  | {sum2: .5f} | {sum3: .5f} |')
print(' ------------------------------------------------------------------------------------- ')

#  Графическое представление
h1, a = h/10, -1
Xlabel, Ylabel1, Ylabel2, Ylabel3, Ylabel4 = np.zeros(101), np.zeros(101), np.zeros(101), np.zeros(101), np.zeros(101)
for i in range(0, 101):
    if i == 50:
        a = 0
    Xlabel[i] = a
    Ylabel1[i] = math.sqrt(Xlabel[i] + 2) + math.tan(Xlabel[i])
    Ylabel2[i] = B[0, 2] + B[1, 2] * Xlabel[i]
    Ylabel3[i] = B1[0, 3] + B1[1, 3] * Xlabel[i] + B1[2, 3] * Xlabel[i]**(2)
    Ylabel4[i] = B2[0, 4] + B2[1, 4] * Xlabel[i] + B2[2, 4] * Xlabel[i]**(2) + B2[3, 4] * Xlabel[i]**(3)
    a = a + h1
plt.plot(Xlabel, Ylabel1)
plt.plot(Xlabel, Ylabel2)
plt.plot(Xlabel, Ylabel3)
plt.plot(Xlabel, Ylabel4)
plt.plot([-1, 1], [0, 0], 'black')
plt.plot([0, 0], [-1, 4], 'black')
for i in range(0, 11):
    plt.plot([Tabulation[0, i], Tabulation[0, i]], [Tabulation[1, i] - 0.3, Tabulation[1, i] + 0.3], 'black')
plt.xticks([-1, -0.5, 0, 0.5, 1])
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Метод наименьших квадратов')
plt.legend(['sqrt(x + 2) + tg(x)', '1.395426 + 1.726465x', '1.414691 + 1.726465x -0.048164x^2',
            '1.414691 + 1.270182x -0.048164x^2 + 0.640847x^3', 'Невязки'])
plt.show()
