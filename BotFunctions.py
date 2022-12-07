# -*- coding: utf-8 -*-
import numpy
import numpy as np
import math
from math import exp, log, sin, cos, tan, atan, sqrt, e, pi


#  Прямой ход метода Гаусса (у элементов матрицы тип обязательно float!)
def gauss(matrix):
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


#  Обратный ход метода Гаусса
def invgauss(matrix):
    for nrow in range(len(matrix)-1, 0, -1):
        row = matrix[nrow]
        for upper_row in matrix[:nrow]:
            factor = upper_row[nrow]
            upper_row -= factor*row
    return matrix


# sign
def sign(z):
    if z < 0:
        return -1
    if z > 0:
        return 1
    if z == 0:
        return 0


#  Мера обусловленности матрицы (p - размер матрицы)
def mcond(matrix, p):
    x = np.random.randint(1, 5, size=(p, 1))  # Произвольный вектор
    a = x.copy()
    s = 0
    for i in range(p):
        s = s + x[i, 0]**2
    NormaX = math.sqrt(s)
    epsilon = 1e-6  # Точность
    y = 1/NormaX * matrix.dot(x)  # Новый вектор (первый) из степенного метода
    z = x.transpose().dot(y)
    z = y - sign(z)*x
    s = 0
    for i in range(p):
        s = s + z[i, 0]**2
    NormaX = math.sqrt(s)  # Норма вектора, сравниваемого с epsilon
    if NormaX > epsilon:
        x = y.copy()
    while NormaX > epsilon:   # СТЕПЕННОЙ МЕТОД
        s = 0
        for i in range(p):
           s = s + x[i, 0]**2
        NormaX = math.sqrt(s)
        y = 1/NormaX * matrix.dot(x)
        z = x.transpose().dot(y)
        z = y - sign(z)*x
        s = 0
        for i in range(p):
            s = s + z[i, 0]**2
        NormaX = math.sqrt(s)  # Норма вектора, сравниваемого с epsilon
        if NormaX > epsilon:
            x = y.copy()
    s = 0
    for i in range(p):
        s = s + x[i, 0]**2
    nA = abs((y[0, 0]/x[0, 0]*math.sqrt(s)))  # Норма матрицы А
    x = a.copy()  # Скопированный вначале вектор
    B = np.linalg.inv(matrix)  # Обратная матрица
    s = 0
    for i in range(p):
        s = s + x[i, 0]**2
    NormaX = math.sqrt(s)
    y = 1/NormaX * B.dot(x)
    z = x.transpose().dot(y)
    z = y - sign(z)*x
    s = 0
    for i in range(p):
        s = s + z[i, 0]**2
    NormaX = math.sqrt(s)
    if NormaX > epsilon:
        x = y.copy()
    while NormaX > epsilon:   # СТЕПЕННОЙ МЕТОД
        s = 0
        for i in range(p):
            s = s + x[i, 0]**2
        NormaX = math.sqrt(s)
        y = 1/NormaX * B.dot(x)
        z = x.transpose().dot(y)
        z = y - sign(z)*x
        s = 0
        for i in range(p):
            s = s + z[i, 0]**2
        NormaX = math.sqrt(s)
        if NormaX > epsilon:
            x = y.copy()
    s = 0
    for i in range(p):
        s = s + x[i, 0]**2
    nAm = abs((y[0, 0]/x[0, 0])*math.sqrt(s))  # Норма обратной матрицы
    return nA*nAm


# Решение СЛАУ с трехдиагональной матрицей методом прогонки (А - матрица + столбец свободных эл; а, b - размер матрицы)
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


# Решение проблемы собственных значений и собственных векторов методом вращения Якоби (А - матрица СИММЕТРИЧНАЯ;
# b - ранг матрицы)
def JacobiRotation(A, b):
        T2 = np.zeros((b, b))  # Будущая матрица собственных векторов
        for i in range(b):
            T2[i, i] = 1.0

        norma = 0  # Норма
        eps = 1e-6  # Epsilon

        for z in range(b):     # Вычисление нормы Фробениуса
            for x in range(b):
                if z != x:
                    norma = norma + A[z, x]**2
        norma = norma**(1/2)
        #A1 = A.copy()
        while norma > eps:      # Метод вращения Якоби
            T1 = np.zeros((b, b))  # Ортогональная матрица U
            for i in range(b):
                T1[i, i] = 1.0
            k = 0  # Номер строки начального главного элемента
            q = 1  # Номер столбца начального главного элемента
            a = abs(A[0, 1])  # Начальный главный элемент
            for i in range(b):    # Поиск главного элемента и его положения
                for j in range(b):
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
            for z in range(b):
                for x in range(b):
                    if z != x:
                        norma = norma + A[z, x]**2
            norma = norma**(1/2)    # Норма Фробениуса

        # Деление каждого собственного вектора на его последний элемент
        for i in range(b):
            for j in range(b):
                T2[i, j] = (T2[i, j])/(T2[3, j])

        for i in range(0, 4):
            print('Собственное число', i + 1, ':', round(A[i, i], 7))
            print('Собственный вектор', i + 1, ':', round(T2[0, i], 7), round(T2[1, i], 7), round(T2[2, i], 7),
                  round(T2[3, i], 7))
            print('')


# Метод наименьших квадратов (кубическая функция) по 11 точкам и значениям
# в них - массив Tabulation (РАЗМЕР Tabulation - 2 строки и 11 столбцов!)
def MIN_SQUARE(Tabulation):
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
    #  Кубическая функция
    B2 = np.zeros((4, 5))
    min_square(4, 5, Tabulation, B2)
    gauss(B2)
    invgauss(B2)
    print(f'ф3(x)={B2[0, 4]: .6f} +{B2[1, 4]: .6f}x {B2[2, 4]: .6f}x^2 +{B2[3, 4]: .6f}x^3')


# Численное интегрирование - метод Симпсона (a, b - начало и конец отрезка интегрирования; formula - подыинтегральная
# функция - ВСЕ ВВОДИТСЯ ПОЛЬЗОВАТЕЛЕМ)
def Simpson(a, b, formula):
    if '^' in formula:
        formula = formula.replace('^', '**')
    def SimpsonInt(h, n, x):  # Метод Симпсона
        Integral = 0
        for i in range(n):
                Integral = Integral + h/3*eval(formula)
                x = x + h
                Integral = Integral + h/3*4*eval(formula)
                x = x + h
                Integral = Integral + h/3*eval(formula)
        return Integral
    eps = 1e-6  # Точность
    n = 8  # Начальное разбиение
    k = 4  # Для прямоугольников и трапеций k=2, для Симпсона k=4
    h = (b-a)/(2*n)  # Начальная длина шага
    Integral1 = SimpsonInt(h, n, a)
    n = n*2
    h = (b-a)/(2*n)
    Integral2 = SimpsonInt(h, n, a)
    while (abs(Integral2 - Integral1)/(2**k - 1) > eps):  # Сравнение точности по Рунге
        Integral1 = Integral2
        n = n*2
        h = (b-a)/(2*n)
        Integral2 = SimpsonInt(h, n, a)
    return Integral2 + abs(Integral2 - Integral1)/(2**k - 1)  # Уточнение по Ричардсону


#  Решение систем линейных алгебраических уравнений с помощью LU-разложения
def LU(matrix):
    # Ax = f при этом matrix = (Af)
    matrix_size = numpy.shape(matrix)
    A = np.zeros((matrix_size[0], matrix_size[0]))
    f = np.zeros((matrix_size[0], 1))
    for i in range(matrix_size[0]):
        f[i] = matrix[i, matrix_size[1] - 1]
        for j in range(matrix_size[1] - 1):
            A[i, j] = matrix[i, j]
    for i in range(matrix_size[0]):
        d = A[:i+1, : i+1]
        if np.linalg.det(d) == 0:
            print('Невозможно применить LU-разложение, так как есть нулевые главные миноры')
            return
    U = np.zeros(numpy.shape(A))
    L = np.zeros(numpy.shape(A))
    for i in range(matrix_size[0]):
        L[i, i] = 1
    U[0, 0] = A[0, 0]
    for i in range(1, matrix_size[0]):
        U[0, i] = A[0, i]
        L[i, 0] = A[i, 0]/U[0, 0]
    for i in range(1, matrix_size[0]):
        lu = 0
        for k in range(i):
            lu = lu + L[i, k]*U[k, i]
        U[i, i] = A[i, i] - lu
    for i in range(1, matrix_size[0]):
        for j in range(i+1, matrix_size[0]):
            lu = 0
            lu1 = 0
            for k in range(i):
                lu = lu + L[i, k]*U[k, j]
                lu1 = lu1 + L[j, k]*U[k, i]
            U[i, j] = A[i, j] - lu
            L[j, i] = 1/U[i, i]*(A[j, i] - lu1)
    y = np.linalg.inv(L).dot(f)
    x = np.linalg.inv(U).dot(y)
    return x


#  Определитель матрицы
def determinant(matrix):
    return np.linalg.det(matrix)


#  Решение СЛАУ методом Крамера
def Kramer(matrix):
    matrix_size = numpy.shape(matrix)
    A = np.zeros((matrix_size[0], matrix_size[0]))
    f = np.zeros((matrix_size[0], 1))
    for i in range(matrix_size[0]):
        f[i] = matrix[i, matrix_size[1] - 1]
        for j in range(matrix_size[1] - 1):
            A[i, j] = matrix[i, j]
    D = np.linalg.det(A)  # Главный определитель
    if D == 0:  # Невозможность решения
        return
    list_of_Dn = []
    for i in range(matrix_size[0]):
        Dn = A.copy()
        Dn[:, i] = f[:, 0]
        list_of_Dn.append(np.linalg.det(Dn))
    solution = []
    for i in range(matrix_size[0]):
        solution.append(list_of_Dn[i]/D)
    return solution


#  Решение СЛАУ с СИММЕТРИЧЕСКОЙ  матрицей методом квадратного корня (аналог LU для сим. матрицы)
def sqrt_method(matrix):
    matrix_size = numpy.shape(matrix)
    A = np.zeros((matrix_size[0], matrix_size[0]))
    f = np.zeros((matrix_size[0], 1))
    for i in range(matrix_size[0]):
        f[i] = matrix[i, matrix_size[1] - 1]
        for j in range(matrix_size[1] - 1):
            A[i, j] = matrix[i, j]
    for i in range(matrix_size[0]):
        d = A[:i+1, : i+1]
        if np.linalg.det(d) == 0:
            print('Невозможно применить, так как есть нулевые главные миноры')
            return
    D = np.zeros((matrix_size[0], matrix_size[0]))
    S = np.zeros((matrix_size[0], matrix_size[0]))
    for i in range(matrix_size[0]):
        sd = 0
        for k in range(i):
            sd = sd + abs(S[k, i])**2*D[k, k]
        D[i, i] = sign(A[i, i] - sd)
        S[i, i] = sqrt(abs(A[i, i] - sd))
        for j in range(i + 1, matrix_size[0]):
            ssd = 0
            for k in range(i):
                ssd = ssd + S[k, i]*S[k, j]*D[k, k]
            S[i, j] = (A[i, j] - ssd)/(S[i, i]*D[i, i])
    y = np.linalg.inv(np.transpose(S)).dot(f)
    x = np.linalg.inv(D.dot(S)).dot(y)
    return x


#  Количество решений СЛАУ с квадртаной матрицей
def value_of_solution(matrix):
    matrix_size = numpy.shape(matrix)
    A = np.zeros((matrix_size[0], matrix_size[0]))
    f = np.zeros((matrix_size[0], 1))
    for i in range(matrix_size[0]):
        f[i] = matrix[i, matrix_size[1] - 1]
        for j in range(matrix_size[1] - 1):
            A[i, j] = matrix[i, j]
    D = np.linalg.det(A)  # Главный определитель
    list_of_Dn = []
    for i in range(matrix_size[0]):
        Dn = A.copy()
        Dn[:, i] = f[:, 0]
        list_of_Dn.append(np.linalg.det(Dn))
    summa = 0
    for i in list_of_Dn:
        summa = summa + abs(i)
    if summa == 0 and D == 0:
        print('Система имеет бесконечно много решений')
    elif D == 0 and summa != 0:
        print('Система не имеет решений')
    else:
        print('Система имеет единственное решение')


#  Решение СЛАУ с неквадартаными матрицами, у которых могут быть: 0 и бесконечно много
#  решений (теорема Кронекера-Капелли)
def non_square_SLAU(matrix):
    size = np.shape(matrix)
    matrix_size = numpy.shape(matrix)
    A = np.zeros((matrix_size[0], matrix_size[1] - 1))
    for i in range(matrix_size[0]):
        for j in range(matrix_size[1] - 1):
            A[i, j] = matrix[i, j]
    rank1 = np.linalg.matrix_rank(A)
    rank2 = np.linalg.matrix_rank(matrix)
    # Теорема Кронекера-Капелли]
    if rank1 != rank2:
        print('Система несовместна')
        return
    for i in range(rank2):
        a1 = matrix[i, i]
        if a1 != 0:
            for j in range(size[1]):
                matrix[i, j] = matrix[i, j]/a1
            for i1 in range(i+1, size[0]):
                a = matrix[i1, i]
                for j in range(i, size[1]):
                    matrix[i1, j] = matrix[i1, j] - a*matrix[i, j]
        else:
            j = i + 1
            while matrix[i, j] == 0:
                j = j + 1
            a1 = matrix[i, j]
            for j in range(size[1]):
                matrix[i, j] = matrix[i, j]/a1
            for i1 in range(i+1, size[0]):
                a = matrix[i1, i]
                for j in range(i, size[1]):
                    matrix[i1, j] = matrix[i1, j] - a*matrix[i, j]

    for i in range(rank1):
        for i1 in range(i + 1, rank1):
            a = matrix[i, i1]
            for j in range(i1, size[1]):
                matrix[i, j] = matrix[i, j] - a*matrix[i1, j]
    solution = np.zeros((size[1] - 1, size[1] - rank1))
    for i in range(size[0]):
        for j in range(rank1, size[1]):
            solution[i, j - rank1] = matrix[i, j]
    for i in range(size[1] - 1):
        for j in range(size[1] - 1 - rank1):
            solution[i, j] = -1*solution[i, j]
    for i in range(size[1] - rank1 - 1):
        solution[i + rank1, i] = 1
    return solution


m = np.array([[1.0, -2.0, 5.0, -3.0, 0.0],
              [2.0, -3.0, 8.0, -5.0, 1.0],
              [0.0, 1.0, -2.0, 1.0, 1.0],
              [-1.0, 1.0, -3.0, 2.0, -1.0]])
m1 = np.array([[2.0, 1.0, 1.0, -2.0, 4.0, 1.0],
               [13.0, 8.0, 4.0, -3.0, 6.0, 9.0],
               [5.0, 4.0, 2.0, -3.0, 6.0, 3.0],
               [3.0, 2.0, 1.0, -1.0, 2.0, 2.0]])

m2 = np.array([[4.0, 2.0, 5.0],
              [1.0, 3.0, 1.0]])
print(non_square_SLAU(m))
formula = 'sin(x)*cos(x)*tan(x)'
print(f'{Simpson(0, 1, formula): .6f}')
