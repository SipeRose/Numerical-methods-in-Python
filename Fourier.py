# -*- coding: utf-8 -*-
import numpy as np
from math import sin, pi, cos
from scipy.fft import fft, ifft    # Встроенная функция быстрого преобразования Фурье
from matplotlib import pyplot as plt


def X(t, T, RND):  # Данный сигнал
    #x = sin(2 * pi * 50 * T[t]) + sin(2 * pi * 120 * T[t]) + 2*RND[t]
    x = sin(2*pi*50*(4 % 3 + 1)*T[t]) + sin(2*pi*120*(4 % 2 + 1)*T[t]) + sin(2*pi*30*(4 % 5 + 1)*T[t]) \
        + (4 % 3 + 1)*RND[t]
    return x


def fFouriert(k, Mas3):  # Дискретное быстрое преобразование Фурье
    N = 2**(k+1)
    for i in range(0, 512 - 2**(k+1) + 1, 2**(k+1)):
        z = 0
        while z < 2**k:
            Mas3[i+z, k] = Mas3[i+z, k-1] + (cos(-2*pi*z/N) + 1j*sin(-2*pi*z/N))*Mas3[i+z+2**k, k-1]
            Mas3[i+2**k+z, k] = Mas3[i+z, k-1] - (cos(-2*pi*z/N) + 1j*sin(-2*pi*z/N))*Mas3[i+z+2**k, k-1]
            z = z + 1


RND = np.random.randn(512)  # Случайная величина
T = np.zeros(512)  # Массив аргументов
for i in range(0, 512):
    T[i] = 0.001*i
Mas1 = np.zeros(512)  # Массив значений исходного сигнала
Mas2 = np.zeros((512, 2))  # Первый столбец - номер точки; второй - двоичный код номера
str1 = ''  # Заготовки под двоичный код
str2 = ''
for i in range(0, 512):
    Mas1[i] = X(i, T, RND)  # Заполнение массива значений исходного сигнала
    Mas2[i, 0] = i  # Номер точки
    str1 = bin(i)  # Двоичный код номера
    for k in range(-1, -1*(len(str1)-1), -1):  # Инверсия двоичного кода
        str2 = str2 + str1[k]
    while len(str2) < 9:  # Двоичный код содержит 9 (для 512 точек) цифр из {0, 1}
        str2 = str2 + '0'
    Mas2[i, 1] = int(str2, 2)  # Перевод двоичного кода в номер точки
    str2 = ''

# Дискретное преобразование Фурье
print('')
print('Дискретное преобразование Фурье')
for k in range(0, 512):
    x = 0
    for n in range(0, 512):
        x = x + Mas1[n] * (cos(-2*pi*k*n/512) + 1j*sin(-2*pi*k*n/512))
    if k in range(5) or k in range(507, 512):  # Вывод превых и последних 5 значений преобразоания
        print(f'X_{k} =', x)

# Быстрое преобразование Фурье
Mas3 = np.zeros((512, 9), dtype='complex_')  # Таблица для 9 шагов БПФ
for i in range(1, 512, 2):  # 1 шаг: вычисление X(0) и X(1)
    Mas3[i-1, 0] = Mas1[int(Mas2[i-1, 1])] + Mas1[int(Mas2[i, 1])]
    Mas3[i, 0] = Mas1[int(Mas2[i-1, 1])] - Mas1[int(Mas2[i, 1])]

for i in range(1, 9): # БПФ
    fFouriert(i, Mas3)
print('')
print('Быстрое преобразование Фурье')
for i in range(5):  # Вывод первых и последних пяти значений БПФ
    print(f'X_{i} =', Mas3[i, 8])
for i in range(507, 512):
    print(f'X_{i} = ', Mas3[i, 8])

# Быстрое преобразование Фурье через встроенную функцию
print('')
print('Быстрое преобразование Фурье через встроенную функцию')
Mas4 = fft(Mas1)
for i in range(5):
    print(f'X_{i} = ', Mas4[i])
for i in range(507, 512):
    print(f'X_{i} = ', Mas4[i])

# Исходный сигнал
print('')
print('Исходный сигнал')
for i in range(5):
    print(f'{i} {Mas1[i]}')
for i in range(507, 512):
    print(f'{i} {Mas1[i]}')

# Обратное дискретное преобразование Фурье
print('')
print('Исхдный сигнал через обратное дискретное преобразование Фурье')
for k in range(0, 512):
    x = 0
    for n in range(0, 512):
        x = x + (Mas3[n, 8] * (cos(2*pi*k*n/512) + 1j*sin(2*pi*k*n/512)))/512
    if k in range(5) or k in range(507, 512):
        print(f'x_{k} =', x)

# Обратное быстрое преобразование Фурье
Mas5 = ifft(Mas3[:, 8])
print('')
print('Исходный сигнал через обратное преобразование Фурье')
for i in range(5):
    print(f'x_{i} = ', Mas5[i])
for i in range(507, 512):
    print(f'x_{i} = ', Mas5[i])

# Обратное быстрое преобразование Фурье преобразования, полученного через встроенную функцию
Mas6 = ifft(Mas4)
print('')
print('Исходный сигнал через обратное преобразование Фурье с помощью встроенной функции')
for i in range(5):
    print(f'x_{i} = ', Mas6[i])
for i in range(507, 512):
    print(f'x_{i} = ', Mas6[i])

# Графическое представление
Xlabel = np.zeros(512)
for i in range(512):
    Xlabel[i] = i
plt.title('Преобразование Фурье')
plt.plot(Xlabel, Mas4)
plt.show()
for i in range(512):
    Mas4[i] = abs(Mas4[i])
plt.title('Спектральная плотность')
plt.plot(Xlabel, Mas4)
plt.show()
plt.title('Исходный сигнал')
plt.plot(Xlabel, Mas1)
plt.show()
