import numpy as np
import Spline
mat=np.array([[1,-9,0,0,63], #Исходная матрица
              [5,4,1,0,24],
              [0,-6,3,1,40],
              [0,0,-2,1,-11]])
alpha1=(-1)*mat[0, 1] #альфа(1)
beta1=mat[0,4]       #бета(1)
AlphaBeta = np.array([[alpha1, 0.0, 0.0], #Массив альфа(i) и бета(i)
                    [beta1, 0.0, 0.0]])
#Прямой ход метода прогонки
for i in range(1,3):
     AlphaBeta[0,i]=(-1)*(mat[i,i+1]/(mat[i,i-1]*AlphaBeta[0,i-1]+mat[i,i])) #Вычисление альфа(i)
     AlphaBeta[1,i]=(mat[i,4]-mat[i,i-1]*AlphaBeta[1,i-1])/(mat[i,i-1]*AlphaBeta[0,i-1]+mat[i,i]) #Вычисление бета(i)
#Обратный ход метода прогонки
y=[0.0,0.0,0.0,0.0]  #Вектор решений
k2=-mat[3,2] #Капа2
mu2=mat[3,4] #Мю2
y[3]=(k2*AlphaBeta[1,2]+mu2)/(1-k2*AlphaBeta[0,2]) #y(3)
for i in range(2, -1, -1):        #Вычисление y(i)
     y[i]=AlphaBeta[0,i]*y[i+1]+AlphaBeta[1,i]
print(y) #Вывод решения

U = Spline.progonka(mat, 4, 5)
print(U)