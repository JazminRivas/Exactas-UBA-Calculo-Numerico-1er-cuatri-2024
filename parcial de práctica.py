# -*- coding: utf-8 -*-
"""
Created on Sat Jun 15 16:45:46 2024

@author: pc
"""
#PARCIALITO PRACTICA
#EJ 1
#A
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from typing import Optional, Tuple
import time


def metodo_sor(A: np.ndarray, b: np.ndarray, tolerancia,w, x0: Optional[np.ndarray] = None, max_iter = 1000) -> np.ndarray:
  D = np.diag(np.diag(A))
  L = np.tril(A,-1)
  U = np.triu(A,1)
  D_inv = np.linalg.inv(D+(w*L))
  Bsor = D_inv @ ((1-w)*D - w*U)
  Csor = D_inv @ (w*b)
  if x0 is None:
      x0 = np.zeros(len(b))
  x = x0.copy()
  for i in range(max_iter):
        x_nuevo = Bsor @ x + Csor
        if np.linalg.norm(x_nuevo - x, 2) < tolerancia:
            return x_nuevo, i+1   #como arranca desde el 0, le sumo 1 y es la iteracion que hizo que converja
        x = x_nuevo
    
  return x, max_iter

  
x0 = np.array([1,1,1,1])
b= np.array([1,0,3,4])
valores_alpha = [1,3]
valores_w = [1,3/2,3]
for alpha in valores_alpha:
    A = np.array([[2*alpha,0,0,0],[-1,8,0,1],[-1,2,-6,1],[1,0,0,alpha]])
    for w in valores_w:
        resultado, iteraciones = metodo_sor(A,b,1e-10,w,x0)
        print("Para w=",w,"y alpha=",alpha,"el resultado del metodo sor es", resultado,"y la cantidad de iteraciones son",iteraciones)
        print("Para w=",w,"y alpha",alpha,"el resultado es",np.linalg.solve(A,b))
#justamente para w=3 este método NO converge, tira overflow.

#Parte B, necesito calcular el radio espectral de la matriz de iteracion en cada paso 
def matriz_sor(A: np.ndarray, b: np.ndarray,w):
    D = np.diag(np.diag(A))
    L = np.tril(A,-1)
    U = np.triu(A,1)
    D_inv = np.linalg.inv(D+(w*L))
    Bsor = D_inv @ ((1-w)*D - w*U)
    return Bsor
for alpha in valores_alpha:
    for w in valores_w:
        A = np.array([[2*alpha,0,0,0],[-1,8,0,1],[-1,2,-6,1],[1,0,0,alpha]])
        matriz_iteracion = matriz_sor(A,b,w)
        autovalores = np.linalg.eigvals(matriz_iteracion)
        modulos = abs(autovalores)
        print("el radio espectral de la matriz de iteracion con w=",w,"y alpha=",alpha,"es",max(modulos))

#Ej 2
#parte A
def productoria(Xs,i,X_Eval):
    producto = 1
    for j in range (0,len(Xs)):
       if j != i:
           pol = (X_Eval - Xs[j])/(Xs[i] - Xs[j])
           producto = producto*pol
    return producto

def polinomio_lagrange(Xs,Ys,X_Eval):
    listasuma = []
    for i in range (0,len(Xs)):
       sumatoria =  productoria(Xs,i,X_Eval)*Ys[i]
       listasuma.append(sumatoria)
    return sum(listasuma) 
def f(x):
    return np.exp(-(x-5)**2)
#defino la grilla donde voy a graficar la funcion original
x_ej2 = np.linspace(0,10,100)
y_ej2 = f(x_ej2)
valores_n = [5,10,15]
for n in valores_n:
    grilla_puntos = np.linspace(0,10,n+1)
    y = f(grilla_puntos)
    #interpolo usando lagrange
    coefs = [polinomio_lagrange(grilla_puntos,y,x_eval) for x_eval in x_ej2] #la funcion devuelve al polinomio interpolador de lagrange evaluado en un punto, y lo evaluo en cada uno de los 100 puntos de la grilla original para comaparar
    plt.scatter(grilla_puntos,y, label = "Datos a interpolar con (n="+str(n)+")")
    plt.plot(x_ej2,coefs,label='Interpolador Lagrange con (n=' + str(n) + ')')
#grafico la funcion original:
plt.plot(x_ej2, y_ej2, 'k--', label='Función original')
plt.title('Interpolación de f(x), ej 2 parte A')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
plt.show()
#a medida que aumenta n, s aproxima mejor en el intervalo [4,6], aunque fuera de ese intervalo oscila y aproxima mal

#Parte B, toca interpolar usando Tchebychev
def zeros_chebyshev(a, b, n):
    return (a + b)/2 + (b - a)/2 * np.cos((2*np.arange(n) + 1) * np.pi / (2*(n))) #el arange crea un vector que va del 0 al n-1 [0,...,n-1] y evalua en todos

a=0
b=10
for n in valores_n:
    ceros_tchebychev = zeros_chebyshev(a,b,n+1)
    print(ceros_tchebychev)
    #ahora evaluo a la f en los ceros de este polinomio
    y_tche = f(ceros_tchebychev)
    interpolo_tcheby = np.poly1d(np.polyfit(ceros_tchebychev, y_tche, n+1))(x_ej2)
    plt.scatter(ceros_tchebychev,y_tche, label = "Ceros de tchebychev a interpolar con n="+str(n))
    plt.plot(x_ej2,interpolo_tcheby,label='Interpolador en los ceros de Tchebychev con n=' + str(n))
#grafico la funcion original:
plt.plot(x_ej2, y_ej2, 'k--', label='Función original')
plt.title('Interpolación de f(x), ej 2 parte B', fontsize=8)
plt.xlabel('x')
plt.ylabel('f(x)')
plt.legend()
plt.grid(True)
plt.show()