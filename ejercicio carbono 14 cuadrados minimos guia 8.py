# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 20:27:57 2024

@author: pc
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
data = pd.read_csv("carbono14.csv")
A = data.to_numpy()
print(A)
def cuadrados_minimos(x,y,n):
    A = np.vander(x,n+1) #tiene n+1 columnas porque tiene n+1 coeficientes (n es el grado)
    ATA = A.T@A
    ATb = A.T@y
    coeficientes = np.linalg.solve(ATA, ATb)
    return coeficientes
#la primera columna es el tiempo desde la muerte y la segunda es como decae el carbono 14
x = A[:,0]
y = A[:,1]
pol = cuadrados_minimos(x,y,1)

#defino la grilla para graficar el polinomio
x_grafico = np.linspace(0,14020,10000)
plt.plot(x_grafico,np.poly1d(pol)(x_grafico),color = "red", label ="polinomio lineal aproximador " )
plt.scatter(x,y, label = "datos")
plt.title('Carbono 14')
plt.xlabel('tiempo')
plt.ylabel("decadencia de C14")
plt.legend()
plt.grid(True)
plt.show() #es relativamente una buena aproximacion


#parte B, tengo que ajustar los datos para f(x) = ae^bx. si f(x) = ae^bx --> Ln(f(x))=Ln(a) + bx que es una recta
#le tomo logaritmo a todo
y_log = np.log(y)
print(y_log)
pol_parteB = cuadrados_minimos(x,y_log,1) #es una recta que devuelve Ln(a) + bx -- [b,ln(a)] la tupla
a = np.exp(pol_parteB[1]) #a Ln(a) lo vuelvo una exponencial, e^ln(a) = a y con eso recupero a
y_interpolado = a*np.exp(pol_parteB[0]*x_grafico)
plt.scatter(x, y, label='Datos originales')
plt.plot(x_grafico, y_interpolado, label="Aproximacion exponencial por cuadrados minimos", color='red')

# Configuración de la gráfica
plt.title('Ajuste exponencial usando mínimos cuadrados,parte B')
plt.xlabel('tiempo')
plt.ylabel('decadencia de C14')
plt.legend()
plt.grid(True)
plt.show()

def error(x,y,f): #la f sería tu polinomio aproximador
    error_total = 0
    for i in range(len(y)):
        error = (y[i] - f(x[i]))**2
        error_total += error
    return error_total
def f(x): #es mi polinomio interpolador
    return a*np.exp(pol_parteB[0]*x)
print("el error es", error(x,y,f))

#parte C, tengo que graficarlos al mismo tiempo
plt.scatter(x, y, label='Datos originales')
plt.plot(x_grafico,np.poly1d(pol)(x_grafico),color = "red", label ="polinomio lineal aproximador " )
plt.plot(x_grafico, y_interpolado, color = "green", label="Aproximacion exponencial por cuadrados minimos")
plt.title('Comparación de métodos,parte C')
plt.xlabel('tiempo')
plt.ylabel('decadencia de C14')
plt.legend()
plt.grid(True)
plt.show()