# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 20:11:33 2024

@author: pc
"""

import pandas as pd
import numpy as np
def cuadrados_minimos(x,y,n):
    A = np.vander(x,n+1) #tiene n+1 columnas porque tiene n+1 coeficientes (n es el grado)
    ATA = A.T@A
    ATb = A.T@y
    coeficientes = np.linalg.solve(ATA, ATb)
    return coeficientes
data = pd.read_csv("notas.csv")
A = data.to_numpy()
x = A[:,0]
y = A[:,1]
lineal = cuadrados_minimos(x,y,1)
cuadratico = cuadrados_minimos(x,y,2)
print("el ajuste lineal es",np.poly1d(lineal))
print("el ajuste cuadratico es", np.poly1d(cuadratico))
#Para encontrar el porcentaje necesario,tengo que resolver cuando y=8.
#Empiezo con el lineal, la ecuacion es 0.7253x+3.389 = 8
print("el porcentaje recomendable según el ajuste lineal es",(8-lineal[1])/lineal[0])
#Para el segundo tengo que resolver la ecuacion 0.0007126x^2 + 0.001366x+4.881 = 8 --> 0.0007126x^2 + 0.001366x -3.119 = 0. Busco las raices y devuelvo la que sea mayor a 0
raices = np.roots([cuadratico[0],cuadratico[1],cuadratico[2]-8])
for raiz in raices:
    if raiz > 0 :
        print("el porcentaje recomendable según el ajuste cuadratico es",raiz)
