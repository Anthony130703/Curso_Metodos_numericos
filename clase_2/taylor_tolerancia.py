import numpy as np
import matplotlib.pyplot as plt
from modelos import factorial

#Definiendo las variables a utilizar
x:float = 1.3
tolerancia = 1e-6
valor_real = np.sin(x)
aprox:float = 0.0
n:int = 0
error:float = 1.0

#Código para poder hacer una tabla con cada iteracion
print(f"Valor real de sin(1.3) = {valor_real}\n")
print(f"{'n':<5} | {'Aproximación':<15} | {'Error Real':<15}")
print("-" * 40)

#Desarrollando la sumatoria y verificando si la condicion se cumple para luego encontrar el valor de n 
while error >= tolerancia:
    #Realizando el calculo
    k = (2*n) + 1
    termino = (((-1)**n)*((x)**k)/ factorial(k))
    aprox += termino
    #Calculando el error
    error = np.abs(valor_real - aprox)
    #Mostrando los resultados de la iteracion en la tabla
    print(f"{n:<5} | {aprox:<15.8f} | {error:.10f}")
    if error >= tolerancia:
        n += 1
        
print("\n--- CONCLUSIÓN ---")
print(f"El bucle computacional se detuvo en n = {n}")