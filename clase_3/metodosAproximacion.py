import numpy as np
from modelos import biseccion

def f(x):
    return 5 * np.exp(-x) + x - 5


x_lower = 1.6
x_upper = 6
tolerancia = 1e-6

raiz, errores, interaciones = biseccion(f, x_lower, x_upper, tolerancia)

print(f"Raíz encontrada (x): {raiz}")
if errores:
    print(f"Error final alcanzado: {errores[-1]}")
print(f"El numero de iteraciones: {interaciones}")