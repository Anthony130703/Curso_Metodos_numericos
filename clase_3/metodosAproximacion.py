import numpy as np
from modelos import biseccion, falsaPosicion, puntoFijo, newtonRaphson

def f(x):
    return 5 * np.exp(-x) + x - 5

def g(x):
    return 5 - 5 * np.exp(-x)

def df(x):
    return -5 * np.exp(-x) + 1

x_lower = 1.6
x_upper = 6
tolerancia = 1e-6


print("=== COMPARACIÓN DE MÉTODOS ===")

#Prueba de Bisección
raiz_bis, errores_bis, iter_bis = biseccion(f, x_lower, x_upper, tolerancia)
print("\n--- Bisección ---")
print(f"Raíz (x): {raiz_bis}")
print(f"Iteraciones: {iter_bis}")
if errores_bis:
    print(f"Error final: {errores_bis[-1]}")

#Prueba de Falsa Posición
raiz_fp, errores_fp, iter_fp = falsaPosicion(f, x_lower, x_upper, tolerancia)
print("\n--- Falsa Posición ---")
print(f"Raíz (x): {raiz_fp}")
print(f"Iteraciones: {iter_fp}")
if errores_fp:
    print(f"Error final: {errores_fp[-1]}")

#Prueba de Punto Fijo 
raiz_pf, errores_pf, iter_pf = puntoFijo(g, x_lower, tolerancia)
print("\n--- Punto Fijo ---")
print(f"Raíz (x): {raiz_pf}")
print(f"Iteraciones: {iter_pf}")
if errores_pf:
    print(f"Error final: {errores_pf[-1]}")

#Prueba de Newton-Raphson
raiz_nr, errores_nr, iter_nr = newtonRaphson(f, df, x_upper, tolerancia)
print("\n--- Newton-Raphson ---")
print(f"Raíz (x): {raiz_nr}")
print(f"Iteraciones: {iter_nr}")
if errores_nr:
    print(f"Error final: {errores_nr[-1]}")