import numpy as np
import matplotlib.pyplot as plt
from modelos import biseccion, falsaPosicion, puntoFijo, newtonRaphson, secante

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
    
#Prueba de la secante
raiz_sec, errores_sec, iter_sec = secante(f, 4.0, 6.0, tolerancia)
print("\n--- Secante ---")
print(f"Raíz (x): {raiz_sec}")
print(f"Iteraciones: {iter_sec}")
if errores_sec:
    print(f"Error final: {errores_sec[-1]}")
    
# === GRÁFICA DE CONVERGENCIA ===
# Configuramos el tamaño de la ventana
plt.figure(figsize=(10, 6))

# Trazamos cada línea de error vs la iteración (el eje X es el rango de la longitud de la lista)
plt.plot(range(1, len(errores_bis) + 1), errores_bis, marker='o', label='Bisección')
plt.plot(range(1, len(errores_fp) + 1), errores_fp, marker='s', label='Falsa Posición')
plt.plot(range(1, len(errores_pf) + 1), errores_pf, marker='^', label='Punto Fijo')
plt.plot(range(1, len(errores_nr) + 1), errores_nr, marker='D', label='Newton-Raphson')
plt.plot(range(1, len(errores_sec) + 1), errores_sec, marker='x', label='Secante')

# Configuración visual
plt.yscale('log')  # Escala logarítmica para ver la caída del error claramente
plt.title('Velocidad de Convergencia de los Métodos Numéricos', fontsize=14, fontweight='bold')
plt.xlabel('Número de Iteraciones', fontsize=12)
plt.ylabel('Error Relativo Aproximado (Ea)', fontsize=12)
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.legend()

# Guardamos la imagen en alta calidad (DPI 300) ideal para reportes y documentos
#plt.savefig('grafica_convergencia.png', dpi=300, bbox_inches='tight')

plt.show()