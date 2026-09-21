import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return np.sin(x)

def derivada_exacta(x):
    return np.cos(x)

x0 = 1.0
valor_exacto = derivada_exacta(x0)
valores_h = np.logspace(-1, -17, 100)

errores_adelante = []
errores_centrada = []

for h in valores_h:
    # Método 1: Diferencia hacia adelante O(h)
    der_adelante = (f(x0 + h) - f(x0)) / h
    errores_adelante.append(abs(der_adelante - valor_exacto))
    
    # Método 2: Diferencia centrada O(h^2)
    der_centrada = (f(x0 + h) - f(x0 - h)) / (2 * h)
    errores_centrada.append(abs(der_centrada - valor_exacto))

# Gráfica comparativa
plt.figure(figsize=(10, 6))

plt.plot(valores_h, errores_adelante, color='red', marker='.', linestyle='', label="Hacia adelante O(h)")
plt.plot(valores_h, errores_centrada, color='blue', marker='.', linestyle='', label="Centrada O(h^2)")

plt.xscale('log')
plt.yscale('log')
plt.gca().invert_xaxis()

plt.axvline(x=1e-8, color='red', linestyle='--', alpha=0.5, label="h óptimo Adelante (~10^-8)")
plt.axvline(x=1e-5, color='blue', linestyle='--', alpha=0.5, label="h óptimo Centrada (~10^-5)")

plt.title("Diferencia Hacia Adelante vs Centrada", fontsize=14)
plt.xlabel("Tamaño de h (Se hace más pequeño hacia la derecha →)", fontsize=12)
plt.ylabel("Error Absoluto", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend()
plt.show()