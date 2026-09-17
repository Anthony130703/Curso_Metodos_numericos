import matplotlib.pyplot as plt
from modelos import X

r = 1
n = 1000
x0 = [0.45, 0.5, 0.55]

Historial = {}

for i in x0:
    Historial[i] = X(n, i, r)
    
plt.figure(figsize=(12, 6))
for i, valores in Historial.items():
    eje_x = list(range(len(valores)))
    plt.scatter(eje_x, valores, s=2, label=f"x_0 = {i}")
    
plt.title("Comparación del Mapa Logístico para diferentes x0", fontsize=14)
plt.xlabel("Número de Iteración (n)", fontsize=12)
plt.ylabel("Valor de X", fontsize=12)
plt.yscale("log")
plt.grid(True, linestyle='--', alpha = 0.5)
plt.legend()
plt.show()