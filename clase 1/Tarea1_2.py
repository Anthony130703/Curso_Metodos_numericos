from modelos import X
import matplotlib.pyplot as plt

x0 = 0.5
r = 1.0
n = 3000
Historial = {}
while r <= 4.0:
    r_limpio = round(r, 1)
    Historial[r_limpio] = X(n, x0, r_limpio)
    r += 0.1
    
plt.figure(figsize=(12, 6))

todos_los_x = []
todos_los_y = []

for i, valores in Historial.items():
    puntos_usados = valores[1000:]
    eje_x = [i] * len(puntos_usados)
    todos_los_x.extend(eje_x)
    todos_los_y.extend(puntos_usados)

plt.scatter(todos_los_x, todos_los_y, s=0.1, color = "black")
    
plt.title("Diagrama de Bifurcación del Mapa Logístico ($x_0 = 0.5$)", fontsize=14)
plt.xlabel("Parámetro de control (r)", fontsize=12)
plt.ylabel("Valores estables de X (Iteraciones 1000 a 3000)", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)

plt.show()