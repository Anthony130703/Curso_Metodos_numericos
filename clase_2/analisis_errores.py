import numpy as np
from modelos import factorial
import matplotlib.pyplot as plt

#Definiendo la suma truncada para el seno(x)
def serieTaylor(n:int, x:float) -> float:
    sumaTruncada = 0.0
    for i in range(0, n + 1):
        termino = (((-1)**i) * ((x)**(2*i + 1)))/factorial(2*i + 1)
        sumaTruncada += termino
    return sumaTruncada

#Calculando el error real ∣ aprox n − sin x ∣ para n = 1,2, . . . , 20
x:float = float(input('Ingrese el valor de x: '))
#Usamos int para asegurarnos que sea un numero entero
n_max = int(input("¿Hasta qué valor de n deseas evaluar?: "))
valorReal = np.sin(x)

#Pre-asginamos memoria de una tamaño y lo llenamos de ceros
n_valores = np.zeros(n_max)
errores = np.zeros(n_max)

for i in range (n_max):
    n = i + 1 #Con esto n ira desde 1 hasta n_max
    #errorCalculado = np.abs(valorReal - serieTaylor(n, x)) + 1e-16
    #errorCalculado = np.maximum(np.abs(valorReal - serieTaylor(n, x)), 1e-16)
    errorCalculado = np.abs(valorReal - serieTaylor(n, x))

    if errorCalculado < 1e-16:
        errorCalculado = 1e-16

    #Guardando los valores
    n_valores[i] = n
    errores[i] = errorCalculado

# Configuración de la gráfica semilogarítmica
plt.figure(figsize=(10, 6))
plt.semilogy(n_valores, errores, marker='o', color='red')

# Etiquetas y diseño visual
plt.title('Error de Truncamiento vs Redondeo (Serie de Taylor)')
plt.xlabel('Número de términos (n)')
plt.ylabel('Error Real |aprox - sin(x)| (Escala Logarítmica)')

# Activamos la cuadrícula tanto para las líneas mayores como menores
plt.grid(True, which="both", linestyle="--", alpha=0.7) 

# Lanzamos la ventana emergente con el gráfico
plt.show()