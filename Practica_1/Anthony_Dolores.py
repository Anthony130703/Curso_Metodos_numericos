import matplotlib.pyplot as plt
import numpy as np
from numpy import loadtxt

#Soluciones de la primera practica calificada

#Problema 1:
# a) primero obtengo un epsilon incial igual a 1 (ei = 1), tambien creo otras variables como epsilon nuevo (ef), test y 
#    una lista vacia.
#    mientras test es distinto de 1 entonces
#       ef = ei / 2.0
#       test = 1 + ef
#       iteracion += 1
#       añado a la lista el valor del ef
#    finalmente imprimo los valores de iteracion y el epsilon que dio antes de que test = 1.

# b) implementacion del pseudocodigo 
ei:float = 1.0
ef:float = 0.0
iteracion:int = 0
test:float = 0.0
valores = []

#Implementando el bucle
while test != 1.0:
    ef = ei / 2.0
    test = 1.0 + ef
    iteracion += 1
    valores.append(ef)
    ei = ef
    
print(f"El epsilon maquina obtenido es: {valores[-1]}")
print(f"Iteraciones realizadas: {iteracion}")

# c) Haciendo la grafica del R_k
epsilon:float = 0.0
Rk:float = 0.0
max_k:int = 60
valoresR_k = np.zeros(max_k +1) #Pre-asiganamos memoria del tamaño de 60 y lo llenamos de ceros
valoresK = np.zeros(max_k + 1)   #Pre-asiganamos memoria del tamaño de 60 y lo llenamos de ceros
for i in range(max_k + 1):
    #Calculando el valor del epsilon para ese valor de k
    epsilon = 2**(-i)
    
    #Obteniendo el valor de R_k
    Rk = ((1 + epsilon) - 1) / epsilon
    
    #Guardando los valores obtenidos en sus respectivas listas para luego usarlo en la grafica
    valoresR_k[i] = Rk
    valoresK[i] = i
    
plt.figure(figsize=(12, 6))
plt.scatter(valoresK, valoresR_k, s=2, color = "black")
plt.title("$R_k vs$ $k$", fontsize=14)
plt.xlabel("Parámetro de control (k)", fontsize=12)
plt.ylabel("Valores de $R_k$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()

# d)El valor de R_k deja de valor 1 porque este ya supero el valor de 52 bits que puede almacenar la matisa
#   y por ende solo almacene los valores de 0.000 y los demas numeros se pierden.
#   La relacion de la matisa con el numero de iteraciones, es que por cada iteracion se va consumiendo un bit
#   mas y asi sucesivamente, y por lo tanto solo puede almacenar hasta 52 bits, para la iteracion numero 53
#   estos valores ya no se almacenan.

#Problema 4:

#Creando la funcion para la posicion X de la pelota en funcion del tiempo y su angulo
def xPelota(angulo:float, t:float):
    xf = 37.0 * np.cos((angulo * np.pi) / 180.0) * t
    return xf
   
#Creando la funcion para la posicion Y de la pelota en funcion del tiempo y su angulo 
def yPelota(angulo:float, t:int):
    yf = 37.0 * np.sin((angulo * np.pi) / 180.0)*t  - ((0.5)*(9.81)* (t**2))
    return yf

#Creando la funcion que pondremos a iterar con el metodo de la biseccion
def f(tiempo:float):
    return (0.15)*((37.0 * np.cos((67.5 * np.pi) / 180.0) * tiempo)**2) * np.exp(-(0.04)*(37.0 * np.cos((67.5 * np.pi) / 180.0) * tiempo)) - (37.0 * np.sin((67.5 * np.pi) / 180.0)*tiempo  - ((0.5)*(9.81)* (tiempo**2)))

#Metodo de la biseccion
def biseccion(funcion, xl:float, xu:float, tolerancia, max_iter = 100):
    #Para verificar si los puntos escogidos estan bien
    if funcion(xl) * funcion(xu) > 0:
        print("Error: El intervalo inicial no garantiza una raíz.")
        return None, None, None

    #Declarando las varibales a ultizar
    iteracion:int = 0
    ea:float = 1.0
    xr_old:float = 0.0
    errores_historial = []

    while ea > tolerancia and iteracion < max_iter:
        #calculando el punto medio del intervalo
        xr = (xl + xu)/2.0
        
        #Obteniendo el error relativo
        if iteracion > 0:
            ea = abs((xr - xr_old)/xr) * 100
            errores_historial.append(ea)
        
        #determinando en que lado del intervalo se encuentra la raiz
        test = funcion(xl) * funcion(xr)
        
        if test < 0:
            xu = xr
        elif test > 0:
            xl = xr
        else: 
            ea = 0.0
            
        xr_old = xr
        iteracion += 1
    
    return xr, errores_historial, iteracion

t_lower = 5.0
t_upper = 6.0
tolerancia = 0.05

t_bis, errores_bis, iter_bis = biseccion(f, t_lower, t_upper, tolerancia)
raiz_bis = xPelota(67.5, t_bis)

# a) Encontrando el x_imp con un error menor del 5%
print(f"El tiempo es: {t_bis}")
print(f"Raíz (x): {raiz_bis}")
print(f"Iteraciones: {iter_bis}")
if errores_bis:
    print(f"Error final: {errores_bis[-1]}")

# b) Realizando la grafica para x_imp en funcion de theta
valores_x = np.zeros(20)
valores_theta = np.zeros(20)

for i in range (20):
    theta = i + 50.0
    valores_x[i] = xPelota(theta, t_bis)
    valores_theta[i] = theta
    
plt.figure(figsize=(10, 6))
plt.scatter(valores_theta, valores_x, s=2, color = "black")
plt.title('$X_{imp} vs  Angulo$', fontsize=14)
plt.xlabel("$Angulo$", fontsize=12)
plt.ylabel("Valores de $X_{imp}$", fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()
    
#Problema 2:
def f(x):
    return (0.1) * np.cos((2*np.pi) * x)

def derivada_exacta(x):
    return (-0.1) * np.sin((2*np.pi) * x)

x0 = 0.3
valor_exacto = derivada_exacta(x0)
valores_h = np.zeros(16)

for i in range(16):
    n = i + 1
    h = 10**(-n)
    valores_h[i] = h

errores_adelante = []
errores_centrada = []

for h in valores_h:
    # Método 1: Diferencia hacia adelante O(h)
    der_adelante = (f(x0 + h) - f(x0)) / h
    errores_adelante.append(abs(der_adelante - valor_exacto)/ valor_exacto)
    
    # Método 2: Diferencia centrada O(h^2)
    der_centrada = (f(x0 + h) - f(x0 - h)) / (2 * h)
    errores_centrada.append(abs(der_centrada - valor_exacto) / valor_exacto)
    
matriz = np.array([
    errores_adelante,
    errores_centrada
])

np.savetxt("Datos_item_a_problema_2", matriz, fmt ="%.4f")