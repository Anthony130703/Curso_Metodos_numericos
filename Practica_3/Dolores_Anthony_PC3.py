"""
Universidad Nacional de Ingeniería
Facultad de Ciencias - Escuela Profesional de Física
CF060 Métodos Numéricos - Ciclo 2026-II
Práctica Calificada N.° 3 - Miércoles 18.11.2026, 12:00-16:00 h
Profesor: PhD. Alejandro Paredes

Alumno: Dolores Narciso Anthony Jesus
Código: 20242212D
"""
#!/usr/bin/python3

from numpy import *

import matplotlib.pyplot as plt

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import matplotlib.pyplot as plt
import math

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
            ea = abs((xr - xr_old)/xr)
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

def sustitucion_forward(L, b):
    #Resolviendo Lx = b para L triangular inferior
    n = len(b)
    x = np.zeros(n, dtype=float)
    for k in range(n):
        if abs(L[k][k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = 0
        for i in range(k):
            suma = suma + L[k][i]*x[i]
        x[k] = (b[k] - suma) / L[k][k]
    return x

def cholesky(A):
    n = len(A)
    #Creando una matriz G llena de ceros
    G = [[0.0] * n for _ in range(n)]
    #Recorremos fila por fila
    for i in range(n):
        #PASO 1: Calculando el elemento de la diagonal G[i][i]
        suma_cuadrados = 0.0
        #haciendo que k recorra las filas estrictamente por encima de i
        for k in range(i):
            suma_cuadrados += G[k][i] * G[k][i]
        valor_raiz = A[i][i] - suma_cuadrados
        #Validando: si es negativo o cero, no es Positiva Definida
        if valor_raiz <= 0:
            raise ValueError("La matriz no es Simétrica Positiva Definida.")
        G[i][i] = math.sqrt(valor_raiz)
        #PASO 2: Calculamos el resto de los elementos en la fila i (a la derecha)
        for j in range(i + 1, n):
            suma_cruzada = 0.0
            for k in range(i):
                #Multiplicando el elemento de la columna i por el de la columna j
                suma_cruzada += G[k][i] * G[k][j]    
            G[i][j] = (A[i][j] - suma_cruzada) / G[i][i]
            
    return G

def transponer_matriz(G):
    n = len(G)
    #Creando una nueva matriz llena de ceros
    G_T = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            #Asignamos invirtiendo los indices
            G_T[i][j] = G[j][i]
    return G_T

def resolver_cholesky(A, b):
    #Obtenemos G usando nuestro algoritmo de Cholesky
    G = cholesky(A)
    #Calculamos la transpuesta G^T
    G_T = transponer_matriz(G)
    #Resolvemos (G^T) * d = b
    d = sustitucion_forward(G_T, b)
    #Resolvemos G * x = d
    x = sustitucion_backward(G, d)
    
    return x

def factorizacion_lu_crout(A):
    n = len(A)
    #Creamos matrices L y U llenas de ceros
    L = [[0.0] * n for _ in range(n)]
    U = [[0.0] * n for _ in range(n)]
    #llenamos la diagonal de L con unos
    for i in range(n):
        L[i][i] = 1.0
    #Recorremos cada columna j, para resolver lA sumatoria
    for j in range(n):
        #PASO 1: Calculando los elementos de U (filas i <= j)
        for i in range(j + 1):
            suma = 0.0
            #La sumatoria va desde k = 0 hasta i - 1
            for k in range(i):
                suma += L[i][k] * U[k][j]    
            U[i][j] = A[i][j] - suma
        #Comprobación de seguridad: Si el pivote en U es cero, el método falla
        if abs(U[j][j]) < 1e-15:
            raise ValueError(f"Pivote cero detectado en U[{j}][{j}].")
        #PASO 2: Calculando los elementos de L (filas i > j)
        for i in range(j + 1, n):
            suma = 0.0
            # La sumatoria va desde k = 0 hasta j - 1
            for k in range(j):
                suma += L[i][k] * U[k][j]
            L[i][j] = (A[i][j] - suma) / U[j][j]
    #Finalemnte retornamos las matrices L y U        
    return L, U

def resolver_lu_simple(A, b):
    #Factorizamos la matriz A
    L, U = factorizacion_lu_crout(A)
    #Resolvemos L * d = b (Sustitución Forward)
    # Como L es triangular inferior, usamos sustitución hacia adelante
    d = sustitucion_forward(L, b)
    #Resolviendo U * x = d (Sustitución Backward)
    # Como U es triangular superior, usamos sustitución hacia atrás
    x = sustitucion_backward(U, d)
    
    return x

def es_simetrica(A, tol=1e-12):
    n = len(A)

    for i in range(n):
        for j in range(i + 1, n):
            if abs(A[i][j] - A[j][i]) > tol:
                return False

    return True

def matriz_por_vector(A, x):
    n = len(x)
    resultado = np.zeros(n, dtype = float)
    for i in range(n):
        suma = 0.0
        for j in range(n):
            suma = suma + A[i][j]*x[j]
        resultado[i] = suma
    return resultado

def residuo(A, x, b):
    n = len(b)
    Ax = matriz_por_vector(A, x)
    r = np.zeros(n, dtype=float)
    for i in range(n):
        r[i] = Ax[i] - b[i]
    return r

def calcular_delta(x_nuevo, x_anterior):
    n = len(x_nuevo)
    delta = 0.0
    for i in range(n):
        diferencia = abs(x_nuevo[i] - x_anterior[i])
        if diferencia > delta:
            delta = diferencia
    return delta

def sor(A, b, omega =1.2, tol= 1e-8, max_iter=1000):
    n = len(b)
    x = np.zeros(n, dtype=float)
    #Verificando el parametro w
    if omega<=0 or omega >=2:
        raise ValueError("Omega debe estar entre 0 y 2")
    #verificando la diagonal
    for i in range(n):
            if abs(A[i][i]) < 1e-15:
                raise ValueError("Elemento diagonal nulo")
    #Comenzando las iteraciones
    for k in range(max_iter):
        # Guardamos los valores anteriores
        x_anterior = x.copy()
        for i in range (n):
            suma = 0.0
            for j in range(n):
                if j != i:
                    suma = suma + A[i][j]*x[j]
            #Calculamos el valor de Gauss-Seidel
            x_gs = (b[i] - suma)/A[i][i]
            #Aplicamos la relajacion
            x[i] = (1 - omega)*x[i] + omega*x_gs
        # Calculamos delta
        delta = calcular_delta(x, x_anterior)
        # Calculando el residuo
        r = residuo(A, x, b)
        #comprobando la convergencia
        if delta < tol and norma_vec(r, 'inf') < tol:
            return x, k + 1, delta

    raise ValueError("SOR no alcanzo la convergencia")

def gauss_seidel(A, b, tol=1e-8, max_iter=1000):
    return sor(A, b, omega=1.0, tol=tol, max_iter=max_iter)

def sustitucion_backward(U, b):
    ##Resolviendo Ux = b para U triangular superior
    n = len(b)
    x = np.zeros(n, dtype=float)
    for k in range(n - 1, -1, -1):
        if abs(U[k][k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = 0
        for i in range (k+1, n):
            suma = suma + U[k][i] * x[i]
        x[k] = (b[k] - suma) / U[k][k]
    return x

def gauss_pivoteo (A, b):
    n = len(b)
    A_copy = [fila[:] for fila in A]
    b_copy = b[:]

    #Etapa 1: Eliminacion hacia delante con pivoteo parcial
    for k in range(n-1):
        #Asumiendo inicialmente que el valor maximo esta en la fila actual k
        p = k 
        max_valor = abs(A_copy[k][k])
        #Buscando si hay un valor mas grande debajo en la misma columna k
        for i in range(k+1, n):
            if abs(A_copy[i][k]) > max_valor:
                max_valor = abs(A_copy[i][k])
                p = i
        #Si el maximo no esta en la fila k, intercambiamos la fila k con la fila p
        if p != k:
            #Intercambiamos de la fila completa en A
            for j in range(k, n): #Solo necesitamos intercambiar desde la columna k en adelante
                temp = A_copy[k][j]
                A_copy[k][j] = A_copy[p][j]
                A_copy[p][j] = temp
            #Intercambiamos la columna b tambien
            temp_b = b_copy[k]
            b_copy[k] = b_copy[p]
            b_copy[p] = temp_b
        #Verificando si es que toda la columan es cero
        if abs(A_copy[k][k]) < 1e-15:
            raise ValueError("El sistema no tiene solucion única")
        #Aplicando la eliminación gausiana
        for i in range(k+1, n):
            factor = A_copy[i][k]/A_copy[k][k]
            for j in range(k+1, n):
                A_copy[i][j] = A_copy[i][j] - factor*A_copy[k][j]
            b_copy[i] = b_copy[i] - factor*b_copy[k]
    #Finalemnte usamos el metodo sustitucion backward
    return sustitucion_backward(A_copy, b_copy)

def norma_vec(v, tipo=2):
    #Calculando la norma vectorial 1, 2 o infinito ('inf')
    v = np.asarray(v, dtype=float)
    if tipo == 1:
        return float(np.sum(np.abs(v)))
    elif tipo == 2:
        return float(np.sqrt(np.sum(v ** 2)))
    elif tipo == 'inf':
        return float(np.max(np.abs(v)))
    else:
        raise ValueError("Tipo de norma no soportado. Use 1, 2 o 'inf'.")

def calcular_inversa(A):
    n = len(A)
    #Creando una matriz vacia para guardar la inversa
    A_inv = [[0.0] * n for _ in range(n)]
    #Resolvemos el sistema n veces (una por cada columna)
    for j in range(n):
        #Construimos la columna j de la matriz Identidad (vector I_j)
        I_j = [0.0] * n
        I_j[j] = 1.0
        #Usando nuestro metodo de Gauss para resolver A * x = I_j
        columna_solucion = gauss_pivoteo(A, I_j)
        #Guardamos esa solución como la columna j de nuestra matriz inversa
        for i in range(n):
            A_inv[i][j] = columna_solucion[i]
            
    return A_inv

def norma_mat(A, tipo='F'):
    A = np.asarray(A, dtype=float)
    if tipo == 'F':
        #Suma el cuadrado de absolutamente todos los elementos y saca la raiz
        return float(np.sqrt(np.sum(A ** 2)))   
    elif tipo == 'inf':
        #1. np.abs(A) vuelve todo positivo.
        #2. np.sum(..., axis=1) suma cada fila individualmente, creando un vector de sumas.
        #3. np.max() saca el valor más grande de ese vector de sumas.
        return float(np.max(np.sum(np.abs(A), axis=1)))  
    elif tipo == 1:
        #Igual, pero axis=0 suma verticalmente (cada columna).
        return float(np.max(np.sum(np.abs(A), axis=0)))    
    else:
        raise ValueError("Tipo de norma no soportado. Use 'F', 'inf' o 1.") 

def numero_condicion(A, tipo='F'):
    #Calculando kappa(A) = ||A|| * ||A^{-1}|| en la norma indicada
    A_inv = calcular_inversa(A)
    return norma_mat(A, tipo=tipo) * norma_mat(A_inv, tipo=tipo)

def es_diagonal_dominante(A):
    n = len(A)
    dominante = True
    for i in range(n):
        suma = 0.0
        for j in range(n):
            if i != j:
                suma = suma + abs(A[i][j])
        if abs(A[i][i]) <= suma:
            dominante = False
    return dominante

def solucionador_automatico(A, b):
    condicion = numero_condicion(A)
    
    if condicion > 10000:
        return gauss_pivoteo(A, b)
        
    simetrica = es_simetrica(A)
    
    if simetrica:
        es_positiva = True
        try:
            cholesky(A)
        except:
            es_positiva = False
            
        if es_positiva:
            return resolver_cholesky(A, b)
            
    dominante = es_diagonal_dominante(A)
    
    if dominante:
        return gauss_seidel(A, b)
        
    return resolver_lu_simple(A, b)

class InterpolacionMonomios:
    def __init__(self, x_datos, y_datos):
        self.x = x_datos
        self.y = y_datos
        self.n = len(x_datos)

        #Construyendo la matriz de Vandermonde
        self.A = [[0.0]* self.n for _ in range (self.n)]
        for i in range(self.n):
            for j in range(self.n):
                self.A[i][j] = (self.x[i]**j)

        #Resolviendo el sistema A*lambda = y
        self.lambdas = solucionador_automatico(self.A, self.y)
        
    def evaluar(self, t):
        resultado = 0.0
        for j in range(self.n):
            resultado = resultado + self.lambdas[j] * (t ** j)
        return resultado
    
    def graficar(self, x_min=None, x_max=None):
        if x_min is None and x_max is None:
            margen = (max(self.x) - min(self.x)) * 0.1
            x_min = min(self.x) - margen
            x_max = max(self.x) + margen
        else:
            if x_min is None:
                x_min = min(self.x)
            if x_max is None:
                x_max = max(self.x)

        t_vector = np.linspace(x_min, x_max, 100)
        y_vector = np.zeros(100)
        for i in range(100):
            y_vector[i] = self.evaluar(t_vector[i])
            
        plt.figure(figsize=(8, 5))
        plt.plot(t_vector, y_vector, 'b-', label="Polinomio interpolante P5")
        plt.plot(self.x, self.y, 'ro', markersize=8, label="Datos")
        plt.xlim(x_min, x_max)
        plt.xlabel("Amplitud (grados)")
        plt.ylabel("T / T0")
        plt.title("Período del péndulo de gran amplitud")
        plt.grid(True)
        plt.legend()
        plt.show()

class InterpolacionLagrange:
    def __init__(self, x_datos, y_datos):
        self.x = x_datos
        self.y = y_datos
        self.n = len(x_datos)
        self.w = []
        
        log_pesos = []
        signos = []
        for j in range(self.n):
            suma = 0.0
            signo = 1.0
            for k in range(self.n):
                if j != k:
                    dif = self.x[j] - self.x[k]
                    if dif < 0.0:
                        signo = signo * -1.0
                    suma = suma + np.log(abs(dif))
            log_pesos.append(-suma)
            signos.append(signo)
            
        max_log = max(log_pesos)
        for j in range(self.n):
            peso_real = signos[j] * np.exp(log_pesos[j] - max_log)
            self.w.append(peso_real)
        

    def evaluar(self, t):
        for i in range(self.n):
            if abs(t - self.x[i]) < 1e-12:
                return self.y[i]
        
        num = 0.0
        den = 0.0
        for j in range(self.n):
            termino = self.w[j] / (t - self.x[j])
            num = num + termino * self.y[j]
            den = den + termino
            
        return num / den
    
    def graficar(self):
        x_min = min(self.x)
        x_max = max(self.x)
        margen = (x_max - x_min) * 0.1
        t_vector = np.linspace(x_min - margen, x_max + margen, 100)
        y_vector = np.zeros(100)
        for i in range(100):
            y_vector[i] = self.evaluar(t_vector[i])
            
        plt.figure(figsize=(8, 5))
        plt.plot(t_vector, y_vector, color='orange', linestyle='-', label="Lagrange")
        plt.plot(self.x, self.y, 'ro', markersize=8)
        plt.grid(True)
        plt.legend()
        plt.show()


def lagrange_k7(x_datos, y_datos, omega):
    interpolante = InterpolacionLagrange(x_datos, y_datos)
    return interpolante.evaluar(omega)


class RegresionLineal:
    def __init__(self, x_datos, y_datos):
        self.x = x_datos
        self.y = y_datos
        self.n = len(x_datos)

        self.a0 = 0.0
        self.a1 = 0.0
        self.r2 = 0.0
        self.syx = 0.0

        self._calcularAjuste()

    def _calcularAjuste(self):
        suma_x = 0.0
        suma_y = 0.0
        suma_xy = 0.0
        suma_x2 = 0.0

        for i in range(self.n):
            suma_x = suma_x + self.x[i]
            suma_y = suma_y + self.y[i]
            suma_xy = suma_xy + (self.x[i] * self.y[i])
            suma_x2 = suma_x2 + (self.x[i] * self.x[i])

        promedio_x = suma_x / self.n
        promedio_y = suma_y / self.n

        #Aplicando las formulas para hallar los coeficientes
        numerador_a1 = (self.n * suma_xy) - (suma_x * suma_y)
        denominador_a1 = (self.n * suma_x2) - (suma_x * suma_x)
        self.a1 = numerador_a1 / denominador_a1
        self.a0 = promedio_y - (self.a1 * promedio_x)

        #Calculando los errores
        st = 0.0
        sr = 0.0
        for i in range(self.n):
            st = st + (self.y[i] - promedio_y)**2
            modelo_y = self.a0 + self.a1 * self.x[i]
            sr = sr + (self.y[i] - modelo_y)**2

        self.syx = np.sqrt(sr /(self.n - 2))
        self.r2 = (st - sr) / st

    def evaluar(self, x_val):
        return self.a0 + self.a1 * x_val
    
class RegresionPolinomial:
    def __init__(self, x_datos, y_datos, grado):
        self.x = x_datos
        self.y = y_datos
        self.n = len(x_datos)
        self.m = grado
        self.coeficientes = []
        self.syx = 0.0
        self.r2 = 0.0
        self._calcularAjuste()

    def _calcularAjuste(self):
        filas = self.m + 1
        matriz = []
        for i in range(filas):
            fila_actual = []
            for j in range(filas):
                fila_actual.append(0.0)
            matriz.append(fila_actual)

        vector_b = []
        for i in range(filas):
            vector_b.append(0.0)

        for i in range(filas):
            for j in range(i + 1):
                k = i + j
                suma = 0.0
                for l in range(self.n):
                    suma = suma + (self.x[l] ** k)
                matriz[i][j] = suma
                matriz[j][i] = suma

            suma_y = 0.0
            for l in range(self.n):
                suma_y = suma_y + self.y[l] * (self.x[l] ** i)
            vector_b[i] = suma_y

        self.coeficientes = solucionador_automatico(matriz, vector_b)

        st = 0.0
        sr = 0.0
        suma_y_total = 0.0
        for i in range(self.n):
            suma_y_total = suma_y_total + self.y[i]
        promedio_y = suma_y_total / self.n

        for i in range(self.n):
            st = st + (self.y[i] - promedio_y)**2
            estimado = 0.0
            for j in range(filas):
                estimado = estimado + self.coeficientes[j] * (self.x[i]**j)
            sr = sr + (self.y[i] - estimado)**2

        self.syx = np.sqrt(sr / (self.n - (self.m + 1)))
        self.r2 = (st - sr) / st

    def evaluar(self, x_val):
        resultado = 0.0
        for j in range(self.m + 1):
            resultado = resultado + self.coeficientes[j] * (x_val ** j)
        return resultado

def ejercicio1():

    theta_datos = [0, 30, 60, 90, 120, 150]
    razones_periodo = [1.00000, 1.01741, 1.07318, 1.18034, 1.37288, 1.76220]
    interpolante = InterpolacionMonomios(theta_datos, razones_periodo)

    print("\n Item(a):")
    for i, coeficiente in enumerate(interpolante.lambdas):
        print(f"a{i} = {coeficiente:.12e}")

    print("\n Item(b):")
    valores_exactos = [(45, 1.03997), (100, 1.23223)]
    for theta, valor_exacto in valores_exactos:
        valor_aproximado = interpolante.evaluar(theta)
        error_relativo = abs(valor_aproximado - valor_exacto) / abs(valor_exacto)
        print(f"P5({theta} grados) = {valor_aproximado:.8f} ")
        print(f"valor exacto = {valor_exacto:.5f}")
        print(f"error relativo = {error_relativo:.6%}")

    print("\n Item(c):")
    valores_extrapolacion = [(165, 2.18544), (175, 2.87766)]
    for theta, valor_exacto in valores_extrapolacion:
        valor_aproximado = interpolante.evaluar(theta)
        error_relativo = abs(valor_aproximado - valor_exacto) / abs(valor_exacto)
        print(f"P5({theta} grados) = {valor_aproximado:.8f}")
        print(f"valor exacto = {valor_exacto:.5f} ")
        print(f"error relativo = {error_relativo:.6%}")


    print("\n Item(d):")
    print("Cuando la amplitud se acerca a 180 grados, el péndulo pasa cada vez ")
    print("más tiempo cerca de la posición vertical invertida, que es inestable ")
    print("y donde su velocidad angular se aproxima a cero. Por ello, el período ")
    print("tiende a infinito. Un polinomio permanece finito en 180 grados y no ")
    print("puede representar correctamente esa divergencia.")
    
    print("Se grafica P5 junto con los datos en el intervalo [0, 175] grados.")
    interpolante.graficar(0, 175)
    
def ejercicio2():
    print("\n Item(a):")

    def amplitud(omega):
        denominador = ((1 - omega**2) ** 2 + 0.16 * omega**2) ** 0.5
        return 1 / denominador

    omega_grafica = np.linspace(0, 2, 401)
    amplitud_exacta = []
    for omega in omega_grafica:
        amplitud_exacta.append(amplitud(omega))

    plt.figure(figsize=(9, 6))
    plt.plot(omega_grafica, amplitud_exacta, "k-", linewidth=2, label="A(omega)")

    grados = [4, 8, 12]
    colores = ["b", "r", "g"]
    for indice in range(len(grados)):
        n = grados[indice]
        nodos = []
        for i in range(n + 1):
            nodos.append(2 * i / n)

        amplitudes_nodos = []
        for omega in nodos:
            amplitudes_nodos.append(amplitud(omega))

        polinomio = []
        for omega in omega_grafica:
            polinomio.append(lagrange_k7(nodos, amplitudes_nodos, omega))

        plt.plot(omega_grafica, polinomio, color=colores[indice], label=f"P{n}(omega), n = {n}",)
        plt.plot(nodos, amplitudes_nodos, marker="o", linestyle="None", color=colores[indice], markersize=3,)

    plt.xlabel("Frecuencia omega")
    plt.ylabel("Amplitud A(omega)")
    plt.title("Interpolación de Lagrange de la amplitud")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
    
    print("\n Item(b):")
    def puntosChebyshev(n):
        puntos = []
        for j in range(n+1):
            angulo = (((2*j) + 1) * np.pi)/(2*(n + 1))
            x_j = (1.0 + np.cos(angulo))
            puntos.append(x_j)
        return puntos

    plt.figure(figsize=(9, 6))
    plt.plot(omega_grafica, amplitud_exacta, "k-", linewidth=2, label="A(omega)")

    for indice in range(len(grados)):
        n = grados[indice]
        nodos = puntosChebyshev(n)

        amplitudes_nodos = []
        for omega in nodos:
            amplitudes_nodos.append(amplitud(omega))

        polinomio = []
        for omega in omega_grafica:
            polinomio.append(lagrange_k7(nodos, amplitudes_nodos, omega))

        plt.plot(omega_grafica, polinomio,color=colores[indice], label=f"P{n}(omega), n = {n} (Chebyshev)",)
        plt.plot(nodos, amplitudes_nodos, marker="o", linestyle="None", color=colores[indice], markersize=3,)

    plt.xlabel("Frecuencia omega")
    plt.ylabel("Amplitud A(omega)")
    plt.title("Interpolación de Lagrange con nodos de Chebyshev")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

def ejercicio3():
    print("\n Item(a):")
    a_datos = [0.387, 0.723, 1.000, 1.524, 5.203, 9.537, 19.19, 30.07]
    T_datos = [0.241, 0.615, 1.000, 1.881, 11.862, 29.457, 84.01, 164.8]
    n = len(a_datos)
    
    ln_a = []
    ln_T = []
    for i in range(n):
        ln_a.append(math.log(a_datos[i]))
        ln_T.append(math.log(T_datos[i]))
        
    modelo1 = RegresionLineal(ln_a, ln_T)
    m = modelo1.a1
    r2 = modelo1.r2
    print("Potencia m:", m)
    print("R2:", r2)
    print("Como el r^2 es muy cercano a 1, esto comprueba")
    print("que el modelo1 es muy compatible con la tercera ley de kepler,")
    print("además que el valor del exponenente es muy cercado a 1.5")
    
    print("\n Item(b):")
    t_datos = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
    y_datos = [50.00, 49.76, 47.09, 42.00, 34.42, 24.46, 12.00] 
    n = len(t_datos)
    
    modelo2 = RegresionPolinomial(t_datos, y_datos, 2)
    coeficientes = modelo2.coeficientes
    print("y_0 (m):", coeficientes[0])
    print("v_0(m/s):", coeficientes[1])
    print("Aceleracion de la gravedad g(m/s^2):", (-2)*coeficientes[2])
    
    print("\n Item(c):")
    def funcion(x):
        return coeficientes[0] + coeficientes[1] * x + coeficientes[2]*(x**2)
    
    ts_aproximado, errorests, iter = biseccion(funcion, 2, 4, 1e-6)
    print("t_s(s):", ts_aproximado)
    velocidad_impacto = coeficientes[1] + 2 * coeficientes[2] * ts_aproximado
    print("Velocidad de impacto (m/s):", velocidad_impacto)

    t_curva = np.linspace(0, ts_aproximado, 100)
    y_curva = []
    for i in range(len(t_curva)):
        y_curva.append(funcion(t_curva[i]))

    plt.figure(figsize=(8, 5))
    plt.plot(t_datos, y_datos, "ro", markersize=7, label="Datos")
    plt.plot(t_curva, y_curva, "b-", label="Ajuste cuadrático")
    plt.plot(ts_aproximado, 0, "kx", markersize=10, label="Impacto")
    plt.xlabel("Tiempo (s)")
    plt.ylabel("Altura y (m)")
    plt.title("Altura de la esfera y punto de impacto")
    plt.grid(True)
    plt.legend()
    plt.show()

        
ejercicio1()
ejercicio2()
ejercicio3()