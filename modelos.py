import numpy as np
#Formula para el mapa logisitico
def X(n, x0, r):
    Datos = [x0]
        
    for i in range(n):
        siguiente = r * (Datos[-1] * (1 - Datos[-1]))
        Datos.append(siguiente)
    return Datos

#Formula para calcular el factorial de un numero
def factorial(n: int) -> int:
    resultado = 1
    # Multiplica secuencialmente desde 2 hasta n
    for i in range(2, n + 1):
        resultado *= i
    return resultado

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

#Metodo de la falsa posicion
def falsaPosicion(funcion, xl:float, xu:float, tolerancia, max_iter = 100):
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
        #para reducir calculo inecesario
        fl = funcion(xl)
        fu = funcion(xu)

        # Nueva fórmula para hallar xr (la recta secante)
        xr = xu - (fu * (xl - xu)) / (fl - fu)
        
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

#Metodo del punto fijo
def puntoFijo(funcion, x0:float , tolerancia ,max_iter = 100):
    #Variables a utilizar
    iteracion:int = 0
    ea:float = 1.0
    xr:float = x0
    errores_historial = []

    #implementando un bucle para hacer el metodo
    while ea > tolerancia and iteracion < max_iter:
        #Guardando el primer punto para luego medir el error
        xr_old = xr 

        #Asignando el punto siguiente mediante la funcion
        xr = funcion(xr_old)

        #Calculando el error desde la primera vuelta
        if xr != 0: #Esto lo hacemos para cuidarnos de la division por cero
            ea = abs((xr - xr_old)/ xr)
            errores_historial.append(ea)

        iteracion += 1

    return xr, errores_historial, iteracion

#Metodo de Newton-Raphson
def newtonRaphson(funcion, Dfuncion , x0:float , tolerancia ,max_iter = 100):
    #Variables a utilizar
    iteracion:int = 0
    ea:float = 1.0
    xr:float = x0
    errores_historial = []

    #implementando un bucle para hacer el metodo
    while ea > tolerancia and iteracion < max_iter:
        #Guardando el primer punto para luego medir el error
        xr_old = xr 
        df = Dfuncion(xr_old)
        
        #Verificando que la derivada no sea 0 para evitar errores
        if df == 0:
            print("Error: La derivada es cero. El método ha colapsado.")
            return None, None, None  #Abortamos sin usar break

        #Asignando el punto siguiente mediante la regla
        xr = xr_old - ((funcion(xr_old))/(df))

        #Calculando el error desde la primera vuelta
        if xr != 0: #Esto lo hacemos para cuidarnos de la division por cero
            ea = abs((xr - xr_old)/ xr)
            errores_historial.append(ea)

        iteracion += 1

    return xr, errores_historial, iteracion

#Metodo de la secante
#IDEA: Es similar al metodo de la falsa posicion, solo que utilizando un metodo abierto
def secante(funcion, x0:float, x1:float, tolerancia, max_iter = 100):
    iteracion: int = 0
    ea: float = 1.0
    errores_historial = []
    
    while ea > tolerancia and iteracion < max_iter:
        f0 = funcion(x0)
        f1 = funcion(x1)
        
        #Evitando la división por cero 
        if (f0 - f1) == 0.0:
            print("Error: División por cero. Pendiente nula.")
            return None, None, None

        # Fórmula de la secante
        x2 = x1 - (f1 * (x0 - x1)) / (f0 - f1)

        # Calculando el error
        if x2 != 0:
            ea = abs((x2 - x1) / x2)
            errores_historial.append(ea)

        # Deslizando los puntos para la siguiente vuelta
        x0 = x1
        x1 = x2
        
        iteracion += 1

    return x2, errores_historial, iteracion


#Métodos para solucionar sistemas lineales

#Sustitución forward
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

#sustitucion backward
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

#Metodo de gauss con pivote
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

#metodo de factorizacion LU (crout)
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

#Ahora creamos la funcion para resolver el LU
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

#Haciendo el metodo de cholesky
import math

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

#Codigo para la transpuesta de una matriz
def transponer_matriz(G):
    n = len(G)
    #Creando una nueva matriz llena de ceros
    G_T = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            #Asignamos invirtiendo los indices
            G_T[i][j] = G[j][i]
    return G_T

#Resolviendo cholesky
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

#Algoritmo para calcular la inversa de una matriz
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

#metodo para la norma de un vector
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

#Metodo para la norma de una matriz
#axis=1: Suma horizontalmente (por filas)
#axis=0: Suma verticalmente (por columnas)
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

#Calculando el numero de condicion
def numero_condicion(A, tipo='F'):
    #Calculando kappa(A) = ||A|| * ||A^{-1}|| en la norma indicada
    A_inv = calcular_inversa(A)
    return norma_mat(A, tipo=tipo) * norma_mat(A_inv, tipo=tipo)