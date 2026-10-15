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
    #Resuelve Lx = b para L triangular inferior
    n = len(b)
    x = np.zeros(n, dtype=float)
    for k in range(n):
        if abs(L[k, k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = 0
        for i in range(k):
            suma = suma + L[k , i]*x[i]
        x[k] = (b[k] - suma) / L[k, k]
    return x

#sustitucion backward
def sustitucion_backward(U, b):
    """Resuelve Ux = b para U triangular superior."""
    n = len(b)
    x = np.zeros(n, dtype=float)
    for k in range(n - 1, -1, -1):
        if abs(U[k, k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = 0
        for i in range (K+1, n):
            suma = suma + U[k,i] * x[i]
        x[k] = (b[k] - suma) / U[k, k]
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