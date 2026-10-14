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
    """Resuelve Lx = b para L triangular inferior."""
    n = len(b)
    x = np.zeros(n, dtype=float)
    for k in range(n):
        if abs(L[k, k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = np.dot(L[k, :k], x[:k])
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
        suma = np.dot(U[k, k+1:], x[k+1:])
        x[k] = (b[k] - suma) / U[k, k]
    return x

#Método de gauss pivote
def gauss_pivoteo(A, b, usar_pivoteo=True):
    """
    Resuelve Ax = b mediante Eliminación de Gauss (con o sin pivoteo parcial).
    Trabaja sobre copias para no alterar las matrices originales.
    """
    A_copia = np.array(A, dtype=float, copy=True)
    b_copia = np.array(b, dtype=float, copy=True)
    n = len(b_copia)

    # 1. Eliminación hacia adelante
    for k in range(n - 1):
        if usar_pivoteo:
            # Buscar el máximo en valor absoluto en la columna k desde la fila k
            p = k + np.argmax(np.abs(A_copia[k:, k]))
            if p != k:
                A_copia[[k, p], :] = A_copia[[p, k], :]
                b_copia[[k, p]] = b_copia[[p, k]]

        if abs(A_copia[k, k]) < 1e-15:
            raise ValueError("Matriz singular o requiere pivoteo.")

        for i in range(k + 1, n):
            factor = A_copia[i, k] / A_copia[k, k]
            A_copia[i, k:] -= factor * A_copia[k, k:]
            b_copia[i] -= factor * b_copia[k]

    # 2. Sustitución inversa
    return sustitucion_backward(A_copia, b_copia)

#Factorizacion LU
def factorizacion_lu_crout(A):
    """
    Descompone A = L * U según el algoritmo de Crout (L con 1s en la diagonal).
    Retorna (L, U).
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    L = np.eye(n, dtype=float)
    U = np.zeros((n, n), dtype=float)

    for j in range(n):
        # Procedimiento 1: Calcular elementos de U (i <= j)
        for i in range(j + 1):
            suma = np.dot(L[i, :i], U[:i, j])
            U[i, j] = A[i, j] - suma

        if abs(U[j, j]) < 1e-15:
            raise ValueError("Pivote nulo encontrado en factorización LU.")

        # Procedimiento 2: Calcular elementos de L (i > j)
        for i in range(j + 1, n):
            suma = np.dot(L[i, :j], U[:j, j])
            L[i, j] = (A[i, j] - suma) / U[j, j]

    return L, U

def resolver_lu(A, b):
    """Resuelve Ax = b usando factorización LU de Crout."""
    L, U = factorizacion_lu_crout(A)
    d = sustitucion_forward(L, b)
    x = sustitucion_backward(U, d)
    return x

#método cholesky
def factorizacion_cholesky(A):
    """
    Descompone una matriz SPD como A = G^T * G, donde G es triangular superior.
    Retorna G.
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    G = np.zeros((n, n), dtype=float)

    for i in range(n):
        suma_diag = np.dot(G[:i, i], G[:i, i])
        val = A[i, i] - suma_diag
        if val <= 0:
            raise ValueError("La matriz no es positiva definida.")
        G[i, i] = np.sqrt(val)

        for j in range(i + 1, n):
            suma = np.dot(G[:i, i], G[:i, j])
            G[i, j] = (A[i, j] - suma) / G[i, i]

    return G

def resolver_cholesky(A, b):
    """Resuelve Ax = b (donde A = G^T G) usando Cholesky."""
    G = factorizacion_cholesky(A)
    d = sustitucion_forward(G.T, b)
    x = sustitucion_backward(G, d)
    return x

# 2. MÉTODOS ITERATIVOS ESTACIONARIOS

#metodo jacobi
def jacobi(A, b, x0=None, tol=1e-6, max_iter=1000):
    """
    Método de Jacobi.
    Retorna: (x, iteraciones, historial_errores)
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n, dtype=float) if x0 is None else np.array(x0, dtype=float, copy=True)
    historial_err = []

    for k in range(1, max_iter + 1):
        x_new = np.zeros(n, dtype=float)
        for i in range(n):
            # Suma de a_ij * x_j para j != i
            suma = np.dot(A[i, :i], x[:i]) + np.dot(A[i, i+1:], x[i+1:])
            x_new[i] = (b[i] - suma) / A[i, i]

        delta = norma_vec(x_new - x, tipo=2)
        historial_err.append(delta)
        x = x_new

        if delta < tol:
            return x, k, historial_err

    return x, max_iter, historial_err

#metodo sor
def sor(A, b, omega=1.2, x0=None, tol=1e-6, max_iter=1000):
    """
    Método Successive Over-Relaxation (SOR).
    Si omega = 1.0, es exactamente Gauss-Seidel.
    Retorna: (x, iteraciones, historial_errores)
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n, dtype=float) if x0 is None else np.array(x0, dtype=float, copy=True)
    historial_err = []

    for k in range(1, max_iter + 1):
        x_old = x.copy()
        for i in range(n):
            # Usa x[:i] ya actualizados (k+1) y x[i+1:] antiguos (k)
            suma_nuevos = np.dot(A[i, :i], x[:i])
            suma_antiguos = np.dot(A[i, i+1:], x_old[i+1:])
            x_gs = (b[i] - suma_nuevos - suma_antiguos) / A[i, i]
            x[i] = (1.0 - omega) * x_old[i] + omega * x_gs

        delta = norma_vec(x - x_old, tipo=2)
        historial_err.append(delta)

        if delta < tol:
            return x, k, historial_err

    return x, max_iter, historial_err

#metodo gauss seidel
def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=1000):
    """
    Método de Gauss-Seidel (caso particular de SOR con omega = 1.0).
    Retorna: (x, iteraciones, historial_errores)
    """
    return sor(A, b, omega=1.0, x0=x0, tol=tol, max_iter=max_iter)

# 3. MÉTODOS DE OPTIMIZACIÓN (PARA MATRICES SPD)

#Maximo descenso
def maximo_descenso(A, b, x0=None, tol=1e-6, max_iter=1000):
    """
    Método de Máximo Descenso (Steepest Descent).
    Sigue exactamente la convención r_k = A*x_k - b de las diapositivas.
    Retorna: (x, iteraciones, historial_errores)
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n, dtype=float) if x0 is None else np.array(x0, dtype=float, copy=True)
    historial_err = []

    for k in range(1, max_iter + 1):
        r = A @ x - b  # Gradiente: nabla Phi(x_k)
        rr = np.dot(r, r)
        if np.sqrt(rr) < tol:
            break

        Ar = A @ r
        alpha = rr / np.dot(r, Ar)
        x_new = x - alpha * r

        delta = norma_vec(x_new - x, tipo=2)
        historial_err.append(delta)
        x = x_new

        if delta < tol:
            return x, k, historial_err

    return x, max_iter, historial_err

#gradiente conjugado
def gradiente_conjugado(A, b, x0=None, tol=1e-6, max_iter=1000):
    """
    Método de Gradiente Conjugado.
    Sigue exactamente el pseudocódigo de las diapositivas 18 y 19:
    r_0 = A*x_0 - b, p_0 = -r_0.
    Retorna: (x, iteraciones, historial_errores)
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    x = np.zeros(n, dtype=float) if x0 is None else np.array(x0, dtype=float, copy=True)

    r = A @ x - b
    p = -r.copy()
    rr = np.dot(r, r)
    historial_err = []

    if np.sqrt(rr) < tol:
        return x, 0, [0.0]

    for k in range(1, max_iter + 1):
        Ap = A @ p
        pAp = np.dot(p, Ap)
        if abs(pAp) < 1e-16:
            break

        alpha = rr / pAp
        x_new = x + alpha * p
        r_new = r + alpha * Ap

        delta = norma_vec(x_new - x, tipo=2)
        historial_err.append(delta)

        rr_new = np.dot(r_new, r_new)
        if delta < tol or np.sqrt(rr_new) < tol:
            return x_new, k, historial_err

        beta = rr_new / rr
        p = -r_new + beta * p

        x = x_new
        r = r_new
        rr = rr_new

    return x, max_iter, historial_err

# 4. INVERSA, NORMAS Y CONDICIONAMIENTO (SIN NUMPY.LINALG)

#inversa gauss
def inversa_gauss(A):
    """
    Calcula A^{-1} sin numpy.linalg resolviendo A * x_j = e_j para cada columna j,
    tal como solicita el problema de la densidad de la varilla.
    """
    A = np.array(A, dtype=float)
    n = A.shape[0]
    A_inv = np.zeros((n, n), dtype=float)
    I = np.eye(n, dtype=float)

    for j in range(n):
        A_inv[:, j] = gauss_pivoteo(A, I[:, j], usar_pivoteo=True)

    return A_inv

#inversa norma vectorial
def norma_vec(v, tipo=2):
    """Calcula la norma vectorial 1, 2 o infinito ('inf') sin numpy.linalg."""
    v = np.asarray(v, dtype=float)
    if tipo == 1:
        return float(np.sum(np.abs(v)))
    elif tipo == 2:
        return float(np.sqrt(np.sum(v ** 2)))
    elif tipo == 'inf':
        return float(np.max(np.abs(v)))
    else:
        raise ValueError("Tipo de norma no soportado. Use 1, 2 o 'inf'.")
    
#Calcula la norma matricial
def norma_mat(A, tipo='F'):
    """
    Calcula la norma matricial sin numpy.linalg:
    - 'F': Norma de Frobenius
    - 'inf': Máxima suma absoluta de filas
    - 1: Máxima suma absoluta de columnas
    """
    A = np.asarray(A, dtype=float)
    if tipo == 'F':
        return float(np.sqrt(np.sum(A ** 2)))
    elif tipo == 'inf':
        return float(np.max(np.sum(np.abs(A), axis=1)))
    elif tipo == 1:
        return float(np.max(np.sum(np.abs(A), axis=0)))
    else:
        raise ValueError("Use 'F', 'inf' o 1.")
    
#Numero de condiciones
def numero_condicion(A, tipo='F'):
    """Calcula kappa(A) = ||A|| * ||A^{-1}|| en la norma indicada."""
    A_inv = inversa_gauss(A)
    return norma_mat(A, tipo=tipo) * norma_mat(A_inv, tipo=tipo)