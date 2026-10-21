import math
import random
import numpy as np
import matplotlib.pyplot as plt

##########################################
# FUNCIONES AUXILIARES

# Norma de un vector
def norma_vec(v, tipo=2):
    if tipo == 1:
        suma = 0.0
        for i in range(len(v)):
            suma += abs(v[i])
        return suma

    elif tipo == 2:
        suma = 0.0
        for i in range(len(v)):
            suma += v[i]*v[i]
        return math.sqrt(suma)

    elif tipo == 'inf':
        mayor = 0.0
        for i in range(len(v)):
            if abs(v[i]) > mayor:
                mayor = abs(v[i])
        return mayor

    raise ValueError("Tipo de norma no soportado")


# Norma de Frobenius de una matriz
def norma_frobenius(A):
    n = len(A)
    m = len(A[0])
    suma = 0.0

    for i in range(n):
        for j in range(m):
            suma += A[i][j]*A[i][j]

    return math.sqrt(suma)


# Sustitucion backward
def sustitucion_backward(U, b):
    n = len(b)
    x = np.zeros(n, dtype=float)

    for k in range(n - 1, -1, -1):
        if abs(U[k][k]) < 1e-15:
            raise ValueError("Elemento diagonal nulo")

        suma = 0.0

        for j in range(k + 1, n):
            suma += U[k][j]*x[j]

        x[k] = (b[k] - suma)/U[k][k]

    return x


# Eliminacion de Gauss con pivoteo parcial
def gauss_pivoteo(A, b):
    n = len(b)

    M = np.array(A, dtype=float, copy=True)
    v = np.array(b, dtype=float, copy=True)

    for k in range(n - 1):

        # Buscando el pivote
        p = k
        mayor = abs(M[k][k])

        for i in range(k + 1, n):
            if abs(M[i][k]) > mayor:
                mayor = abs(M[i][k])
                p = i

        if mayor < 1e-15:
            raise ValueError("Matriz singular o casi singular")

        # Intercambiando filas
        if p != k:
            for j in range(n):
                temporal = M[k][j]
                M[k][j] = M[p][j]
                M[p][j] = temporal

            temporal = v[k]
            v[k] = v[p]
            v[p] = temporal

        # Eliminacion hacia adelante
        for i in range(k + 1, n):
            factor = M[i][k]/M[k][k]

            M[i][k] = 0.0

            for j in range(k + 1, n):
                M[i][j] -= factor*M[k][j]

            v[i] -= factor*v[k]

    return sustitucion_backward(M, v)


# Inversa de una matriz por eliminacion de Gauss
def inversa_gauss(A):
    n = len(A)
    A_inv = np.zeros((n, n), dtype=float)

    for j in range(n):
        unidad = np.zeros(n, dtype=float)
        unidad[j] = 1.0

        columna = gauss_pivoteo(A, unidad)

        for i in range(n):
            A_inv[i][j] = columna[i]

    return A_inv


# Numero de condicion de Frobenius
def condicion_frobenius(A):
    inversa = inversa_gauss(A)

    return norma_frobenius(A)*norma_frobenius(inversa)


# Calculando el residuo r = Ax - b
def residuo(A, x, b):
    n = len(b)
    r = np.zeros(n, dtype=float)

    for i in range(n):
        suma = 0.0

        for j in range(n):
            suma += A[i][j]*x[j]

        r[i] = suma - b[i]

    return r

##########################################
# METODOS ITERATIVOS

# Metodo Jacobi
def jacobi(A, b, tol=1e-6, max_iter=50):
    n = len(b)
    x = np.zeros(n, dtype=float)
    errores = []

    # Verificando la diagonal
    for i in range(n):
        if abs(A[i][i]) < 1e-15:
            raise ValueError("Elemento diagonal nulo")

    # Comenzando las iteraciones
    for k in range(1, max_iter + 1):

        x_nuevo = np.zeros(n, dtype=float)

        for i in range(n):
            suma = 0.0

            for j in range(n):
                if j != i:
                    suma += A[i][j]*x[j]

            x_nuevo[i] = (b[i] - suma)/A[i][i]

        # Calculamos la diferencia de iteraciones
        delta = np.zeros(n, dtype=float)

        for i in range(n):
            delta[i] = x_nuevo[i] - x[i]

        # Error relativo porcentual
        norma_nueva = norma_vec(x_nuevo, 2)

        if norma_nueva == 0:
            if norma_vec(delta, 2) == 0:
                error = 0.0
            else:
                error = float('inf')
        else:
            error = 100.0*norma_vec(delta, 2)/norma_nueva

        errores.append(error)

        # Actualizamos
        x = x_nuevo

        if error < tol:
            return x, k, errores

    return x, max_iter, errores


# Metodo SOR
def sor(A, b, omega=1.0, tol=1e-6, max_iter=50):
    n = len(b)
    x = np.zeros(n, dtype=float)
    errores = []

    if omega <= 0 or omega >= 2:
        raise ValueError("Omega debe estar entre 0 y 2")

    for i in range(n):
        if abs(A[i][i]) < 1e-15:
            raise ValueError("Elemento diagonal nulo")

    # Comenzando las iteraciones
    for k in range(1, max_iter + 1):

        x_anterior = x.copy()

        for i in range(n):
            suma = 0.0

            for j in range(n):
                if j != i:
                    suma += A[i][j]*x[j]

            # Calculamos Gauss-Seidel
            x_gs = (b[i] - suma)/A[i][i]

            # Aplicamos relajacion
            x[i] = (1.0 - omega)*x[i] + omega*x_gs

        # Calculando delta
        delta = np.zeros(n, dtype=float)

        for i in range(n):
            delta[i] = x[i] - x_anterior[i]

        # Error relativo porcentual
        norma_x = norma_vec(x, 2)

        if norma_x == 0:
            if norma_vec(delta, 2) == 0:
                error = 0.0
            else:
                error = float('inf')
        else:
            error = 100.0*norma_vec(delta, 2)/norma_x

        errores.append(error)

        if error < tol:
            return x, k, errores

    return x, max_iter, errores


# Gauss-Seidel como caso particular de SOR
def gauss_seidel(A, b, tol=1e-6, max_iter=50):
    return sor(A, b, omega=1.0, tol=tol, max_iter=max_iter)

##########################################
# PROBLEMA 1: ESFERAS CONDUCTORAS

# Construccion de la matriz de potenciales
def matriz_potenciales(a, d=1.0):
    c = 8.99e9
    n = 6

    P = np.zeros((n, n), dtype=float)

    for i in range(n):
        for j in range(n):

            if i == j:
                P[i][j] = c/a

            else:
                # Distancia entre los vertices
                pasos = abs(i - j)

                if pasos > 3:
                    pasos = 6 - pasos

                if pasos == 1:
                    distancia = d

                elif pasos == 2:
                    distancia = math.sqrt(3.0)*d

                else:
                    distancia = 2.0*d

                P[i][j] = c/distancia

    return P


# Grafica de convergencia
def graficar_convergencia(err_j, err_gs, a):

    iter_j = list(range(1, len(err_j) + 1))
    iter_gs = list(range(1, len(err_gs) + 1))

    plt.figure(figsize=(9, 5))

    plt.semilogy(iter_j, err_j, 'o-', label="Jacobi")
    plt.semilogy(iter_gs, err_gs, 's-', label="Gauss-Seidel")

    plt.xlabel("Numero de iteracion")
    plt.ylabel("Error relativo porcentual (%)")
    plt.title(f"Convergencia con a = {a} m")

    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    nombre = f"PC2_problema1_a_{str(a).replace('.', '_')}.png"
    plt.savefig(nombre, dpi=150)
    plt.show()


# Resolviendo los incisos b y c
def ejecutar_problema_1():

    print("\nPROBLEMA 1: ESFERAS CONDUCTORAS")

    V = np.array([100.0, 0.0, 50.0, 25.0, 75.0, 0.0])

    resultados_gs = []
    tol = 1e-6

    # Probando los dos radios
    for a in [0.1, 0.3]:

        P = matriz_potenciales(a)

        q_j, it_j, errores_j = jacobi(
            P, V, tol=tol, max_iter=50
        )

        q_g, it_g, errores_g = gauss_seidel(
            P, V, tol=tol, max_iter=50
        )

        # Verificando convergencia
        converge_j = errores_j[-1] < tol
        converge_g = errores_g[-1] < tol

        print(f"\nRadio a = {a} m")

        print("Jacobi:")
        print("Iteraciones:", it_j)

        if converge_j:
            print("Cargas (nC):", q_j*1e9)
        else:
            print("Jacobi no convergio")
            print("No utilizar las cargas como solucion")

        print("\nGauss-Seidel:")
        print("Iteraciones:", it_g)
        print("Convergio:", converge_g)
        print("Cargas (nC):", q_g*1e9)

        error_residuo = norma_vec(residuo(P, q_g, V), 'inf')
        print("Residuo (V):", error_residuo)

        resultados_gs.append(q_g.copy())

        # Generando la grafica
        graficar_convergencia(errores_j, errores_g, a)

    # Comparando los signos de las cargas
    print("\nComparacion de signos:")

    for i in range(6):

        q_01 = resultados_gs[0][i]
        q_03 = resultados_gs[1][i]

        if q_01*q_03 < 0:
            print(f"Esfera {i+1} cambia de signo")

            print("Carga con a=0.1:", q_01*1e9, "nC")
            print("Carga con a=0.3:", q_03*1e9, "nC")

    print("\nConclusion:")
    print("La esfera 4 cambia de signo.")
    print("Al aumentar el radio, la influencia relativa")
    print("de las demas esferas es mayor.")
    
    ##########################################
# PROBLEMA 2: SISTEMA DE RESORTES

# Construyendo la matriz K
def matriz_resortes(k, kc):

    K = np.zeros((2, 2), dtype=float)

    K[0][0] = k + kc
    K[0][1] = -kc

    K[1][0] = -kc
    K[1][1] = k + kc

    return K


def ejecutar_problema_2():

    print("\nPROBLEMA 2: SISTEMA DE RESORTES")

    # Inciso a
    print("Ecuaciones de equilibrio:")
    print("(k+kc)x1 - kc*x2 = F1")
    print("-kc*x1 + (k+kc)x2 = F2")

    ######################################
    # Inciso b: numero de condicion

    k = 10.0
    lista_kc = [1.0, 10.0, 100.0, 1000.0, 10000.0]

    eje_teorico = []
    eje_numerico = []

    print("\nTabla de numeros de condicion")
    print(f"{'kc/k':>10} | {'kappa_F':>15} | {'r+1/r':>15}")
    print("-"*48)

    for kc in lista_kc:

        K = matriz_resortes(k, kc)

        relacion = kc/k

        r = 1.0 + 2.0*relacion

        teorico = r + 1.0/r

        numerico = condicion_frobenius(K)

        print(
            f"{relacion:10.1f} | "
            f"{numerico:15.8f} | "
            f"{teorico:15.8f}"
        )

        eje_teorico.append(teorico)
        eje_numerico.append(numerico)

    # Grafica log-log
    plt.figure(figsize=(8, 5))

    plt.loglog(
        eje_teorico, eje_numerico,
        'o-', label="Kappa F numerico"
    )

    plt.loglog(
        eje_teorico, eje_teorico,
        '--', label="Referencia y = x"
    )

    plt.xlabel("r + 1/r")
    plt.ylabel("Numero de condicion de Frobenius")
    plt.title("Condicionamiento del sistema")

    plt.grid(True, which="both", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    plt.savefig("PC2_problema2_condicion.png", dpi=150)
    plt.show()

    print("\nConclusion:")
    print("El condicionamiento empeora cuando kc/k aumenta.")

    ######################################
    # Inciso c: perturbaciones aleatorias

    kc = 1000.0

    K = matriz_resortes(k, kc)

    F = np.array([1.0, -1.0])

    # Resolviendo el sistema original
    x = gauss_pivoteo(K, F)

    N = 1000
    factores = np.zeros(N, dtype=float)

    # Para obtener resultados reproducibles
    generador = random.Random(2026)

    norma_F = norma_vec(F, 2)
    norma_x = norma_vec(x, 2)

    for i in range(N):

        # Angulo aleatorio
        alfa = generador.uniform(0.0, 2.0*math.pi)

        # Construimos delta_F
        delta_F = np.zeros(2, dtype=float)

        delta_F[0] = 0.01*norma_F*math.cos(alfa)
        delta_F[1] = 0.01*norma_F*math.sin(alfa)

        # Nuevas fuerzas
        F_perturbada = F + delta_F

        # Resolviendo el sistema perturbado
        x_perturbada = gauss_pivoteo(K, F_perturbada)

        # Variacion de la solucion
        delta_x = x_perturbada - x

        # Errores relativos
        error_rel_x = norma_vec(delta_x, 2)/norma_x
        error_rel_F = norma_vec(delta_F, 2)/norma_F

        # Factor de amplificacion
        factores[i] = error_rel_x/error_rel_F

    # Buscando el mayor factor de amplificacion
    maximo = factores[0]

    for i in range(1, N):
        if factores[i] > maximo:
            maximo = factores[i]

    kappa_F = condicion_frobenius(K)

    print("\nResultados del inciso c:")
    print("Solucion x:", x)
    print("Factor maximo:", maximo)
    print("Numero de condicion:", kappa_F)

    # Histograma
    plt.figure(figsize=(9, 5))

    plt.hist(
        factores,
        bins=np.arange(0, 211, 30),
        edgecolor="black",
        label="1000 perturbaciones"
    )

    plt.axvline(
        kappa_F,
        color="red",
        linestyle="--",
        label="Kappa F"
    )

    plt.xlabel("Factor de amplificacion")
    plt.ylabel("Frecuencia")
    plt.title("Amplificacion de errores")

    plt.xlim(0, 210)
    plt.xticks(np.arange(0, 211, 30))

    plt.grid(True, axis="y", linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    plt.savefig("PC2_problema2_histograma.png", dpi=150)
    plt.show()
    
    ##########################################
# PROBLEMA 3: CIRCUITO DE RESISTENCIAS

# Construyendo A*I = b para n celdas
def construir_sistema_circuito(n, R=2.0, V=1.5):

    if n < 2:
        raise ValueError("n debe ser al menos 2")

    # Numero de corrientes
    m = 2*n - 1

    # Matriz y vector del sistema
    A = np.zeros((m, m), dtype=float)
    b = np.zeros(m, dtype=float)

    fila = 0

    ######################################
    # Primera malla:
    # 2R*I1 + R*I2 = V

    A[fila][0] = 2.0*R
    A[fila][1] = R

    b[fila] = V
    fila += 1

    ######################################
    # Ecuaciones de nodo:
    # I_(2k-1) - I_(2k) - I_(2k+1) = 0

    for k in range(1, n):

        A[fila][2*k - 2] = 1.0
        A[fila][2*k - 1] = -1.0
        A[fila][2*k] = -1.0

        fila += 1

    ######################################
    # Mallas interiores:
    # R*I_(2k-2) - 2R*I_(2k-1) - R*I_(2k) = 0

    for k in range(2, n):

        A[fila][2*k - 3] = R
        A[fila][2*k - 2] = -2.0*R
        A[fila][2*k - 1] = -R

        fila += 1

    ######################################
    # Ultima malla:
    # R*I_(2n-2) - 3R*I_(2n-1) = 0

    A[fila][2*n - 3] = R
    A[fila][2*n - 2] = -3.0*R

    fila += 1

    # Verificando el numero de ecuaciones
    if fila != m:
        raise ValueError("Numero de ecuaciones incorrecto")

    return A, b


# Imprimiendo matrices ordenadamente
def imprimir_matriz_y_vector(A, b):

    n = len(b)

    for i in range(n):
        print("[", end="")

        for j in range(n):
            print(f"{A[i][j]:6.1f}", end=" ")

        print(f"]  |  {b[i]:6.2f}")


# Resolviendo los incisos
def ejecutar_problema_3():

    print("\nPROBLEMA 3: CIRCUITO DE RESISTENCIAS")

    ######################################
    # Inciso b: matriz para n=6

    A6, b6 = construir_sistema_circuito(n=6)

    print("\nMatriz para n=6:")
    print("Dimension:", len(A6), "x", len(A6))

    imprimir_matriz_y_vector(A6, b6)

    ######################################
    # Inciso c: resolver para n=12

    A12, b12 = construir_sistema_circuito(n=12)

    # Resolviendo mediante eliminacion de Gauss
    corrientes = gauss_pivoteo(A12, b12)

    print("\nCorrientes para n=12:")

    for i in range(len(corrientes)):
        print(f"I_{i+1:2d} = {corrientes[i]:.12f} A")

    # Comprobando el resultado
    error = norma_vec(
        residuo(A12, corrientes, b12),
        'inf'
    )

    print("\nNorma infinito del residuo:", error)

    ######################################
    # Grafica de las corrientes

    valores_k = list(range(1, len(corrientes) + 1))

    plt.figure(figsize=(10, 5))

    plt.plot(
        valores_k,
        corrientes,
        'o-',
        label="Corrientes I_k"
    )

    plt.xlabel("Indice k")
    plt.ylabel("Corriente I_k (A)")
    plt.title("Corrientes del circuito para n=12")

    plt.xticks(valores_k)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.legend()
    plt.tight_layout()

    plt.savefig("PC2_problema3_corrientes.png", dpi=150)
    plt.show()
    
##########################################
# EJECUTANDO LOS TRES PROBLEMAS

ejecutar_problema_1()

ejecutar_problema_2()

ejecutar_problema_3()