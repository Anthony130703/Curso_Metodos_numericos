import numpy as np
import matplotlib.pyplot as plt
#Solucion item b) de la pregunta 1
c = 8.99 * 1e9
A_iter = [[ 10*c,  c, (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c, c],
          [ c, 10*c,  c, (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c],
          [ (1/np.sqrt(3))*c, c, 10*c,  c, (1/np.sqrt(3))*c, c/2],
          [ c/2, (1/np.sqrt(3))*c, c, 10*c,  c, (1/np.sqrt(3))*c],
          [ (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c, c, 10*c,  c],
          [  c, (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c, c, 10*c]]
b_iter = [100, 0, 50, 25, 75, 0]

def norma_vec(v, tipo=2):
    #Calculando la norma vectorial para los errores
    v = np.asarray(v, dtype=float)
    if tipo == 1:
        return float(np.sum(np.abs(v)))
    elif tipo == 2:
        return float(np.sqrt(np.sum(v ** 2)))
    elif tipo == 'inf':
        return float(np.max(np.abs(v)))
    else:
        raise ValueError("Tipo de norma no soportado. Use 1, 2 o 'inf'.")
    
def jacobi(A, b, x0=None, tol=1e-6, max_iter=1000):
    
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
        norma_x_new = norma_vec(x_new, tipo=2)
        error_relativo = (
            delta / norma_x_new * 100
            if norma_x_new != 0
            else (0.0 if delta == 0 else np.inf)
        )
        historial_err.append(error_relativo)
        x = x_new

        if error_relativo < tol:
            return x, k, historial_err

    return x, max_iter, historial_err

def sor(A, b, omega=1.2, x0=None, tol=1e-6, max_iter=1000):
    
    # Si omega = 1.0, es exactamente Gauss-Seidel
    # Retorna: (x, iteraciones, historial_errores)
    
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
        norma_x = norma_vec(x, tipo= 2)
        error_relativo = (
            delta / norma_x * 100
            if norma_x != 0
            else (0.0 if delta == 0 else np.inf)
        )
        historial_err.append(error_relativo)
        
        if error_relativo < tol:
            return x, k, historial_err


    return x, max_iter, historial_err

def gauss_seidel(A, b, x0=None, tol=1e-6, max_iter=1000):
    # caso particular de SOR con omega = 1.0)
    # Retorna: (x, iteraciones, historial_errores)
    
    return sor(A, b, omega=1.0, x0=x0, tol=tol, max_iter=max_iter)

q_jac, it_jac, erro_jac = jacobi(A_iter, b_iter, tol=1e-6, max_iter=50)
q_gs,  it_gs,  erro_gs = gauss_seidel(A_iter, b_iter, tol=1e-6, max_iter=50)
print(f"Jacobi ({it_jac} iter):      {q_jac}")
#print(f"Error relativo porcentual por iteración: {erro_jac}")
print(f"Gauss-Seidel ({it_gs} iter):      {q_gs}")
#print(f"Error relativo porcentual por iteración: {erro_gs}")

iteraciones = np.arange(1, len(erro_jac) + 1)
iteraciones_gs = np.arange(1, len(erro_gs) + 1)
plt.plot(iteraciones, erro_jac, marker="o", label="Jacobi")
plt.plot(iteraciones_gs, erro_gs, marker="s", label="Gauss-Seidel")
plt.xlabel("Número de iteración")
plt.ylabel("Error relativo aproximado (%)")
plt.title("Convergencia de Jacobi y Gauss-Seidel")
plt.yscale("log")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

#Solucion item c) de la pregunta 1

A_iter_c = [[ (10/3)*c,  c, (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c, c],
          [ c, (10/3)*c,  c, (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c],
          [ (1/np.sqrt(3))*c, c, (10/3)*c,  c, (1/np.sqrt(3))*c, c/2],
          [ c/2, (1/np.sqrt(3))*c, c, (10/3)*c,  c, (1/np.sqrt(3))*c],
          [ (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c, c, (10/3)*c,  c],
          [  c, (1/np.sqrt(3))*c, c/2, (1/np.sqrt(3))*c, c, (10/3)*c]]

q_jac_c, it_jac_c, erro_jac_c = jacobi(A_iter_c, b_iter, tol=1e-6, max_iter=50)
q_gs_c,  it_gs_c,  erro_gs_c = gauss_seidel(A_iter_c, b_iter, tol=1e-6, max_iter=50)
print(f"Jacobi ({it_jac_c} iter):      {q_jac_c}")
#print(f"Error relativo porcentual por iteración: {erro_jac}")
print(f"Gauss-Seidel ({it_gs_c} iter):      {q_gs_c}")
#print(f"Error relativo porcentual por iteración: {erro_gs}")

for i, (elem1, elem2) in enumerate(zip(q_gs, q_gs_c)):
    # Multiplicamos los elementos correspondientes
    producto = elem1 * elem2
    
    # Verificamos si NO son mayores que 0 (es decir, menores o iguales a 0)
    if producto <= 0:
        print(f"La carga numero {i +1} ({elem1}, {elem2}) cambio de signo")
        

iteraciones = np.arange(1, len(erro_jac) + 1)
iteraciones_gs = np.arange(1, len(erro_gs) + 1)
plt.plot(iteraciones, erro_jac, marker="o", label="Jacobi")
plt.plot(iteraciones_gs, erro_gs, marker="s", label="Gauss-Seidel")
plt.xlabel("Número de iteración")
plt.ylabel("Error relativo aproximado (%)")
plt.title("Convergencia de Jacobi y Gauss-Seidel")
plt.yscale("log")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.show()

#Problema 2

#Solucion del item b)
k_c = [1, 10, 100, 1000, 10000]
k = 10
matrices = []

for kc in k_c:
    k_matriz = np.array([[k + kc, -kc],
                         [-kc, k + kc]], dtype=float)
    matrices.append(k_matriz)
        
def sustitucion_backward(U, b):
    n = len(b)
    x = np.zeros(n, dtype=float)
    for k in range(n - 1, -1, -1):
        if abs(U[k, k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = np.dot(U[k, k+1:], x[k+1:])
        x[k] = (b[k] - suma) / U[k, k]
    return x
        
def gauss_pivoteo(A, b, usar_pivoteo=True):
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

def inversa_gauss(A):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    A_inv = np.zeros((n, n), dtype=float)
    I = np.eye(n, dtype=float)

    for j in range(n):
        A_inv[:, j] = gauss_pivoteo(A, I[:, j], usar_pivoteo=True)

    return A_inv

def norma_mat(A, tipo='F'):
    A = np.asarray(A, dtype=float)
    if tipo == 'F':
        return float(np.sqrt(np.sum(A ** 2)))
    elif tipo == 'inf':
        return float(np.max(np.sum(np.abs(A), axis=1)))
    elif tipo == 1:
        return float(np.max(np.sum(np.abs(A), axis=0)))
    else:
        raise ValueError("Use 'F', 'inf' o 1.")
    
def numero_condicion(A, tipo='F'):
    #Calculando kappa(A) = ||A|| * ||A^{-1}|| en la norma indicada
    A_inv = inversa_gauss(A)
    return norma_mat(A, tipo=tipo) * norma_mat(A_inv, tipo=tipo)

resultados = []
for kc, k_matriz in zip(k_c, matrices):
    kc_sobre_k = kc / k
    r = 1 + 2 * kc_sobre_k
    kappa_f_numerico = numero_condicion(k_matriz, tipo='F')
    r_mas_inverso = r + 1 / r
    resultados.append((kc_sobre_k, kappa_f_numerico, r_mas_inverso))

print("\nTabla de números de condición:")
print(f"{'kc/k':>10} | {'κF numérico':>16} | {'r + 1/r':>16}")
print("-" * 48)
for kc_sobre_k, kappa_f_numerico, r_mas_inverso in resultados:
    print(f"{kc_sobre_k:10.1f} | {kappa_f_numerico:16.8f} | {r_mas_inverso:16.8f}")

valores_r = [resultado[2] for resultado in resultados]
valores_kappa_f = [resultado[1] for resultado in resultados]
plt.loglog(valores_r, valores_kappa_f, marker="o")
plt.xlabel("r + 1/r")
plt.ylabel("κF numérico")
plt.title("Número de condición de Frobenius")
plt.grid(True, which="both", linestyle="--", alpha=0.5)
plt.tight_layout()
plt.show()