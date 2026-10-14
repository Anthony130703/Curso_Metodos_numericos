from modelos import jacobi, gauss_seidel, resolver_cholesky, gradiente_conjugado, gauss_pivoteo, numero_condicion
import numpy as np

# Ejemplo 1: Sistema 3x3 de las diapositivas de Jacobi / Gauss-Seidel (Solución: [-1, 3, 2])
A_iter = [[ 4,  2, 3],
          [ 3, -5, 2],
          [-2,  3, 8]]
b_iter = [8, -14, 27]

x_jac, it_jac, _ = jacobi(A_iter, b_iter, tol=1e-7)
x_gs,  it_gs,  _ = gauss_seidel(A_iter, b_iter, tol=1e-7)
print(f"Jacobi ({it_jac} iter):      { np.round(x_jac, 5) }")
print(f"Gauss-Seidel ({it_gs} iter): { np.round(x_gs, 5) }")

# Ejemplo 2: Masas y resortes acoplados (k=10, m1=m2=m3=1 kg, g=9.81 m/s^2)
# Matriz SPD -> Podemos usar Cholesky o Gradiente Conjugado
k, g = 10.0, 9.81
A_masas = [[ 3*k, -2*k,    0],
           [-2*k,  3*k,   -k],
           [   0,   -k,    k]]
b_masas = [1*g, 1*g, 1*g]

x_chol = resolver_cholesky(A_masas, b_masas)
x_cg, it_cg, _ = gradiente_conjugado(A_masas, b_masas, tol=1e-10)
print(f"Posiciones estacionarias (Cholesky): {np.round(x_chol, 4)} m")
print(f"Posiciones estacionarias (CG, {it_cg} it):  {np.round(x_cg, 4)} m")

# Ejemplo 3: Matriz de Hilbert H_n para la Densidad de la Varilla (Diapositivas 44-45)
def matriz_hilbert(n):
    return np.array([[1.0 / (i + j + 1) for j in range(n)] for i in range(n)])

H2 = matriz_hilbert(2)
M2 = np.array([2.5, 2.5 * 0.6])
c2 = gauss_pivoteo(H2, M2)
kappa_H2 = numero_condicion(H2, tipo='F')
print(f"Coeficientes de densidad (n=2): c = {np.round(c2, 4)}, kappa_F(H2) = {kappa_H2:.4f}")