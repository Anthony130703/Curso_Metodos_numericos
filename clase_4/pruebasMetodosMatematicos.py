from modelos import (
    jacobi, gauss_seidel, sor,
    maximo_descenso, gradiente_conjugado,
    cholesky, es_simetrica
)

# Probando Jacobi, Gauss-Seidel y SOR
A = [
    [4, 2, 3],
    [3, -5, 2],
    [-2, 3, 8]
]

b = [8, -14, 27]

x1, k1, delta1 = jacobi(A, b)
x2, k2, delta2 = gauss_seidel(A, b)
x3, k3, delta3 = sor(A, b, omega=1.2)

print("Jacobi:")
print("Solucion:", x1)
print("Iteraciones:", k1)
print("Delta:", delta1)

print("\nGauss-Seidel:")
print("Solucion:", x2)
print("Iteraciones:", k2)
print("Delta:", delta2)

print("\nSOR:")
print("Solucion:", x3)
print("Iteraciones:", k3)
print("Delta:", delta3)


# Probando maximo descenso y gradiente conjugado
A = [
    [4, 1],
    [1, 3]
]

b = [1, 2]

x4, k4, delta4 = maximo_descenso(A, b)
x5, k5, delta5 = gradiente_conjugado(A, b)

print("\nMaximo descenso:")
print("Solucion:", x4)
print("Iteraciones:", k4)
print("Delta:", delta4)

print("\nGradiente conjugado:")
print("Solucion:", x5)
print("Iteraciones:", k5)
print("Delta:", delta5)


# Comprobando si la matriz es simetrica definida positiva
if not es_simetrica(A):
    raise ValueError("La matriz no es simetrica")

G = cholesky(A)

print("\nLa matriz es simetrica definida positiva")