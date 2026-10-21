from modelos import jacobi, gauss_seidel, sor, maximo_descenso, gradiente_conjugado
A = [
    [4, 2, 3],
    [3, -5, 2],
    [-2, 3, 8]
]

b = [8, -14, 27]

x1, k1 = jacobi(A, b)
x2, k2 = gauss_seidel(A, b)
x3, k3 = sor(A, b, omega=1.2)

print("Jacobi:", x1, "Iteraciones:", k1)
print("Gauss-Seidel:", x2, "Iteraciones:", k2)
print("SOR:", x3, "Iteraciones:", k3)

A = [
    [4, 1],
    [1, 3]
]

b = [1, 2]

x1, k1 = maximo_descenso(A, b)
x2, k2 = gradiente_conjugado(A, b)

print("Maximo descenso:", x1, "Iteraciones:", k1)
print("Gradiente conjugado:", x2, "Iteraciones:", k2)