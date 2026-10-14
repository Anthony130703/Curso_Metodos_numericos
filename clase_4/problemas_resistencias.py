def armar_sistema_resistores(N, V_plus=5.0):
    A = np.zeros((N, N), dtype=float)
    W = np.zeros(N, dtype=float)
    
    for i in range(N):
        # Diagonal principal: 3 en los extremos, 4 en el interior
        A[i, i] = 3.0 if (i == 0 or i == N - 1) else 4.0
        # Primeros vecinos
        if i - 1 >= 0: A[i, i - 1] = -1.0
        if i + 1 < N:  A[i, i + 1] = -1.0
        # Segundos vecinos
        if i - 2 >= 0: A[i, i - 2] = -1.0
        if i + 2 < N:  A[i, i + 2] = -1.0
        
    W[0] = V_plus
    W[1] = V_plus
    return A, W

# Caso N = 6 (Método directo: Cholesky o Gauss)
A6, W6 = armar_sistema_resistores(6, V_plus=5.0)
V6 = resolver_cholesky(A6, W6)
print("Potenciales para N = 6 (V):", np.round(V6, 4))
# Salida: [3.7209, 3.4884, 2.6744, 2.3256, 1.5116, 1.2791]

# Caso N = 1000 (Método iterativo: Gradiente Conjugado)
A1000, W1000 = armar_sistema_resistores(1000, V_plus=5.0)
V1000, it_cg, _ = gradiente_conjugado(A1000, W1000, tol=1e-8, max_iter=2000)
print(f"N = 1000 resuelto en {it_cg} iteraciones con CG.")
print("Primeros 4 nodos (V):", np.round(V1000[:4], 4))
print("Últimos 4 nodos (V): ", np.round(V1000[-4:], 4))