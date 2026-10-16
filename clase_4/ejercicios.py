import math
from modelos import sustitucion_backward, sustitucion_forward, gauss_pivoteo, calcular_inversa, norma_mat, numero_condicion,norma_vec
#ejercicio a)
U_a = [
    [1.0,  2.0,  1.0],
    [0.0, -4.0,  1.0],
    [0.0,  0.0, -2.0]
]
b_a = [5.0, 2.0, 4.0]

x_a = sustitucion_backward(U_a, b_a)
print("Solución (a) Backward: ", x_a)

#ejercicio b)
L_b = [
    [2.0, 0.0, 0.0],
    [1.0, 4.0, 0.0],
    [4.0, 3.0, 3.0]
]
b_b = [4.0, 2.0, 5.0]

x_b = sustitucion_forward(L_b, b_b)
print("Solución (b) Forward:  ", x_b)

#Ejercicio de las resistencias
def resolver_red_resistores(N, V_plus=5.0):
    #Creamos la matriz A (NxN) y el vector W (N) llenos de ceros
    A = [[0.0] * N for _ in range(N)]
    W = [0.0] * N
    #Algoritmo para completar la matriz en base a lo encontrado manualmente
    for i in range(N):
        #La diagonal principal es 4, excepto en los extremos que es 3
        if i == 0 or i == N - 1:
            A[i][i] = 3.0
        else:
            A[i][i] = 4.0
        #Reglas para las condiciones con sus potenciales vecinos (asignamos -1.0 si existen dentro del rango)
        if i + 1 < N: A[i][i+1] = -1.0
        if i - 1 >= 0: A[i][i-1] = -1.0
        if i + 2 < N: A[i][i+2] = -1.0
        if i - 2 >= 0: A[i][i-2] = -1.0
    #Condiciones de frontera (El voltaje entra por V1 y V2)
    W[0] = V_plus
    W[1] = V_plus
    # Resolvemos usando gauss_pivoteo
    return gauss_pivoteo(A, W)

#Solución para N = 6 
V_6 = resolver_red_resistores(6, 5.0)
print("Potenciales para N=6:")
for i, v in enumerate(V_6):
    print(f"V_{i+1} = {v:.4f} V")

# Solución para N = 1000 
# Solo calculamos e imprimimos los primeros y últimos para no inundar la pantalla
V_1000 = resolver_red_resistores(1000, 5.0)
print("\nPotenciales para N=1000 (Muestra):")
print(f"V_1 = {V_1000[0]:.4f} V")
print(f"V_2 = {V_1000[1]:.4f} V")
print(f"... decae a lo largo de la red ...")
print(f"V_999 = {V_1000[998]:.6f} V")
print(f"V_1000 = {V_1000[999]:.6f} V")

#Solucion a el sistema de basas resorte

def resolver_sistema_masa_resorte(k, m1, m2, m3, g=9.81):
    #Construimos la matriz de rigidez K
    K = [
        [3*k, -2*k,  0.0],
        [-2*k, 3*k, -1.0*k],
        [0.0, -1.0*k, 1.0*k]
    ]
    #Construimos el vector de fuerzas (pesos)
    F = [
        m1 * g,
        m2 * g,
        m3 * g
    ]
    #Resolvemos el sistema Kx = F usando nuestra función de Gauss
    desplazamientos = gauss_pivoteo(K, F)
    
    return desplazamientos

k_val = 10.0      # N/m
m_val = 1.0       # kg
gravedad = 9.81   # m/s^2
x_solucion = resolver_sistema_masa_resorte(k_val, m_val, m_val, m_val, gravedad)
print("Soluciones Estacionarias del Sistema Masa-Resorte:")
print(f"x_1 = {x_solucion[0]:.4f} m")
print(f"x_2 = {x_solucion[1]:.4f} m")
print(f"x_3 = {x_solucion[2]:.4f} m")

#Solucion del problema densidad de una varilla

def crear_matriz_hilbert(n):
    #Creando la matriz de Hilbert H_n según h_ij = 1 / (i + j - 1)
    H = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            # Usamos +1 porque nuestros índices de Python empiezan en 0
            H[i][j] = 1.0 / ((i + 1) + (j + 1) - 1.0)
    return H


print(" PARTE (B): DENSIDAD LINEAL (N=2)")


H2 = crear_matriz_hilbert(2)
M_real = [2.5, 1.5]

#Encontrando los coeficientes c
c_real = gauss_pivoteo(H2, M_real)
print(f"Coeficientes reales c: {c_real}")

#Inversa y Condicionamiento 
H2_inv = calcular_inversa(H2)
norma_H2 = norma_mat(H2, tipo='F')
norma_H2_inv = norma_mat(H2_inv, tipo='F')
kappa_H2 = numero_condicion(H2, tipo='F')

print(f"Norma ||H2||_F: {norma_H2:.4f}")
print(f"Norma ||H2^-1||_F: {norma_H2_inv:.4f}")
print(f"Número de condición kappa_F: {kappa_H2:.4f}")

#Perturbación (Error de 1 mm en x_cm)
M_err = [2.5, 1.5025]
c_err = gauss_pivoteo(H2, M_err)
print(f"Coeficientes con error c_err: {c_err}")

#Verificación de la cota
delta_c = [c_err[0] - c_real[0], c_err[1] - c_real[1]]
delta_M = [M_err[0] - M_real[0], M_err[1] - M_real[1]]

error_relativo_c = norma_vec(delta_c, tipo=2) / norma_vec(c_real, tipo=2)
error_relativo_M = norma_vec(delta_M, tipo=2) / norma_vec(M_real, tipo=2)

print(f"\nError relativo en M: {error_relativo_M:.5f}")
print(f"Error relativo en c: {error_relativo_c:.5f}")
print(f"¿Se cumple la cota ({error_relativo_c:.5f} <= {kappa_H2 * error_relativo_M:.5f})?: {error_relativo_c <= (kappa_H2 * error_relativo_M)}")


print("\n")
print(" PARTE (C): DENSIDAD CUADRÁTICA (N=3)")


H3 = crear_matriz_hilbert(3)
M3_real = [5.5, 3.25, 2.35]

c3_real = gauss_pivoteo(H3, M3_real)
print(f"Coeficientes reales c (n=3): {c3_real}")

# Condicionamiento de H3
kappa_H3 = numero_condicion(H3, tipo='F')
print(f"Número de condición kappa_F (H3): {kappa_H3:.4f}")
print("-> ¡Nota cómo creció masivamente frente al 19 de H2!")

# Perturbación del 1% en M2
M3_err = [5.5, 3.25, 2.3735]
c3_err = gauss_pivoteo(H3, M3_err)
print(f"\nCoeficientes perturbados c_err (n=3): {c3_err}")

delta_c3 = [c3_err[0] - c3_real[0], c3_err[1] - c3_real[1], c3_err[2] - c3_real[2]]
err_rel_c3 = norma_vec(delta_c3, tipo=2) / norma_vec(c3_real, tipo=2)

print(f"Error relativo en la densidad final (n=3): {err_rel_c3 * 100:.2f}%")
print("-> Conclusión: La inestabilidad matemática arruinó la reconstrucción.")

#Solucion del sistema de masas
print(" SISTEMA MASA-RESORTE")

# Datos del problema
k1 = 1.0      # N/m (blando)
k2 = 1000.0   # N/m (rígido)

# a) Matriz de Rigidez K
K = [
    [k1 + k2, -k2],
    [-k2,     k2]
]
print("Matriz de Rigidez K:")
print(f"[{K[0][0]}, {K[0][1]}]")
print(f"[{K[1][0]}, {K[1][1]}]")

# b) Condicionamiento
K_inv = calcular_inversa(K)
kappa_K = numero_condicion(K, tipo='F')
cifras_perdidas = math.log10(kappa_K)

print(f"\nNúmero de condición kappa_F: {kappa_K:.2f}")
print(f"Cifras significativas perdidas: ~{cifras_perdidas:.2f}")

# c) Resolución con fuerzas F = (-1, 1)
F = [-1.0, 1.0]
x_real = gauss_pivoteo(K, F)

print("\nResultados con F = (-1, 1) N:")
print(f"x1 = {x_real[0]:.4f} m")
print(f"x2 = {x_real[1]:.4f} m")

# d) Resolución con perturbación en F1
F_err = [-0.999, 1.0]
x_err = gauss_pivoteo(K, F_err)

print("\nResultados con F1 perturbada a -0.999 N:")
print(f"x1 = {x_err[0]:.4f} m")
print(f"x2 = {x_err[1]:.4f} m")

# Cálculo de errores relativos
delta_x = [x_err[0] - x_real[0], x_err[1] - x_real[1]]
delta_F = [F_err[0] - F[0], F_err[1] - F[1]]

err_rel_F = norma_vec(delta_F, tipo=2) / norma_vec(F, tipo=2)
err_rel_x = norma_vec(delta_x, tipo=2) / norma_vec(x_real, tipo=2)

print(f"\nError relativo en la Fuerza: {err_rel_F * 100:.4f}%")
print(f"Error relativo en la Posición: {err_rel_x * 100:.2f}%")

cota_teorica = kappa_K * err_rel_F
print(f"Verificación de la Cota: {err_rel_x:.4f} <= {cota_teorica:.4f}")