import numpy as np
from modelos import gauss_pivoteo, inversa_gauss, numero_condicion, norma_vec
# --- PROBLEMA 2: DENSIDAD DE LA VARILLA ---
def hilbert(n):
    return np.array([[1.0 / (i + j + 1) for j in range(n)] for i in range(n)])

# (b) Caso n = 2
H2 = hilbert(2)
M2 = np.array([2.5, 2.5 * 0.6])
c2 = gauss_pivoteo(H2, M2)
H2_inv = inversa_gauss(H2)
kF_H2 = numero_condicion(H2, tipo='F')

M2_pert = np.array([2.5, 2.5 * 0.601])
c2_pert = gauss_pivoteo(H2, M2_pert)
err_rel_c2 = norma_vec(c2_pert - c2, 2) / norma_vec(c2, 2)
cota_c2 = kF_H2 * (norma_vec(M2_pert - M2, 2) / norma_vec(M2, 2))
print(f"[Varilla n=2] c = {c2}, kappa_F = {kF_H2:.4f}, Err_rel = {err_rel_c2:.5f} <= Cota = {cota_c2:.5f}")

# (c) Caso n = 3
H3 = hilbert(3)
M3 = np.array([5.5, 3.25, 2.35])
c3 = gauss_pivoteo(H3, M3)
kF_H3 = numero_condicion(H3, tipo='F')

M3_pert = np.array([5.5, 3.25, 2.35 * 1.01])
c3_pert = gauss_pivoteo(H3, M3_pert)
print(f"[Varilla n=3] c original = {np.round(c3, 3)}, c perturbado = {np.round(c3_pert, 3)}, kappa_F = {kF_H3:.2f}")

# --- PROBLEMA 3: SISTEMA MASA-RESORTE ---
K = np.array([[1001.0, -1000.0],
              [-1000.0, 1000.0]])
F = np.array([-1.0, 1.0])
dF = np.array([0.001, 0.0])

x = gauss_pivoteo(K, F)
x_pert = gauss_pivoteo(K, F + dF)
dx = x_pert - x

kF_K = numero_condicion(K, tipo='F')
err_F = norma_vec(dF, 2) / norma_vec(F, 2)
err_x = norma_vec(dx, 2) / norma_vec(x, 2)
print(f"[Resortes] x = {x}, dx = {np.round(dx, 5)}, Amplificación = {err_x/err_F:.1f}, kappa_F = {kF_K:.2f}")