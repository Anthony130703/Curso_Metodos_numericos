import numpy as np
e = np.zeros((3, 3, 3))
e[0,1,2] = e[1,2,0] = e[2,0,1] = 1.0
e[0,2,1] = e[2,1,0] = e[1,0,2] = -1.0
d = np.eye(3)

# identidad e-delta: e_{abg} e_{mng} = d_am d_bn - d_an d_bm
lhs = np.einsum('abg, mng -> abmn', e, e)
rhs = np.einsum('am, bn -> abmn', d, d) - np.einsum('an, bm -> abmn', d, d)
print(np.allclose(lhs, rhs))

#BAC-CAB con vectores al azar
A, B, C = np.random.randn(3, 3)
print(np.allclose(np.cross(A, np.cross(B, C)), B*(A@C) - C*(A@B))) # true

# Promedio angular <n_a n_b> = delta_ab / 3
u = np.random.randn(200000, 3); n = u/np.linalg.norm(u, axis = 1, keepdims = True)
print(np.round(np.einsum('ka, kb -> ab', n, n)/len(n), 3))  # ~ I/3