import numpy as np
import matplotlib.pyplot as plt
import math

##########################################
# PROBLEMA 1: EPSILON DE MAQUINA

# a) Pseudocodigo
# 1. Inicializamos epsilon = 1
# 2. Guardamos el ultimo epsilon valido
# 3. Dividimos epsilon entre 2
# 4. Si 1 + epsilon == 1, terminamos
# 5. Si no, actualizamos epsilon valido y repetimos

# b) Implementacion

epsilon = 1.0
ultimo_epsilon_valido = epsilon
iteraciones = 0

while True:
    epsilon = epsilon / 2.0
    iteraciones += 1

    if 1.0 + epsilon == 1.0:
        break

    ultimo_epsilon_valido = epsilon

print("PROBLEMA 1")
print("Epsilon maquina:", ultimo_epsilon_valido)
print("Divisiones realizadas:", iteraciones)
print("Divisiones con epsilon valido:", iteraciones - 1)


# c) Calculando y graficando R_k

max_k = 60
valores_k = np.zeros(max_k + 1)
valores_R = np.zeros(max_k + 1)

for k in range(max_k + 1):
    epsilon_k = 2.0**(-k)

    R_k = ((1.0 + epsilon_k) - 1.0)/epsilon_k

    valores_k[k] = k
    valores_R[k] = R_k

plt.figure(figsize=(10, 6))
plt.scatter(valores_k, valores_R, s=20, color="black")
plt.title("R_k vs k")
plt.xlabel("k")
plt.ylabel("R_k")
plt.grid(True)
plt.savefig("problema1_Rk.png")
plt.show()

##########################################
# PROBLEMA 2: DERIVACION NUMERICA

# Funcion de posicion del oscilador
def posicion_oscilador(t):
    return 0.1 * math.cos(2.0*math.pi*t)


# Derivada exacta
def velocidad_exacta(t):
    return -0.1*(2.0*math.pi)*math.sin(2.0*math.pi*t)


# Diferencia hacia adelante
def diferencia_adelante(funcion, t, h):
    return (funcion(t + h) - funcion(t))/h


# Diferencia centrada
def diferencia_centrada(funcion, t, h):
    return (funcion(t + h) - funcion(t - h))/(2.0*h)


# Valores iniciales
t0 = 0.3
v_exacto = velocidad_exacta(t0)

# Matriz de 16 filas y 5 columnas
datos = np.zeros((16, 5), dtype=float)

for i in range(16):
    k = i + 1
    h = 10.0**(-k)

    # Calculando derivadas
    Df = diferencia_adelante(posicion_oscilador, t0, h)
    Dc = diferencia_centrada(posicion_oscilador, t0, h)

    # Errores relativos con signo, segun el enunciado
    error_f = (Df - v_exacto)/v_exacto
    error_c = (Dc - v_exacto)/v_exacto

    # Guardando los resultados
    datos[i][0] = h
    datos[i][1] = Df
    datos[i][2] = Dc
    datos[i][3] = error_f
    datos[i][4] = error_c


# a) Guardando los datos en un archivo .txt
np.savetxt(
    "errores derivada.txt",
    datos,
    fmt="%.15e",
    header="h Df Dc epsilon_f epsilon_c"
)


# b) Graficando errores vs h en escala log-log

valores_h = np.zeros(16)
errores_adelante = np.zeros(16)
errores_centrada = np.zeros(16)

for i in range(16):
    valores_h[i] = datos[i][0]

    # Magnitud de los errores para la grafica
    errores_adelante[i] = abs(datos[i][3])
    errores_centrada[i] = abs(datos[i][4])

plt.figure(figsize=(10, 6))

plt.loglog(
    valores_h,
    errores_adelante,
    "o-",
    label="Diferencia hacia adelante"
)

plt.loglog(
    valores_h,
    errores_centrada,
    "s-",
    label="Diferencia centrada"
)

plt.title("Error relativo vs h")
plt.xlabel("h")
plt.ylabel("Magnitud del error relativo")
plt.legend()
plt.grid(True, which="both")
plt.savefig("problema2_errores.png")
plt.show()


# c) Encontrando el h de menor error

indice_f = 0
indice_c = 0

for i in range(1, 16):

    if errores_adelante[i] < errores_adelante[indice_f]:
        indice_f = i

    if errores_centrada[i] < errores_centrada[indice_c]:
        indice_c = i


# Valores optimos teoricos
h_opt_f = math.sqrt(ultimo_epsilon_valido)
h_opt_c = ultimo_epsilon_valido**(1.0/3.0)


print("\nPROBLEMA 2")

print("Velocidad exacta:", v_exacto)

print("Mejor h hacia adelante:", valores_h[indice_f])
print("Error minimo hacia adelante:", errores_adelante[indice_f])

print("Mejor h centrada:", valores_h[indice_c])
print("Error minimo centrada:", errores_centrada[indice_c])

print("h optimo teorico hacia adelante:", h_opt_f)
print("h optimo teorico centrada:", h_opt_c)

##########################################
# PROBLEMA 3: FRACTAL DE NEWTON

# Funcion que clasifica la raiz alcanzada
def clasificar_newton(z0, raices, tol=1e-6, max_iter=50):
    z = z0

    for i in range(max_iter):

        # Comprobando convergencia a las raices
        for k in range(3):
            if abs(z - raices[k]) < tol:
                return k + 1

        # Evitando division entre cero
        if abs(z) < 1e-14:
            return 0

        # Aplicando Newton-Raphson
        z = z - (z**3 - 1.0)/(3.0*z**2)

    # Comprobacion despues de la ultima iteracion
    for k in range(3):
        if abs(z - raices[k]) < tol:
            return k + 1

    return 0


# Definiendo las tres raices
raiz1 = complex(1.0, 0.0)
raiz2 = complex(-0.5, math.sqrt(3.0)/2.0)
raiz3 = complex(-0.5, -math.sqrt(3.0)/2.0)

raices = [raiz1, raiz2, raiz3]


# Construyendo la malla de 400 x 400 puntos
N = 400

# Matriz para guardar los resultados (x, y, k)
D = np.zeros((N*N, 3), dtype=float)

# Listas de puntos correspondientes a cada raiz
x_rojos = []
y_rojos = []

x_verdes = []
y_verdes = []

x_azules = []
y_azules = []

x_no_converge = []
y_no_converge = []


# Recorriendo la malla
fila = 0

for i in range(N):

    x = -1.5 + i*3.0/(N - 1)

    for j in range(N):

        y = -1.5 + j*3.0/(N - 1)

        # Construyendo el numero complejo
        z0 = complex(x, y)

        # Aplicando Newton
        codigo = clasificar_newton(z0, raices)

        # Guardando los resultados
        D[fila][0] = x
        D[fila][1] = y
        D[fila][2] = codigo

        fila += 1

        # Clasificando los puntos
        if codigo == 1:
            x_rojos.append(x)
            y_rojos.append(y)

        elif codigo == 2:
            x_verdes.append(x)
            y_verdes.append(y)

        elif codigo == 3:
            x_azules.append(x)
            y_azules.append(y)

        else:
            x_no_converge.append(x)
            y_no_converge.append(y)


# a) Guardando los datos
np.savetxt(
    "newton.txt",
    D,
    fmt=["%.8f", "%.8f", "%d"]
)


# b) Graficando el fractal
plt.figure(figsize=(8, 8))

plt.scatter(
    x_rojos, y_rojos,
    s=0.25, color="red", marker="s", linewidths=0
)

plt.scatter(
    x_verdes, y_verdes,
    s=0.25, color="green", marker="s", linewidths=0
)

plt.scatter(
    x_azules, y_azules,
    s=0.25, color="blue", marker="s", linewidths=0
)

# Puntos que no convergieron
if len(x_no_converge) > 0:
    plt.scatter(
        x_no_converge, y_no_converge,
        s=0.25, color="black", marker="s", linewidths=0
    )

plt.title("Fractal de Newton para z^3 - 1")
plt.xlabel("Parte real de z0")
plt.ylabel("Parte imaginaria de z0")

plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)

plt.gca().set_aspect("equal", adjustable="box")

plt.savefig("problema3_fractal_newton.png", dpi=170)
plt.show()


print("\nPROBLEMA 3")
print("Puntos procesados:", fila)
print("Convergen a z1:", len(x_rojos))
print("Convergen a z2:", len(x_verdes))
print("Convergen a z3:", len(x_azules))
print("No convergen:", len(x_no_converge))

##########################################
# PROBLEMA 4: IMPACTO DE UNA PELOTA

# Posicion horizontal de la pelota
def x_pelota(angulo, t):
    theta = angulo*math.pi/180.0

    return 37.0*math.cos(theta)*t


# Posicion vertical de la pelota
def y_pelota(angulo, t):
    theta = angulo*math.pi/180.0

    return 37.0*math.sin(theta)*t - 0.5*9.81*t*t


# Altura de la colina
def altura_colina(x):
    return 0.15*x*x*math.exp(-0.04*x)


# Funcion cuya raiz determina el impacto
def f_impacto(t, angulo):
    x = x_pelota(angulo, t)

    return y_pelota(angulo, t) - altura_colina(x)


# Metodo de biseccion
def biseccion_impacto(angulo, xl, xu,
                      tolerancia=0.05, max_iter=100):

    fl = f_impacto(xl, angulo)
    fu = f_impacto(xu, angulo)

    # Verificando el intervalo inicial
    if fl*fu > 0:
        raise ValueError("El intervalo no tiene cambio de signo")

    if fl == 0:
        return xl, 0, 0.0

    if fu == 0:
        return xu, 0, 0.0

    xr_anterior = None
    error = float("inf")

    # Comenzando las iteraciones
    for k in range(max_iter):

        # Calculamos el punto medio
        xr = (xl + xu)/2.0

        # Evaluamos la funcion
        fr = f_impacto(xr, angulo)

        # Calculamos el error porcentual
        if xr_anterior is not None:
            error = abs((xr - xr_anterior)/xr)*100.0

        # Comprobando raiz exacta
        if fr == 0.0:
            return xr, k + 1, 0.0

        # Actualizando el intervalo
        if fl*fr < 0:
            xu = xr

        else:
            xl = xr
            fl = fr

        # Comprobando la convergencia
        if error < tolerancia:
            return xr, k + 1, error

        xr_anterior = xr

    raise ValueError("La biseccion no alcanzo la tolerancia")


# a) Calculando el impacto para theta = 67.5 grados

angulo_inicial = 67.5

t_impacto, iteraciones, error = biseccion_impacto(
    angulo_inicial,
    5.0,
    6.0,
    tolerancia=0.05
)

x_impacto = x_pelota(angulo_inicial, t_impacto)

print("\nPROBLEMA 4")

print("Tiempo de impacto:", t_impacto)
print("Posicion de impacto:", x_impacto)
print("Iteraciones:", iteraciones)
print("Error aproximado porcentual:", error)

# Encontrando el primer impacto para un angulo
def encontrar_primer_impacto(angulo, tolerancia=0.05):

    theta = angulo*math.pi/180.0

    # Tiempo maximo de vuelo sin considerar la colina
    t_max = 2.0*37.0*math.sin(theta)/9.81

    # Evitando la raiz trivial t = 0
    t_anterior = 1e-6

    f_anterior = f_impacto(t_anterior, angulo)

    pasos = 1500

    # Buscando el primer cambio de signo
    for i in range(1, pasos + 1):

        t_actual = 1e-6 + i*(t_max - 1e-6)/pasos

        f_actual = f_impacto(t_actual, angulo)

        if f_anterior*f_actual <= 0:

            tiempo, iteracion, error = biseccion_impacto(
                angulo,
                t_anterior,
                t_actual,
                tolerancia
            )

            return tiempo, iteracion, error

        t_anterior = t_actual
        f_anterior = f_actual

    raise ValueError("No se encontro impacto")


# b) Graficando x_imp frente a theta

cantidad_angulos = 21

valores_theta = np.zeros(cantidad_angulos)
valores_x = np.zeros(cantidad_angulos)

for i in range(cantidad_angulos):

    theta = 50.0 + i

    # Calculando el tiempo para este angulo
    t, iteracion, error = encontrar_primer_impacto(theta)

    # Calculando la posicion de impacto
    x_imp = x_pelota(theta, t)

    # Guardando los valores
    valores_theta[i] = theta
    valores_x[i] = x_imp


plt.figure(figsize=(10, 6))

plt.plot(
    valores_theta,
    valores_x,
    "o-",
    color="black"
)

plt.title("Posicion de impacto vs Angulo")
plt.xlabel("Angulo (grados)")
plt.ylabel("X_imp (m)")
plt.grid(True)

plt.savefig("problema4_impacto_vs_angulo.png")
plt.show()