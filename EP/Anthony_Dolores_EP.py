import numpy as np
import random
import matplotlib.pyplot as plt
import cmath

#sustitucion backward
def sustitucion_backward(U, b):
    ##Resolviendo Ux = b para U triangular superior
    n = len(b)
    tipo = complex if np.iscomplexobj(U) or np.iscomplexobj(b) else float
    x = np.zeros(n, dtype=tipo)
    for k in range(n - 1, -1, -1):
        if abs(U[k][k]) < 1e-15:
            raise ValueError(f"Elemento diagonal nulo en fila {k}.")
        suma = 0
        for i in range (k+1, n):
            suma = suma + U[k][i] * x[i]
        x[k] = (b[k] - suma) / U[k][k]
    return x

#Metodo de gauss con pivote
def gauss_pivoteo (A, b):
    n = len(b)
    A_copy = [fila[:] for fila in A]
    b_copy = b[:]

    #Etapa 1: Eliminacion hacia delante con pivoteo parcial
    for k in range(n-1):
        #Asumiendo inicialmente que el valor maximo esta en la fila actual k
        p = k 
        max_valor = abs(A_copy[k][k])
        #Buscando si hay un valor mas grande debajo en la misma columna k
        for i in range(k+1, n):
            if abs(A_copy[i][k]) > max_valor:
                max_valor = abs(A_copy[i][k])
                p = i
        #Si el maximo no esta en la fila k, intercambiamos la fila k con la fila p
        if p != k:
            #Intercambiamos de la fila completa en A
            for j in range(k, n): #Solo necesitamos intercambiar desde la columna k en adelante
                temp = A_copy[k][j]
                A_copy[k][j] = A_copy[p][j]
                A_copy[p][j] = temp
            #Intercambiamos la columna b tambien
            temp_b = b_copy[k]
            b_copy[k] = b_copy[p]
            b_copy[p] = temp_b
        #Verificando si es que toda la columan es cero
        if abs(A_copy[k][k]) < 1e-15:
            raise ValueError("El sistema no tiene solucion única")
        #Aplicando la eliminación gausiana
        for i in range(k+1, n):
            factor = A_copy[i][k]/A_copy[k][k]
            for j in range(k+1, n):
                A_copy[i][j] = A_copy[i][j] - factor*A_copy[k][j]
            b_copy[i] = b_copy[i] - factor*b_copy[k]
    #Finalemnte usamos el metodo sustitucion backward
    return sustitucion_backward(A_copy, b_copy)


#PROBLEMA 3: lA CURVA DEL DRAGON
def ejecutar_problema_3():
    print("\nPROBLEMA DEL DRAGON")
    #Inciso a)
    print("\n Inciso a): ")
    C = [
        [0.0,  0.0,  1.0],
        [1.0, 0.0,  1.0],
        [0.5,  0.5, 1.0]
    ]
    b_1 = [0.0, 0.5, 0.0]
    b_2 = [0.0, 0.5, 0.5]
    b_3 = [1.0, 0.5, 0.5]
    b_4 = [0.0, 0.5, 0.0]
    
    u_1 = gauss_pivoteo(C, b_1)
    u_2 = gauss_pivoteo(C, b_2)
    u_3 = gauss_pivoteo(C, b_3)
    u_4 = gauss_pivoteo(C, b_4)

    print("Vector | Coeficiente | Valor")
    print("-----------------------------")
    print(f"u_1    | a_1         | {u_1[0]}")
    print(f"u_1    | b_1         | {u_1[1]}")
    print(f"u_1    | e_1         | {u_1[2]}")
    print(f"u_2    | c_1         | {u_2[0]}")
    print(f"u_2    | d_1         | {u_2[1]}")
    print(f"u_2    | f_1         | {u_2[2]}")
    print(f"u_3    | a_2         | {u_3[0]}")
    print(f"u_3    | b_2         | {u_3[1]}")
    print(f"u_3    | e_2         | {u_3[2]}")
    print(f"u_4    | c_2         | {u_4[0]}")
    print(f"u_4    | d_2         | {u_4[1]}")
    print(f"u_4    | f_2         | {u_4[2]}")
    
    print("\nEl sistema se tuvo que resolver con el metodo de gauss con pivote ya que uno de los coeficientes")
    print("de la primera columna es 0.0")
    
    #Inciso b)
    print("\nInciso b): ")
    def juego_del_caos(numero_iteraciones):
        x = 0.0
        y = 0.0
        x_w1 = []
        y_w1 = []
        x_w2 = []
        y_w2 = []
        #Realizando las iteraciones
        for iteracion in range(numero_iteraciones):
            #Obteniendo los numeros entre 0 y 1
            r = random.random()
            if r < 0.5:
                x_nuevo = u_1[0] * x + u_1[1] * y + u_1[2]
                y_nuevo = u_2[0] * x + u_2[1] * y + u_2[2]
                #para descartar las 20 primeras iteraciones
                if iteracion >= 20:
                    x_w1.append(x_nuevo)
                    y_w1.append(y_nuevo)
            else:
                x_nuevo = u_3[0] * x + u_3[1] * y + u_3[2]
                y_nuevo = u_4[0] * x + u_4[1] * y + u_4[2]
                #para descartar las 20 primeras iteraciones
                if iteracion >= 20:
                    x_w2.append(x_nuevo)
                    y_w2.append(y_nuevo)

            x = x_nuevo
            y = y_nuevo

        return x_w1, y_w1, x_w2, y_w2

    x_5000_w1, y_5000_w1, x_5000_w2, y_5000_w2 = juego_del_caos(5000)
    x_50000_w1, y_50000_w1, x_50000_w2, y_50000_w2 = juego_del_caos(50000)

    plt.figure(figsize=(12, 6))

    plt.subplot(1, 2, 1)
    plt.scatter(x_5000_w1, y_5000_w1, color="red", s=1, label="w1")
    plt.scatter(x_5000_w2, y_5000_w2, color="blue", s=1, label="w2")
    plt.title("Juego del caos: 5 000 iteraciones")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()

    plt.subplot(1, 2, 2)
    plt.scatter(x_50000_w1, y_50000_w1, color="red", s=1, label="w1")
    plt.scatter(x_50000_w2, y_50000_w2, color="blue", s=1, label="w2")
    plt.title("Juego del caos: 50 000 iteraciones")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()

    plt.tight_layout()
    plt.show()
    
def puntoFijo(funcion, x0:float , tolerancia ,max_iter = 100):
    #Variables a utilizar
    iteracion:int = 0
    ea:float = 1.0
    xr:float = x0
    errores_historial = []

    #implementando un bucle para hacer el metodo
    while ea > tolerancia and iteracion < max_iter:
        #Guardando el primer punto para luego medir el error
        xr_old = xr 

        #Asignando el punto siguiente mediante la funcion
        xr = funcion(xr_old)

        #Calculando el error desde la primera vuelta
        if xr != 0: #Esto lo hacemos para cuidarnos de la division por cero
            ea = abs((xr - xr_old)/ xr)
            errores_historial.append(ea)

        iteracion += 1

    return xr, errores_historial, iteracion

#Metodo de Newton-Raphson
def newtonRaphson(funcion, Dfuncion , x0:float , tolerancia ,max_iter = 100):
    #Variables a utilizar
    iteracion:int = 0
    ea:float = 1.0
    xr:float = x0
    errores_historial = []

    #implementando un bucle para hacer el metodo
    while ea > tolerancia and iteracion < max_iter:
        #Guardando el primer punto para luego medir el error
        xr_old = xr 
        df = Dfuncion(xr_old)
        
        #Verificando que la derivada no sea 0 para evitar errores
        if df == 0:
            print("Error: La derivada es cero. El método ha colapsado.")
            return None, None, None  #Abortamos sin usar break

        #Asignando el punto siguiente mediante la regla
        xr = xr_old - ((funcion(xr_old))/(df))

        #Calculando el error desde la primera vuelta
        if xr != 0: #Esto lo hacemos para cuidarnos de la division por cero
            ea = abs((xr - xr_old)/ xr)
            errores_historial.append(ea)

        iteracion += 1

    return xr, errores_historial, iteracion    

#PROBLEMA 1: LA ECUACION DE KEPLER
def ejecutar_problema_1():
    print("\nLA ECUACION DE KEPLER")
    #Inciso a)
    print("\nInciso a): ")
    #Definiendo la funcion a usar en el metodo
    def keplerF(x):
        return 0.1 + (0.967)*np.sin(x)
    
    def dkepler(x):
        return 1.0 - (0.967)*np.cos(x)
    
    def kepler(x):
        return x -(0.967)*np.sin(x) - 0.1
    
    
    E_0 = 0.1
    tolerancia = 1e-10
        
    raiz_pf, errores_pf, iter_pf = puntoFijo(keplerF, E_0, tolerancia)
    print("\n--- Punto Fijo ---")
    print(f"Raíz (x): {raiz_pf}")
    print(f"Iteraciones: {iter_pf}")
    if errores_pf:
        print(f"Error final: {errores_pf[-1]}")
        
    raiz_nr, errores_nr, iter_nr = newtonRaphson(kepler, dkepler, E_0, tolerancia)
    print("\n--- Newton-Raphson ---")
    print(f"Raíz (x): {raiz_nr}")
    print(f"Iteraciones: {iter_nr}")
    if errores_nr:
        print(f"Error final: {errores_nr[-1]}")
        
    print("Punto fijo necesitó más iteraciones porque converge linealmente;")
    print("Newton-Raphson converge más rápido porque su convergencia es cuadrática.")
        
    # === GRÁFICA DE CONVERGENCIA ===
    # Configuramos el tamaño de la ventana
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(errores_pf) + 1), errores_pf, marker='^', label='Punto Fijo')
    plt.plot(range(1, len(errores_nr) + 1), errores_nr, marker='D', label='Newton-Raphson')
    
    plt.yscale('log')  # Escala logarítmica para ver la caída del error claramente
    plt.title('Velocidad de Convergencia de los Métodos Numéricos', fontsize=14, fontweight='bold')
    plt.xlabel('Número de Iteraciones', fontsize=12)
    plt.ylabel('Error Relativo Aproximado (Ea)', fontsize=12)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()
    # plt.show()
    
    #Inciso b)
    print("\nInciso b): ")
    e = 0.967
    a = 17.83
    T = 75.3
    x_orbita = []
    y_orbita = []

    for i in range(1, 201):
        t_i = i * T / 200
        M_i = 2 * np.pi * t_i / T

        def kepler_orbita(E):
            return E - e * np.sin(E) - M_i

        E_i, errores_i, iteraciones_i = newtonRaphson(kepler_orbita, dkepler, np.pi, tolerancia)
        x_orbita.append(a * (np.cos(E_i) - e))
        y_orbita.append(a * np.sqrt(1 - e**2) * np.sin(E_i))

    plt.figure(figsize=(8, 6))
    plt.plot(x_orbita, y_orbita, marker="o", linestyle="-", markersize=3, label="Órbita del cometa")
    plt.plot(0, 0, marker="*", color="orange", markersize=12, label="Sol")
    plt.title("Órbita del cometa Halley")
    plt.xlabel("x (UA)")
    plt.ylabel("y (UA)")
    plt.grid(True)
    plt.legend()
    plt.show()

    print("Los puntos se acumulan cerca del afelio, lejos del Sol.")
    print("El cometa se mueve más lentamente cuando está más lejos del Sol.")
    
    
#PROBLEMA 2: EFECTO TUNEL
def ejecutar_problema_2():
    print("\nPROBLEMA 2: EFECTO TUNEL")
    #Inciso a): 
    print("\nInciso a):")

    E = 0.4
    V0 = 1.0
    a = 0.5
    C = 0.0381

    k = cmath.sqrt(E / C)
    q = cmath.sqrt((E - V0) / C)

    M = [
        [1, -1, -1, 0],
        [(-1 * k), ((-1)*q), q, 0],
        [0, cmath.exp(1j * q * a), cmath.exp(-1j * q * a), -cmath.exp(1j * k * a)],
        [0, q * cmath.exp(1j * q * a), -1* q * cmath.exp(-1j * q * a), -1 * k * cmath.exp(1j * k * a)]
    ]
    b = [-1, -1 * k, 0, 0]

    print("M =")
    for fila in M:
        print(fila)

    print("b =", b)
    
    #Inciso b:
    print("\nInciso b): ")
    u = gauss_pivoteo(M, b)
    print(u)

    r = u[0]
    t = u[3]
    R = abs(r) ** 2
    T_numerica = abs(t) ** 2
    kappa = cmath.sqrt((V0 - E) / C)
    T_analitica = 1 / ( 1 + (V0**2 * np.sinh(kappa * a)**2) / (4 * E * (V0 - E)))

    print(f"r = {r}")
    print(f"t = {t}")
    print(f"R = {R}")
    print(f"T numerica = {T_numerica}")
    print(f"R + T = {R + T_numerica}")
    print(f"T analitica = {T_analitica}")
    print("R + T = 1 significa que se conserva la probabilidad de la partícula.")
    
    
    #Inciso c:
    print("\nInciso c): ")
    energia_minima = 0.05
    energia_maxima = 10.0
    N = 400
    energias = []
    paso = (energia_maxima - energia_minima)/(N-1)
    #obteniendo una lista con los 400 valores de energia equiespaciadas
    for i in range(N):
        valor = energia_minima + (i * paso)
        energias.append(valor)
        
    T_solucion_analitica = []
    T_solucion_numerica = []
        
    for valor in energias:
        if valor < V0:
            kappa = np.sqrt((V0 - valor) / C)
            T_analitica_valor = 1 / (
                1 + (V0**2 * np.sinh(kappa * a)**2) / (4 * valor * (V0 - valor))
            )
        else:
            q_real = np.sqrt((valor - V0) / C)
            T_analitica_valor = 1 / (
                1 + (V0**2 * np.sin(q_real * a)**2) / (4 * valor * (valor - V0))
            )
        T_solucion_analitica.append(T_analitica_valor)

        k = cmath.sqrt(valor / C)
        q = cmath.sqrt((valor - V0) / C)
        M = [
                [1, -1, -1, 0],
                [(-1 * k), ((-1)*q), q, 0],
                [0, cmath.exp(1j * q * a), cmath.exp(-1j * q * a), -cmath.exp(1j * k * a)],
                [0, q * cmath.exp(1j * q * a), -1* q * cmath.exp(-1j * q * a), -1 * k * cmath.exp(1j * k * a)]
            ]
        b = [-1, -1 * k, 0, 0]
        solucion = gauss_pivoteo(M, b)
        T_numerica_valor = abs(solucion[3]) ** 2
        T_solucion_numerica.append(T_numerica_valor)
        
    plt.figure(figsize=(10, 6))
    plt.plot(energias, T_solucion_analitica, label="T(E) analítica")
    plt.scatter(energias, T_solucion_numerica, s=12, color="black",label="T(E) numérica")
    
    plt.title("Transmisión en función de la energía")
    plt.xlabel("E (eV)")
    plt.ylabel("T(E)")
    plt.grid(True)
    plt.legend()
    plt.show()
        
        
        

ejecutar_problema_1()
ejecutar_problema_2()
ejecutar_problema_3()