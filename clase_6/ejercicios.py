from modelos import InterpolacionMonomios, InterpolacionLagrange, InterpolacionNewton, puntosChebyshev, InterpolacionSplineCubica
import numpy as np
import matplotlib.pyplot as plt

def ejercicio1():
    x = [0.0, 1.0, 2.0, 5.5, 11.0, 13.0, 16.0, 18.0]
    y = [0.5, 3.134, 5.3, 9.9, 10.2, 9.35, 7.2, 6.2]

    metodo_mono = InterpolacionMonomios(x, y)
    metodo_lagr = InterpolacionLagrange(x, y)
    metodo_newt = InterpolacionNewton(x, y)

    valor_mono = metodo_mono.evaluar(8.0)
    valor_lagr = metodo_lagr.evaluar(8.0)
    valor_newt = metodo_newt.evaluar(8.0)

    print(valor_mono)
    print(valor_lagr)
    print(valor_newt)

def funcion(x):
    parte1 = np.tanh(20.0 * np.sin(12.0 * x))
    parte2 = (1.0 / 50.0) * np.exp(3.0 * x) * np.sin(300.0 * x)
    return parte1 + parte2

def ejercicio2():

    a = 0.0
    b = 1.0
    n = 99

    #Creando las listas
    x_equi = list(np.linspace(a, b, n+1))
    y_equi = []

    #Completando las lista y
    for i in range (len(x_equi)):
        y_equi.append(funcion(x_equi[i]))

    #Puntos de Chebyshev
    x_cheb = puntosChebyshev(a, b, n)
    y_cheb = []
    for i in range(len(x_cheb)):
        y_cheb.append(funcion(x_cheb[i]))

    lagrange_equi = InterpolacionLagrange(x_equi, y_equi)
    lagrange_cheb = InterpolacionLagrange(x_cheb, y_cheb)

    #Calculo de errores
    x_eval = list(np.linspace(a, b, 5000))
    error_equi = []
    error_cheb = []

    for i in range(len(x_eval)):
        val_real = funcion(x_eval[i])
        val_equid = lagrange_equi.evaluar(x_eval[i])
        error_equi.append(abs(val_real - val_equid))
        val_cheb = lagrange_cheb.evaluar(x_eval[i])
        error_cheb.append(abs(val_real - val_cheb))

    plt.figure(figsize=(10, 5))
    plt.title("Error de Interpolación: Equidistante vs Chebyshev")
    plt.xlabel("Eje X (Intervalo [0, 1])")
    plt.ylabel("Error Absoluto |f(x) - P(x)|")
    plt.plot(x_eval, error_equi, 'r-', label="Error Equidistante")
    plt.plot(x_eval, error_cheb, 'b-', label="Error Chebyshev")
    plt.yscale('log')
    plt.grid(True)
    plt.legend()
    plt.show()

def ejercico2Interpolacion():
    a = 0.0
    b = 1.0
    valores_n = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    errores_maximos = []

    for n in valores_n:
        x_cheb = puntosChebyshev(a, b, n - 1)
        y_cheb = []
        for i in range(len(x_cheb)):
            y_cheb.append(funcion(x_cheb[i]))

        lagrangeCheb = InterpolacionLagrange(x_cheb, y_cheb)

        x_eval = list(np.linspace(a, b, 200))
        error_max = 0.0

        for i in range(len(x_eval)):
            val_real = funcion(x_eval[i])
            val_cheb = lagrangeCheb.evaluar(x_eval[i])
            errorActual = abs(val_real - val_cheb)
            if errorActual > error_max:
                error_max = errorActual

        errores_maximos.append(error_max)

    plt.figure(figsize=(10, 5))
    plt.plot(valores_n, errores_maximos, 'bo-')
    plt.title("Evolución del Error Máximo (Chebyshev)")
    plt.xlabel("Grado del Polinomio (n)")
    plt.ylabel("Error Máximo")
    plt.yscale('log')
    plt.grid(True)
    plt.show()

def verificar_malla():
    a = 0.0
    b = 1.0
    n = 100                                   # prueba con un n, no hace falta todos
    x_cheb = puntosChebyshev(a, b, n - 1)
    y_cheb = [funcion(xi) for xi in x_cheb]
    lag = InterpolacionLagrange(x_cheb, y_cheb)

    for m in [200, 500, 2000, 5000, 10000]:
        x_eval = np.linspace(a, b, m)
        p = np.array([lag.evaluar(t) for t in x_eval])
        print(m, np.max(np.abs(funcion(x_eval) - p)))

def funcion_seno_cuadrado(t):
    return (np.sin(t))**2

def ejercicio3():
    a = 0.0
    b = 2.0 * np.pi
    n_puntos = 8
    grado = n_puntos - 1

    t_equi = list(np.linspace(a, b, n_puntos))
    y_equi = []
    for i in range(len(t_equi)):
        y_equi.append(funcion_seno_cuadrado(t_equi[i]))

    t_cheb = puntosChebyshev(a, b, grado)
    y_cheb = []
    for i in range(len(t_cheb)):
        y_cheb.append(funcion_seno_cuadrado(t_cheb[i]))

    poli_equi = InterpolacionLagrange(t_equi, y_equi)
    spline = InterpolacionSplineCubica(t_equi, y_equi)
    poli_cheb = InterpolacionLagrange(t_cheb, y_cheb)

    t_eval = list(np.linspace(a, b, 200))
    y_real = []
    y_eval_equi = []
    y_eval_spline = []
    y_eval_cheb = []

    for i in range(len(t_eval)):
        y_real.append(funcion_seno_cuadrado(t_eval[i]))
        y_eval_equi.append(poli_equi.evaluar(t_eval[i]))
        y_eval_spline.append(spline.evaluar(t_eval[i]))
        y_eval_cheb.append(poli_cheb.evaluar(t_eval[i]))

    plt.figure(figsize=(12, 6))
    plt.plot(t_eval, y_real, 'k-', linewidth=2, label="Real: sen(t)^2")
    plt.plot(t_eval, y_eval_equi, 'r--', label="Polinomio Grado 7 (Equidistantes)")
    plt.plot(t_eval, y_eval_spline, 'm-', label="Spline Cubica (Equidistantes)")
    plt.plot(t_eval, y_eval_cheb, 'b-.', label="Polinomio Grado 7 (Chebyshev)")
    
    plt.plot(t_equi, y_equi, 'ro', markersize=6)
    plt.plot(t_cheb, y_cheb, 'bs', markersize=6)
    
    plt.title("Comparacion de Metodos Numéricos")
    plt.xlabel("Tiempo (t)")
    plt.ylabel("f(t) = sen(t)^2")
    plt.grid(True)
    plt.legend()
    plt.show()

#ejercicio1()
#ejercicio2()
#ejercico2Interpolacion()
ejercicio3()
#verificar_malla()
