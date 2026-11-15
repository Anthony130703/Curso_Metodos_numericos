from modelos import InterpolacionMonomios, InterpolacionLagrange, InterpolacionNewton, puntosChebyshev
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
    x_eval = list(np.linspace(a, b, 500))
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
    plt.yscale('log')
    plt.grid(True)
    plt.show()



#ejercicio1()
ejercicio2()
#ejercico2Interpolacion()