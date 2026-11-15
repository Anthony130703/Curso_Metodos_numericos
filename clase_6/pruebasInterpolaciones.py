from modelos import InterpolacionMonomios, InterpolacionLagrange, InterpolacionNewton, puntosChebyshev, graficar_comparacion

x_test = [-2.0, 0.0, 1.0]
y_test = [-27.0, -1.0, 0.0]

def pruebaMonomios():
    polinomio = InterpolacionMonomios(x_test, y_test)

    print(polinomio.lambdas)
    print(polinomio.evaluar(1.0))
    polinomio.graficar()

def pruebaLagrange():
    polinomio = InterpolacionLagrange([-2.0, 0.0, 1.0], [-27.0, -1.0, 0.0])
    print("Evaluación en t=1.0:", polinomio.evaluar(1.0))
    polinomio.graficar()
    
def pruebaNewton():
    polinomio = InterpolacionNewton([-2.0, 0.0, 1.0], [-27.0, -1.0, 0.0])
    print(polinomio.lambdas)
    print("Evaluación en t=1.0:", polinomio.evaluar(1.0))
    polinomio.graficar()

def pruebaChebyshev():

    def funcion_runge(x):
        return 1.0 / (1.0 + 25.0 * (x ** 2))

    a = -1.0
    b = 1.0
    n = 10 

    x_cheb = puntosChebyshev(a, b, n)

    y_cheb = []
    for i in range(len(x_cheb)):
        y_cheb.append(funcion_runge(x_cheb[i]))

    graficar_comparacion(x_cheb, y_cheb, ver_mono=False, ver_lagr=False, ver_newt=True)
    
#pruebaMonomios()
#pruebaLagrange()
#pruebaNewton()
pruebaChebyshev()