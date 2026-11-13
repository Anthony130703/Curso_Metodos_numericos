from modelos import InterpolacionMonomios, InterpolacionLagrange, InterpolacionNewton
x_test = [-2.0, 0.0, 1.0]
y_test = [-27.0, -1.0, 0.0]

def pruebaMonomios():
    polinomio = InterpolacionMonomios(x_test, y_test)

    print(polinomio.lambdas)
    print(polinomio.evaluar(1.0))

def pruebaLagrange():
    polinomio = InterpolacionLagrange([-2.0, 0.0, 1.0], [-27.0, -1.0, 0.0])
    print("Evaluación en t=1.0:", polinomio.evaluar(1.0))
    
def pruebaNewton():
    polinomio = InterpolacionNewton([-2.0, 0.0, 1.0], [-27.0, -1.0, 0.0])
    print(polinomio.lambdas)
    print("Evaluación en t=1.0:", polinomio.evaluar(1.0))
    
#pruebaMonomios()
#pruebaLagrange()
pruebaNewton()