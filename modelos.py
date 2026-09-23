#Formula para el mapa logisitico
def X(n, x0, r):
    Datos = [x0]
        
    for i in range(n):
        siguiente = r * (Datos[-1] * (1 - Datos[-1]))
        Datos.append(siguiente)
    return Datos

#Formula para calcular el factorial de un numero
def factorial(n: int) -> int:
    resultado = 1
    # Multiplica secuencialmente desde 2 hasta n
    for i in range(2, n + 1):
        resultado *= i
    return resultado

#Metodo de la biseccion
def biseccion(funcion, xl:float, xu:float, tolerancia, max_iter = 100):
    #Para verificar si los puntos escogidos estan bien
    if funcion(xl) * funcion(xu) > 0:
        print("Error: El intervalo inicial no garantiza una raíz.")
        return None, None, None

    #Declarando las varibales a ultizar
    iteracion:int = 0
    ea:float = 1.0
    xr_old:float = 0.0
    errores_historial = []

    while ea > tolerancia and iteracion < max_iter:
        #calculando el punto medio del intervalo
        xr = (xl + xu)/2.0
        
        #Obteniendo el error relativo
        if iteracion > 0:
            ea = abs((xr - xr_old)/xr)
            errores_historial.append(ea)
        
        #determinando en que lado del intervalo se encuentra la raiz
        test = funcion(xl) * funcion(xr)
        
        if test < 0:
            xu = xr
        elif test > 0:
            xl = xr
        else: 
            ea = 0.0
            
        xr_old = xr
        iteracion += 1
    
    return xr, errores_historial, iteracion

#Metodo de la falsa posicion
def falsaPosicion(funcion, xl:float, xu:float, tolerancia, max_iter = 100):
    #Para verificar si los puntos escogidos estan bien
    if funcion(xl) * funcion(xu) > 0:
        print("Error: El intervalo inicial no garantiza una raíz.")
        return None, None, None

    #Declarando las varibales a ultizar
    iteracion:int = 0
    ea:float = 1.0
    xr_old:float = 0.0
    errores_historial = []

    while ea > tolerancia and iteracion < max_iter:
        #para reducir calculo inecesario
        fl = funcion(xl)
        fu = funcion(xu)

        # Nueva fórmula para hallar xr (la recta secante)
        xr = xu - (fu * (xl - xu)) / (fl - fu)
        
        #Obteniendo el error relativo
        if iteracion > 0:
            ea = abs((xr - xr_old)/xr)
            errores_historial.append(ea)
        
        #determinando en que lado del intervalo se encuentra la raiz
        test = funcion(xl) * funcion(xr)
        
        if test < 0:
            xu = xr
        elif test > 0:
            xl = xr
        else: 
            ea = 0.0
            
        xr_old = xr
        iteracion += 1
    
    return xr, errores_historial, iteracion

#Metodo del punto fijo
def puntoFijo(funcion, x0:float , tolerancia ,max_iter = 100):
    #Variables a utilizar
    iteracion:int = 0
    ea:float = 1.0
    xr:float = x0
    errores_historial = []

    #implementando un bucle para hcer el metodo
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

    #implementando un bucle para hcer el metodo
    while ea > tolerancia and iteracion < max_iter:
        #Guardando el primer punto para luego medir el error
        xr_old = xr 
        df = Dfuncion(xr_old)

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