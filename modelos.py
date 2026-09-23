import numpy as np

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
            ea = np.abs((xr - xr_old)/xr)
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
        