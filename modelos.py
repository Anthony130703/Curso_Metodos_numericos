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