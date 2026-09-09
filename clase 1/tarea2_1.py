import matplotlib.pyplot as plt

#parametros a utilizar
a = 2.0
N = 100
z_inicial = 0.0
eje_x = []
eje_y = []
max_interaciones = 100

#usando diccionario para organizar los puntos a graficar

malla_puntos = {
    "validos":{
        "x": [],
        "y": []
    },
    "no_validos":{
        "x": [],
        "y": []
    }
}

paso = (a - (-a))/(N-1)

def Z(z:complex, c:complex):
    z_actual = (z**2) + c
    return z_actual

for i in range(N):
    valor = -a + (i * paso)
    eje_x.append(valor)
    eje_y.append(valor)
    
for x in eje_x:
    for y in eje_y:
        c = complex(x, y)
        z_actual = z_inicial
        iteracion = 0
        
        # El bucle continúa mientras no alcancemos el máximo de iteraciones
        # Y mientras la magnitud de z_actual no supere 2
        while iteracion < max_interaciones and abs(z_actual) <= 2:
            z_actual = Z(z_actual, c)
            iteracion += 1
            
        # Si al terminar el bucle z no escapó (su magnitud es <= 2), pertenece al conjunto
        if abs(z_actual) <= 2:
            malla_puntos["validos"]["x"].append(x)
            malla_puntos["validos"]["y"].append(y)
        else:
            malla_puntos["no_validos"]["x"].append(x)
            malla_puntos["no_validos"]["y"].append(y)
            
plt.figure(figsize=(8, 8))
plt.scatter(malla_puntos["no_validos"]["x"], malla_puntos["no_validos"]["y"], s=1, color="white")
plt.scatter(malla_puntos["validos"]["x"], malla_puntos["validos"]["y"], s=1, color="black")
plt.title("Visualización del Conjunto de Mandelbrot")
plt.xlabel("Parte Real (x)")
plt.ylabel("Parte Imaginaria (y)")
plt.show()