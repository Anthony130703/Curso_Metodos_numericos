#parametros a utilizar
a = 2.0
N = 100
z_inicial = 0
eje_x = []
eje_y = []
max_interaciones = 100

#usando diccionario para organizar los puntos a graficar

malla_puntos = {
    "validos"
}

paso = (a - (-a))/(N-1)

def Z(z:complex, c:complex):
    z_actual = ((abs(z))**2) + c
    return z_actual

for i in range(N):
    valor = -a + (i * paso)
    eje_x.append(valor)
    eje_y.append(valor)
    
for x in eje_x:
    for y in eje_y:
        c = complex(x, y)
        z_actual = z_inicial
        pertenece_al_conjunto = True
        
        for interaciones in range(max_interaciones):
            z_actual = Z(z_actual, c)
            