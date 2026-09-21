import matplotlib.pyplot as plt

# Haciendo el calculo a mano se obtuvo que para un n=4 el valor es de 0.000000448974, el cual es menor a la cota
# establecida de 1e-6

# Definimos los pasos de la deducción mezclando texto normal y fórmulas LaTeX (entre $)
pasos = [
    r"1. Definición del error (Resto) como la diferencia entre la función y el polinomio:",
    r"$R_n(x) = \sin(x) - P_n(x)$",
    r"",
    r"2. Polinomio de Taylor para el seno truncado en el término $n$ (grado $2n+1$):",
    r"$P_n(x) = x - \frac{x^3}{3!} + \frac{x^5}{5!} - \dots + \frac{(-1)^n x^{2n+1}}{(2n+1)!}$",
    r"",
    r"3. Al ser la derivada par igual a 0, evaluamos el resto en la derivada $2n+3$:",
    r"$R_n(x) = \frac{f^{(2n+3)}(\xi)}{(2n+3)!} x^{2n+3}$",
    r"",
    r"4. Tomamos el valor absoluto en ambos lados para acotar el error máximo:",
    r"$|R_n(x)| = \left| \frac{f^{(2n+3)}(\xi)}{(2n+3)!} x^{2n+3} \right|$",
    r"",
    r"5. Como $f(x) = \sin(x)$, sus derivadas siempre están acotadas por 1:",
    r"$|f^{(2n+3)}(\xi)| \leq 1$",
    r"",
    r"6. Sustituyendo este valor máximo, obtenemos la desigualdad final:",
    r"$|R_n(x)| \leq \frac{|x|^{2n+3}}{(2n+3)!}$"
]

# Configuramos el lienzo de la figura
plt.figure(figsize=(10, 8))

# Encabezado formal con tus datos para la presentación
titulo = "Deducción del Resto de Lagrange - Serie de Taylor para sin(x)"
plt.title(titulo, fontsize=12, loc='left', fontweight='bold', pad=20)

# Imprimimos línea por línea controlando la posición Y
posicion_y = 0.95
for linea in pasos:
    if linea.startswith("$"):
        # Las ecuaciones (empiezan con $) van centradas y con fuente más grande
        plt.text(0.5, posicion_y, linea, fontsize=18, ha='center', va='center')
        posicion_y -= 0.08  # Damos más espacio después de una ecuación
    else:
        # El texto descriptivo va alineado a la izquierda
        plt.text(0.05, posicion_y, linea, fontsize=11, ha='left', va='center')
        posicion_y -= 0.04  # Espacio estándar para texto

# Apagamos los ejes para que parezca una hoja limpia
plt.axis('off')

# Ajustamos márgenes y mostramos la ventana emergente
plt.tight_layout()
plt.show()