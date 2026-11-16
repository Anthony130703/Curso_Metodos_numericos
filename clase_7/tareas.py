import math
import matplotlib.pyplot as plt
from modelos import RegresionLineal, RegresionNoLineal

def ejercicio1_hubble():
    # Datos de la Diapositiva 61
    d = [0.032, 0.034, 0.214, 0.263, 0.275, 0.275, 0.45, 0.5, 0.5, 0.63, 0.8, 0.9, 
         0.9, 0.9, 0.9, 1.0, 1.1, 1.1, 1.4, 1.7, 2.0, 2.0, 2.0, 2.0]
    v = [170, 290, -130, -70, -185, -220, 200, 290, 270, 200, 300, -30, 
         650, 150, 500, 920, 450, 500, 500, 960, 500, 850, 800, 1090]
    n = len(d)

    # --- 1. Usar la clase para el Modelo Clásico (v = a0 + a1*d) ---
    modelo_clasico = RegresionLineal(d, v)
    
    # --- 2. Calcular a mano el Modelo Sin Intercepto (v = H0*d) ---
    sum_dv = 0.0
    sum_d2 = 0.0
    for i in range(n):
        sum_dv = sum_dv + (d[i] * v[i])
        sum_d2 = sum_d2 + (d[i] * d[i])
        
    h0 = sum_dv / sum_d2  # Fórmula teórica
    
    # Error estándar para el modelo sin intercepto (n-1 grados de libertad)
    sr_h0 = 0.0
    for i in range(n):
        v_estimado = h0 * d[i]
        sr_h0 = sr_h0 + (v[i] - v_estimado)**2
    syx_h0 = math.sqrt(sr_h0 / (n - 1))
    
    # --- Reportes en consola ---
    print("--- Modelo 1: Regresión Clásica (v = a0 + a1*d) ---")
    print(f"Intercepto (a0): {modelo_clasico.a0:.4f}")
    print(f"Pendiente (a1): {modelo_clasico.a1:.4f}")
    print(f"Coef. Determinación (r2): {modelo_clasico.r2:.4f}")
    print(f"Error Estándar (Sy/x): {modelo_clasico.syx:.4f}\n")
    
    print("--- Modelo 2: Sin Intercepto (v = H0*d) ---")
    print(f"Constante de Hubble (H0): {h0:.4f} km/s/Mpc")
    print(f"Error Estándar (Sy/x): {syx_h0:.4f}\n")
    
    # --- 4. Cálculo de la Edad del Universo ---
    km_por_Mpc = 3.086e19
    segundos_por_año = 365 * 24 * 3600
    t_H_segundos = (1.0 / h0) * km_por_Mpc
    t_H_años = t_H_segundos / segundos_por_año
    
    print("--- Análisis Físico ---")
    print(f"Edad característica (t_H): {t_H_años:.2e} años")
    
    # --- 5. Gráficas ---
    d_grafica = [0.0, 2.2]
    v_modelo1 = [modelo_clasico.evaluar(x) for x in d_grafica]
    v_modelo2 = [h0 * x for x in d_grafica]
    
    residuos1 = [v[i] - modelo_clasico.evaluar(d[i]) for i in range(n)]
    residuos2 = [v[i] - (h0 * d[i]) for i in range(n)]
    
    plt.figure(figsize=(12, 5))
    
    plt.subplot(1, 2, 1)
    plt.plot(d, v, 'ko', label="Datos de Hubble")
    plt.plot(d_grafica, v_modelo1, 'r-', label="Modelo a0 + a1*d")
    plt.plot(d_grafica, v_modelo2, 'b--', label="Modelo H0*d")
    plt.xlabel("Distancia d (Mpc)")
    plt.ylabel("Velocidad v (km/s)")
    plt.title("Ajuste de la Ley de Hubble")
    plt.legend()
    plt.grid(True)
    
    plt.subplot(1, 2, 2)
    plt.plot(d, residuos1, 'rs', label="Residuos Mod. 1")
    plt.plot(d, residuos2, 'b^', label="Residuos Mod. 2")
    plt.axhline(0, color='k', linestyle='--')
    plt.xlabel("Distancia d (Mpc)")
    plt.ylabel("Residuo")
    plt.title("Análisis de Residuos")
    plt.legend()
    plt.grid(True)
    
    plt.tight_layout()
    plt.show()


def ejercicio2_kepler():
    a_datos = [181.4, 421.7, 671.0, 1070.4, 1882.7, 11461.0]
    t_datos = [0.4982, 1.7691, 3.5512, 7.1546, 16.689, 250.56]
    n = len(a_datos)
    g = 6.674e-11
    pi = math.pi
    m_ref = 1.898e27

    ln_a = []
    ln_t = []
    for i in range(n):
        ln_a.append(math.log(a_datos[i]))
        ln_t.append(math.log(t_datos[i]))

    modelo1 = RegresionLineal(ln_a, ln_t)
    p = modelo1.a1
    r2 = modelo1.r2

    print("--- ITEM 1 ---")
    print("Potencia p:", p)
    print("R2:", r2)
    print("Error p (%):", abs(p - 1.5) / 1.5 * 100)

    x_si = []
    y_si = []
    for i in range(n):
        x_si.append((a_datos[i] * 1000000) ** 3)
        y_si.append((t_datos[i] * 86400) ** 2)

    sum_xy = 0.0
    sum_x2 = 0.0
    for i in range(n):
        sum_xy = sum_xy + x_si[i] * y_si[i]
        sum_x2 = sum_x2 + x_si[i] * x_si[i]

    k_lineal = sum_xy / sum_x2
    m_lineal = (4 * pi * pi) / (g * k_lineal)

    print("\n--- ITEM 2 ---")
    print("Masa Jupiter (Lineal):", m_lineal)

    ln_t2 = []
    ln_a3 = []
    for i in range(n):
        ln_t2.append(math.log(y_si[i]))
        ln_a3.append(math.log(a_datos[i] * 1000000))

    suma_a0 = 0.0
    for i in range(n):
        suma_a0 = suma_a0 + (ln_t2[i] - 3.0 * ln_a3[i])

    a0_log = suma_a0 / n
    k_log = math.exp(a0_log)
    m_log = (4 * pi * pi) / (g * k_log)

    print("\n--- ITEM 3 ---")
    print("Masa Jupiter (Log):", m_log)

    residuos_lin = []
    residuos_log = []
    lunas = [1, 2, 3, 4, 5, 6]

    for i in range(n):
        est_lin = k_lineal * x_si[i]
        res_lin = (y_si[i] - est_lin) / y_si[i]
        residuos_lin.append(res_lin)
        
        est_log = k_log * x_si[i]
        res_log = (y_si[i] - est_log) / y_si[i]
        residuos_log.append(res_log)

    print("\n--- ITEM 4 ---")
    for i in range(n):
        k_ind = y_si[i] / x_si[i]
        m_ind = (4 * pi * pi) / (g * k_ind)
        err_ind = abs(m_ind - m_ref) / m_ref * 100
        print("Luna", i+1, "Masa:", m_ind, "Error (%):", err_ind)

    plt.figure(figsize=(8, 5))
    plt.plot(lunas, residuos_lin, 'rs-', label="Lineal")
    plt.plot(lunas, residuos_log, 'b^-', label="Log")
    plt.axhline(0, color='k')
    plt.title("Comparación de Residuos Relativos (Tercera Ley de Kepler)")
    plt.xlabel("Lunas (1=Amaltea, 6=Himalia)")
    plt.ylabel("Residuo Relativo")
    plt.grid(True)
    plt.legend()
    plt.show()

def ejercicio3_radiactividad():
    t_datos = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 
               5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 
               10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0]
    n_datos = [1626, 1444, 1217, 1117, 935, 815, 736, 625, 578, 501, 
               452, 393, 359, 302, 293, 234, 219, 202, 155, 159, 
               139, 107, 113, 120, 103, 98, 84, 78, 100, 75]
    n = len(t_datos)

    ln_n = []
    for i in range(n):
        ln_n.append(math.log(n_datos[i]))

    modelo_lin = RegresionLineal(t_datos, ln_n)
    
    lam_lin = -modelo_lin.a1
    n0_lin = math.exp(modelo_lin.a0)
    t_half_lin = math.log(2) / lam_lin

    print("--- ITEM 1: Regresion Lineal (Ignorando Fondo B) ---")
    print("N0:", n0_lin)
    print("Lambda:", lam_lin)
    print("t_1/2 (min):", t_half_lin)

    suma_b = 0.0
    for i in range(n - 5, n):
        suma_b = suma_b + n_datos[i]
    b_ini = suma_b / 5.0

    param_ini = [n0_lin, lam_lin, b_ini]

    def funcion_modelo(t, p):
        return p[0] * math.exp(-p[1] * t) + p[2]

    def funcion_derivadas(t, p):
        df_a0 = math.exp(-p[1] * t)
        df_a1 = -p[0] * t * math.exp(-p[1] * t)
        df_a2 = 1.0
        return [df_a0, df_a1, df_a2]

    modelo_nl = RegresionNoLineal(t_datos, n_datos, param_ini, funcion_modelo, funcion_derivadas)

    n0_nl = modelo_nl.parametros[0]
    lam_nl = modelo_nl.parametros[1]
    b_nl = modelo_nl.parametros[2]
    t_half_nl = math.log(2) / lam_nl

    print("\n--- ITEM 2 y 3: Regresion No Lineal (Gauss-Newton) ---")
    print("N0:", n0_nl)
    print("Lambda:", lam_nl)
    print("Fondo B:", b_nl)
    print("t_1/2 (min):", t_half_nl)
    
    error_teorico = abs(t_half_nl - 2.552) / 2.552 * 100
    print("Error respecto a 2.552 min (%):", error_teorico)

#ejercicio1_hubble()
#ejercicio2_kepler()
ejercicio3_radiactividad()