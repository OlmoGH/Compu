import numpy as np
import matplotlib.pyplot as plt

parametros_simulacion = np.loadtxt("C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Parametros.txt", usecols=1)
N = int(parametros_simulacion[0])
L = int(parametros_simulacion[1])
pasos = int(parametros_simulacion[2])

nombre = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Variacion_posicion/Var_posicion_{N}_{pasos}_{L}.csv"
datos = np.loadtxt(nombre, dtype=float, skiprows=1)
n_datos_por_bloque = 1
n_bloques = np.size(datos) // n_datos_por_bloque
datos = datos.reshape(n_bloques, n_datos_por_bloque).mean(axis=1)
plt.plot(np.sqrt(datos))
plt.show()