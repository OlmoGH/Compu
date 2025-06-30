import numpy as np
import matplotlib.pyplot as plt

parametros_simulacion = np.loadtxt("C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Parametros.txt", usecols=1)
N = int(parametros_simulacion[0])
L = int(parametros_simulacion[1])
pasos = int(parametros_simulacion[2])

nombre = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Variacion_posicion/Separacion_{N}_{pasos}_{L}.csv"
datos = np.loadtxt(nombre, dtype=float, skiprows=1)

# El número de datos por bloque nos dice los datos que sobran para poder dividir el array en partes iguales
incremento = 300000
n_datos_por_bloque = 1
datos_sobrantes = np.size(datos) % n_datos_por_bloque
n_datos_utiles = np.size(datos) - datos_sobrantes

n_bloques = n_datos_utiles // n_datos_por_bloque

# Eliminamos los datos sobrantes
if datos_sobrantes != 0 :
    datos = datos[:-datos_sobrantes]
datos = datos.reshape(n_bloques, n_datos_por_bloque).mean(axis=1)
plt.plot(np.sqrt(datos))
for position in np.arange(0, np.size(datos)+1, incremento / n_datos_por_bloque) :
    plt.axvline(position, color='r', linestyle='-')
plt.show()