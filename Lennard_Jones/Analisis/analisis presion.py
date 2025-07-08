import numpy as np
import matplotlib.pyplot as plt

tamaño_bloque = 10
bloques = 100
datos = np.transpose(np.loadtxt("Datos_presion.txt", skiprows=1))

datos_promedio = datos.reshape(3, bloques , tamaño_bloque).mean(axis=2)

velocidades = datos_promedio[0]
temperaturas = datos_promedio[1]
presiones = datos_promedio[2]

x = temperaturas
y = presiones

plt.plot(x, y, 'b-')

dx = x.max() - x.min()
dy = y.max() - y.min()
plt.xlim([x.min() - 0.1 * dx, x.max() + 0.1 * dx])
plt.ylim([y.min() - 0.1 * dy, y.max() + 0.1 * dy])
plt.xlabel("Temperatura del sistema")
plt.ylabel("Presión del sistema")

plt.show()