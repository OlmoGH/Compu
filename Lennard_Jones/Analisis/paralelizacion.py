import matplotlib.pyplot as plt
import numpy as np

num_hilos = 6
datos = np.genfromtxt("Datos_paralelizacion.csv", delimiter=",")
fig, (ax, ax2) = plt.subplots(1, 2, figsize=(10, 5))
for i in range(datos.shape[0] // num_hilos):
    ax.plot(datos[num_hilos * i: num_hilos * i + num_hilos, 1], datos[num_hilos * i: num_hilos * i + num_hilos, 2])

ax.set_xlabel("Hilos")
ax.set_ylabel("Tiempo")
ax.set_title("Dependencia con el número de hilos \n para diferentes números de partículas")

for i in range(num_hilos):
    ax2.plot(datos[i::num_hilos, 0], datos[i::num_hilos, 2])

ax2.set_xlabel("Particulas")
ax2.set_ylabel("Tiempo")
ax2.set_title("Dependencia con el número de partículas \n para diferentes hilos")
plt.tight_layout()
plt.show()