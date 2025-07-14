import numpy as np
import matplotlib.pyplot as plt

datos = np.loadtxt("Paralelizacion.txt", skiprows=1).transpose()

N_values = np.unique(datos[0])

N_threads = [1, 2, 3, 4, 5, 6]

tiempos = np.reshape(datos[2], (N_values.size, -1))

for i, N in enumerate(N_values):
    plt.plot(N_threads, tiempos[i])

plt.show()