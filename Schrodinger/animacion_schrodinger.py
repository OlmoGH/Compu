import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks
import random

skip = 1

N, n_ciclos, lamda, L, sigma, mu, s, pasos, inicio, ancho = np.loadtxt("datos.txt")
potencial = np.loadtxt("V.txt")

datos = np.loadtxt("Estados.txt")

datos_norma = np.loadtxt("Norma.txt")
peaks, _ = find_peaks(datos_norma)
print(f"n Max = {peaks}")

max_norma = datos_norma[peaks[0]]

N_medidas = 1000000
m = 0
for i in range(N_medidas):
    p = random.random()
    if p > max_norma:
        m += 1

print(f"K = {m / N_medidas}")

plt.plot(peaks, datos_norma[peaks], 'o')
plt.plot(datos_norma)
plt.show()

fig, ax = plt.subplots()

funcion, = ax.plot([], [], 'r-')
ax.plot(potencial[::skip], 'k-')
rango = np.arange(N // skip)

y_max = datos.max() / 3
ax.set_xlim([0, N / skip])
ax.set_ylim([0, y_max])

detector = Rectangle([inicio, 0], width=ancho, height=y_max, alpha=0.4, color='purple')
ax.add_patch(detector)

def update(frame):
    funcion.set_data(rango, datos[frame, ::skip])
    return funcion, 

animacion = FuncAnimation(fig, update, frames=int(pasos), blit=True, interval=10)
plt.tight_layout()
plt.show()
#animacion.save("Animacion_schrodinger.mp4")
