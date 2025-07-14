import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Rectangle
from scipy.signal import find_peaks
import random

skip = 5

N, n_ciclos, lamda, L, sigma, mu, s, pasos, inicio, ancho = np.loadtxt("datos.txt")
potencial = np.loadtxt("V.txt")

datos = np.loadtxt("Estados.txt")

datos_norma = np.loadtxt("Norma.txt")

fig, (ax, ax_norma) = plt.subplots(nrows=2, gridspec_kw={'height_ratios':[5, 1]})

funcion, = ax.plot([], [], 'r-')
norma, = ax_norma.plot([], [], 'b-')

rango = np.linspace(0, 1, int(N))
ax.plot(rango, potencial * 10, 'k-')
ax_norma.set_xlim([0, pasos / skip])
ax_norma.set_ylim([0, 2])
ax_norma.set_xlabel("Tiempo")
ax_norma.set_ylabel("Norma")
ax.set_ylim([0, 10])
ax.set_xlim([0, 1])
ax.set_xticks([0, 1], labels=['x = 0', 'x = L'])


def update(frame):
    funcion.set_data(rango, datos[frame * skip])
    norma.set_data(np.arange(frame), datos_norma[:frame * skip:skip])
    print(f"Animando frame {frame} de {pasos/skip}")
    return funcion, norma

animacion = FuncAnimation(fig, update, frames=int(pasos / skip), blit=True, interval=20)
plt.tight_layout()
animacion.save("Animacion_con_norma.gif")
