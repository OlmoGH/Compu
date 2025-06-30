import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

variables = np.loadtxt("datos.txt")
norma = np.loadtxt("Norma.txt")
potencial = np.loadtxt("V.txt")

datos = np.loadtxt("Estados.txt")

fig, (ax, ax_norma) = plt.subplots(ncols=2, figsize=(10, 5))
funcion, = ax.plot([], [], 'r-')
plot_norma, = ax_norma.plot([], [], 'b-')
ax.plot(potencial, 'k-')
rango = np.arange(variables[0])

ax.set_xlim(0, int(variables[0]))
ax.set_ylim(0, 1.5)
ax_norma.set_xlim(0, int(variables[7]))
ax_norma.set_ylim(np.mean(norma) * 0.8, np.mean(norma) * 1.2)

def update(frame):
    funcion.set_data(rango, datos[frame])
    plot_norma.set_data(np.arange(frame), norma[:frame])
    return funcion, plot_norma

animacion = FuncAnimation(fig, update, frames=int(variables[7]), blit=True, interval=10)
plt.tight_layout()
plt.show()
#animacion.save("Animacion_schrodinger.mp4")
