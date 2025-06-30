import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
import pandas as pd

h = 60
intervalo = 10
datos = pd.read_csv("Posiciones.txt", delimiter='\t', header=0, names=['r', 'phi'])
energia = np.loadtxt("Energia.txt", dtype=float, delimiter='\t', skiprows=1)

fig, (ax, ax2) = plt.subplots(ncols = 2, figsize=[15, 5])
x = datos.r * np.cos(datos.phi)
y = datos.r * np.sin(datos.phi)
cohete, = ax.plot([], [], 'b-')
tierra = Circle(xy=[0.0, 0.0], radius=0.0165925078, color='blue', fill=False)
luna = Circle(xy=[1.0, 0.0], radius=0.0165925078, color='red', fill=False)
w = 0.0000026617
ax.add_patch(tierra)
ax.add_patch(luna)

T, = ax2.plot([], [], label="Energia cinetica")
V, = ax2.plot([], [], label="Energia potencial")
H, = ax2.plot([], [], label="Energia mecanica")

def update(frame):
    frame = frame * intervalo
    luna.set_center([np.cos(w * frame * h), np.sin(w * frame * h)])
    cohete.set_data(x[:frame], y[:frame])
    T.set_data(np.arange(frame), energia[:frame, 0])
    V.set_data(np.arange(frame), energia[:frame, 1])
    H.set_data(np.arange(frame), energia[:frame, 2])
    return luna, cohete, T, V, H



ax.set_xlim([-1.5, 1.5])
ax.set_ylim([-1.5, 1.5])
ax2.set_xlim([0, len(x)])
ax2.set_ylim([np.min(energia), np.max(energia)])
ax.set_aspect('equal')
plt.legend()
animacion = FuncAnimation(fig, func=update, frames=len(x) // intervalo, blit=True, interval=10)
plt.show()