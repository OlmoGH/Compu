import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd

# Número de cuerpos (incluyendo el Sol)
cuerpos = 10
factor = 100
posiciones = {}
momentos = {}
energias = {}

# Leer los archivos de datos
for i in range(cuerpos):
    name_pos = f"Datos_cuerpos/cuerpo_{i}.txt"
    posiciones[i] = pd.read_csv(name_pos, delimiter='\t', header=0, names=['x', 'y'])


energia_total = pd.read_csv("Datos_sistema/energia_total.txt", delimiter='\t', header=0).to_numpy()
momento_total = pd.read_csv("Datos_sistema/momento_total.txt", delimiter='\t', header=0).to_numpy()


# Nombres de los cuerpos
nombres = ['Sol', 'Mercurio', 'Venus', 'Tierra', 'Marte',
            'Júpiter', 'Saturno', 'Urano', 'Neptuno', 'Plutón']

# Crear la figura
fig, ax = plt.subplots(figsize=(8, 8))

# Crear puntos (círculos) para cada planeta
trayectorias = []
planetas = []
for i in range(cuerpos):
    (planeta,) = ax.plot([], [], 'o', label=nombres[i])
    planetas.append(planeta)
    (trayectoria,) = ax.plot([], [], '-')
    trayectorias.append(trayectoria)


# Función de inicialización
def init():
    for planeta, trayectoria in zip(planetas, trayectorias):
        planeta.set_data([], [])
        trayectoria.set_data([], [])
    return trayectorias + planetas

# Función de actualización
def update(frame):
    frame = factor * frame
    for i in range(cuerpos):
        x = posiciones[i].x[frame]
        y = posiciones[i].y[frame]
        planetas[i].set_data([x], [y])  # Usamos listas para que funcione
        trayectorias[i].set_data(posiciones[i].x[:frame], posiciones[i].y[:frame])
    return trayectorias + planetas

# Crear animación
frames = len(posiciones[0])//factor
anim = FuncAnimation(fig, update, init_func=init, frames=frames, blit=True, interval=10)


fig2, ax2 = plt.subplots(ncols=2)
energia, = ax2[0].plot([], [], 'r-')
momento, = ax2[1].plot([], [], 'b-')

def update2(frame):
    frame = factor * frame
    energia.set_data(np.arange(frame), energia_total[:frame])
    momento.set_data(np.arange(frame), momento_total[:frame])
    return energia, momento

anim2 = FuncAnimation(fig2, update2,frames=frames, blit=True, interval=10)


# Configurar límites y proporciones
ax2[0].set_xlim(0, len(energia_total))
ax2[1].set_xlim(0, len(momento_total))
ax2[0].set_ylim(min(energia_total[1:]), max(energia_total[1:]))
ax2[1].set_ylim(min(momento_total)-0.01, max(momento_total)+0.01)
ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)
ax.set_aspect('equal')
ax.set_xlabel("x (UA)")
ax.set_ylabel("y (UA)")
ax.set_title("Órbitas de los planetas")
ax.grid(True)
ax.legend(loc='upper right')

# Mostrar animación
plt.show()

