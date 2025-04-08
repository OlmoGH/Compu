import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import pandas as pd
import os

# Número de cuerpos (incluyendo el Sol)
cuerpos = 10
datos = {}

# Leer los archivos de datos
for i in range(cuerpos):
    file_name = f"Datos_cuerpos/cuerpo_{i}.txt"
    if os.path.exists(file_name):
        datos[i] = pd.read_csv(file_name, delimiter='\t', header=0, names=['x', 'y'])
    else:
        print(f"Archivo no encontrado: {file_name}")
        datos[i] = pd.DataFrame({'x': [0], 'y': [0]})  # Valor de emergencia

# Nombres de los cuerpos
nombres = ['Sol', 'Mercurio', 'Venus', 'Tierra', 'Marte',
            'Júpiter', 'Saturno', 'Urano', 'Neptuno', 'Plutón']

# Crear la figura
fig, ax = plt.subplots(figsize=(8, 8))

# Configurar límites y proporciones
ax.set_xlim(-50, 50)
ax.set_ylim(-50, 50)
ax.set_aspect('equal')
ax.set_xlabel("x (UA)")
ax.set_ylabel("y (UA)")
ax.set_title("Órbitas de los planetas")
ax.grid(True)

# Crear puntos (círculos) para cada planeta
trayectorias = []
planetas = []
for i in range(cuerpos):
    (planeta,) = ax.plot([], [], 'o', label=nombres[i])
    planetas.append(planeta)
    (trayectoria,) = ax.plot([], [], '-')
    trayectorias.append(trayectoria)


# Añadir leyenda
ax.legend(loc='upper right')

# Función de inicialización
def init():
    for planeta, trayectoria in zip(planetas, trayectorias):
        planeta.set_data([], [])
        trayectoria.set_data([], [])
    return trayectorias + planetas

# Función de actualización
def update(frame):
    for i in range(cuerpos):
        x = datos[i].x[frame]
        y = datos[i].y[frame]
        planetas[i].set_data([x], [y])  # Usamos listas para que funcione
        trayectorias[i].set_data(datos[i].x[:frame], datos[i].y[:frame])
    return planetas + trayectorias

# Crear animación
frames = len(datos[0])
anim = FuncAnimation(fig, update, init_func=init, frames=frames, blit=True, interval=0)

# Mostrar animación
plt.show()
