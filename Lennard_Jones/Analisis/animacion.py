import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import pandas as pd

# === CONFIGURACIÓN DEL LOS ARCHIVOS DE LECTURA ===
parametros_simulacion = np.loadtxt("C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Parametros.txt", usecols=1)
N = int(parametros_simulacion[0])
Lx = int(parametros_simulacion[1])
Ly = int(parametros_simulacion[2])
pasos = int(parametros_simulacion[3])
v_inicial = int(parametros_simulacion[4])

skip = 5000
inicio = 0

archivo_pos = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Posiciones/Lennard_Jones_{N}_{pasos}_{Lx}_{v_inicial}.csv"
archivo_ene = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Energias/Energias_{N}_{pasos}_{Lx}_{v_inicial}.csv"

# === CARGAR DATOS ===

datos_pos = pd.read_csv(archivo_pos, delimiter=',', header=0).to_numpy()  # asume que hay cabecera
datos_ene = pd.read_csv(archivo_ene, delimiter=',', header=0).to_numpy()
plt.plot(datos_pos[0, ::2], datos_pos[0, 1::2], 'o')
plt.ylim([0, Ly])
plt.xlim([0, Lx])
plt.show()
plt.plot(datos_ene)
plt.plot(datos_ene[:, 0]+datos_ene[:, 1])
plt.show()
print(f"El número de partículas es {N}")
print(f"El valor más alto de la energia es {np.max(datos_ene)}")


# Las filas datos útiles vienen desde 'inicio' en pasos 'skip'
filas = datos_pos[inicio::skip]
energias = datos_ene[inicio::skip]
num_pasos = len(filas)
frames = range(num_pasos)

datos_cinetica = energias[:, 1]
datos_potencial = energias[:, 0]
datos_total = energias.sum(axis=1)


# === CONFIGURAR FIGURA Y SUBPLOTS ===
fig, (ax1, ax_energia) = plt.subplots(1, 2, figsize=(10, 5))

# Subplot de partículas
particles, = ax1.plot([], [], 'bo', ms=6)
ax1.set_xlim(0, Lx)
ax1.set_ylim(0, Ly)
ax1.set_aspect('equal')
ax1.grid(True)

# Subplot del histograma
ax_energia.set_xlim(0, num_pasos)
ax_energia.set_ylim(np.min(datos_ene) * 1.5, np.max(datos_ene) * 1.5)  # Ajusta según se vea mejor
ax_energia.set_xlabel("Paso")
ax_energia.set_ylabel("Energía")
fig.suptitle(f"Simulación de {N} partículas con paso 0.0001 y velocidad inicial {v_inicial}")

cinetica, = ax_energia.plot([], [], 'r-', label="Energía cinética")
potencial, = ax_energia.plot([], [], 'b-', label="Energía potencial")
total, = ax_energia.plot([], [], 'g-', label="Energía total")
prev_cinetica, = ax_energia.plot([], [], 'r-', alpha=0.5)
prev_potencial, = ax_energia.plot([], [], 'b-', alpha=0.5)
prev_total, = ax_energia.plot([], [], 'g-', alpha=0.5)

ax_energia.legend()

# === FUNCIONES DE ANIMACIÓN ===

def init():
    particles.set_data([], [])
    cinetica.set_data([], [])
    potencial.set_data([], [])
    total.set_data([], [])
    prev_cinetica.set_data([], [])
    prev_potencial.set_data([], [])
    prev_total.set_data([], [])
    return [particles, cinetica, potencial, total, prev_cinetica, prev_potencial, prev_total]

def actualizar(frame):
    fila_actual = filas[frame]
    x = fila_actual[::2]
    y = fila_actual[1::2]
    particles.set_data(x, y)

    cinetica_values = datos_cinetica[:frame]
    potencial_values = datos_potencial[:frame]
    total_values = datos_total[:frame]

    tam = np.size(cinetica_values)

    prev_cinetica.set_data(np.arange(num_pasos), datos_cinetica)
    prev_potencial.set_data(np.arange(num_pasos), datos_potencial)
    prev_total.set_data(np.arange(num_pasos), datos_total)
    cinetica.set_data(np.arange(tam), cinetica_values)
    potencial.set_data(np.arange(tam), potencial_values)
    total.set_data(np.arange(tam), total_values)

    print(f"Animando frame {frame} de {num_pasos}")

    return [particles, cinetica, potencial, total, prev_cinetica, prev_potencial, prev_total]


ani = animation.FuncAnimation(
     fig, actualizar, frames=np.size(frames),
     init_func=init, interval=20, blit=True 
     )

plt.show()
# ani.save(f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Animaciones y graficas/Animacion primer calentamiento.mp4", writer='ffmpeg')
