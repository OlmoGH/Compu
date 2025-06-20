import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# === CONFIGURACIÓN DEL LOS ARCHIVOS DE LECTURA ===
parametros_simulacion = np.loadtxt("Parametros.txt", usecols=1, dtype=int)
N = parametros_simulacion[0]
L = parametros_simulacion[1]
pasos = parametros_simulacion[2]

inicio = 0

archivo_pos = f"Datos/Lennard_Jones_{N}_{pasos}_{L}.csv"
archivo_ene = f"Datos/Energias_{N}_{pasos}_{L}.csv"

# === CARGAR DATOS ===

datos_pos = np.genfromtxt(archivo_pos, delimiter=',', skip_header=1)  # asume que hay cabecera
datos_ene = np.genfromtxt(archivo_ene, delimiter=',', skip_header=1)

print(f"El número de partículas es {N}")
print(f"El valor más alto de la energia es {np.max(datos_ene)}")

num_pasos = datos_pos.shape[0]

# === CONFIGURAR FIGURA Y SUBPLOTS ===
fig, (ax1, ax_energia) = plt.subplots(1, 2, figsize=(10, 5))

# Subplot de partículas
particles, = ax1.plot([], [], 'bo', ms=6)
ax1.set_xlim(0, L)
ax1.set_ylim(0, L)
ax1.set_aspect('equal')
ax1.grid(True)

# Subplot del histograma
ax_energia.set_xlim(0, num_pasos)
ax_energia.set_ylim(np.min(datos_ene) * 1.5, np.max(datos_ene) * 1.5)  # Ajusta según se vea mejor
ax_energia.set_title("Energías del sistema")
ax_energia.set_xlabel("Paso")
ax_energia.set_ylabel("Energía")
cinetica, = ax_energia.plot([], [], 'r-', label="Energía cinética")
potencial, = ax_energia.plot([], [], 'b-', label="Energía potencial")
total, = ax_energia.plot([], [], 'g-', label="Energía total")
ax_energia.legend()

# === FUNCIONES DE ANIMACIÓN ===

def init():
    particles.set_data([], [])
    cinetica.set_data([], [])
    potencial.set_data([], [])
    total.set_data([], [])
    return [particles, cinetica, potencial, total]

def actualizar(frame):
    fila = datos_pos[frame + inicio]
    x = fila[::2]
    y = fila[1::2]
    particles.set_data(x, y)
    print(f"Animando frame {frame}", end="\r")


    cinetica.set_data(np.arange(frame + 1), datos_ene[inicio:frame + 1 + inicio, 1])
    potencial.set_data(np.arange(frame + 1), datos_ene[inicio:frame + 1 + inicio, 0])
    total.set_data(np.arange(frame + 1), np.sum(datos_ene, axis=1)[inicio:frame + 1 + inicio])

    return [particles, cinetica, potencial, total]


# ani = animation.FuncAnimation(
#     fig, actualizar, frames=num_pasos - inicio,
#     init_func=init, interval=10, blit=True 
#     )

plt.tight_layout()
#ani.save(f"Animacion_{N}_{pasos}_{L}.mp4", dpi=150, writer="ffmpeg")
# plt.show()
