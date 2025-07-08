import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

# === Parámetros de simulación ===
factor = 10
N, L, pasos, mod_v = np.loadtxt("../Simulacion/Datos/Parametros.txt", usecols=1, converters=float, dtype=int)

archivo_vel = f"../Simulacion/Datos/Velocidades/Velocidades_{N}_{pasos}_{L}.csv"
archivo_pos = f"../Simulacion/Datos/Posiciones/Lennard_Jones_{N}_{pasos}_{L}.csv"

# === Cargar datos ===
datos_vel = np.genfromtxt(archivo_vel, delimiter=',', skip_header=1)
datos_pos = np.genfromtxt(archivo_pos, delimiter=',', skip_header=1)

print(f"Número de partículas: {N}")

# === Análisis ===
## Velocidad
vx = datos_vel[:, 0::2]
vy = datos_vel[:, 1::2]
v_cuadrado = vx**2 + vy**2
velocidades = np.sqrt(v_cuadrado)

## Posición
x = datos_pos[:, 0::2]
y = datos_pos[:, 1::2]

# ================
d = 2  # sistema 2D
T = np.mean(v_cuadrado) / d

print(f"|v| = {mod_v}")
print(f"T = {T:.3f}")

# === Animación ===
fig, (ax_v, ax_anim) = plt.subplots(ncols=2, figsize=(10, 5))

# Calcular curva teórica primero
v_lin = np.linspace(0, velocidades.max(), 500)
dist_teorica = v_lin / T * np.exp(-v_lin**2 / (2 * T))
max_teorico = np.max(dist_teorica)

# Precalcular todos los histogramas para obtener el máximo global
num_bins = 30


# Calcular límites de bins usando todos los datos
_, bins_edges = np.histogram(velocidades.ravel(), bins=num_bins)


# Inicializar histograma vacío
n, bins_edges, barContainer = ax_v.hist([], bins=bins_edges, color='skyblue', 
                                        edgecolor='black', density=True, 
                                        label='Histograma v',
                                        zorder=1)
histograma = list(barContainer)
particulas = plt.scatter([], [], c='b')

# Configurar ejes
ax_v.grid(True)
ax_anim.set_xlim([0, L])
ax_anim.set_ylim([0, L])
ax_v.set_xlim(0, bins_edges.max())
ax_v.set_ylim(0, 0.125)  # FIJAR LÍMITE DEL EJE Y

# Graficar curva teórica (después de establecer límites)
ax_v.plot(v_lin, dist_teorica, 'r-', label='Maxwell-Boltzmann')

ax_v.set_xlabel("Módulo de la velocidad")
ax_v.set_ylabel("Densidad de probabilidad")
fig.suptitle(f"Histograma y distrubución de velocidades junto a la simulación")
ax_v.legend()

def init():
    for bar in histograma:
        bar.set_height(0)

    particulas.set_offsets(np.empty((0, 2)))

    return [particulas] + histograma

def update(frame):
    # Usar histogramas precalculados
    alturas, _ = np.histogram(velocidades[frame:frame+30, :].ravel(), bins=bins_edges, density=True)
    for i, bar in enumerate(histograma):
        bar.set_height(alturas[i])

    posiciones = np.column_stack((x[frame], y[frame]))
    particulas.set_offsets(posiciones)
    print(f"Animando frame {frame} de {pasos}")
    return [particulas] + histograma

animacion = FuncAnimation(fig=fig, func=update, frames=range(pasos), 
                         init_func=init, blit=True, interval=20)
plt.show()
animacion.save(f"../Animaciones y graficas/Animacion velocidades {N} partículas y {pasos} pasos.mp4", writer='ffmpeg')
