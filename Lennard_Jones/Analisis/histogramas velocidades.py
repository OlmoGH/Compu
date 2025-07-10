import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

viridis = mpl.colormaps['viridis'].resampled(4)

# === Parámetros de simulación ===
N, L, pasos, mod_v = np.loadtxt("../Simulacion/Datos/Parametros.txt", usecols=1, dtype=int, converters=float)

inicio = 10000
datos = np.empty((4, pasos - inicio, 2 * N))
for i, v in enumerate([1, 4, 7, 10]):
    nombre = f"../Simulacion/Datos/Velocidades/Velocidades_{N}_{pasos}_{L}_{v}.csv"
    dato_vel = np.genfromtxt(nombre, delimiter=',', skip_header=1)
    datos[i] = dato_vel[inicio:]

vx = datos[:, :, 0::2]
vy = datos[:, :, 1::2]

v_cuadrado = vx**2 + vy**2
velocidades = np.sqrt(v_cuadrado)
v_flat = velocidades.reshape(4, -1)
vx_flat = vx.reshape(4, -1)
d = 2  # sistema 2D
T = np.mean(v_cuadrado, axis=(1, 2)) / d

print(f"T = {T}")

# === Gráfica ===
fig, (ax_v, ax_vx) = plt.subplots(ncols=2, figsize=(10, 5))
ax_v.grid(True)
ax_vx.grid(True)

for i in range(4):
    v_lin = np.linspace(0, 20, 500)
    vx_lin = np.linspace(-20, 20, 500)
    
    ax_v.plot(v_lin, v_lin / T[i] * np.exp(-v_lin ** 2 / (2 * T[i])), color=viridis(i), linewidth=2, zorder=i+2)
    ax_vx.plot(vx_lin, 1 / np.sqrt(2 * np.pi * T[i]) * np.exp(-vx_lin ** 2 / (2 * T[i])), color=viridis(i), linewidth=2, zorder=i+2)

    ax_v.hist(v_flat[i], bins=30, density=True, color=viridis(i), label=f"T = {T[i]:.2f}", rwidth=0.8, alpha=0.8, zorder=i+2)
    ax_vx.hist(vx_flat[i], bins=30, density=True, color=viridis(i), label=f'T = {T[i]:.2f}', rwidth=0.8, alpha=0.8, zorder=i+2)

ax_v.set_xlabel("Módulo de la velocidad")
ax_vx.set_xlabel("Velocidad horizontal")
ax_v.set_ylabel("Densidad de probabilidad")
ax_vx.set_ylabel("Densidad de probabilidad")
ax_v.set_xlim([0, 10])
ax_vx.set_xlim([-10, 10])
fig.suptitle(
    fr"Distribuciones de velocidades para {N} partículas" + "\n" + 
    fr"t $\in$ [20, {int(pasos * 0.002)}] $\quad$ $v_0 =  \{{1, 4, 7, 10 \}}$"
)
ax_v.legend()
ax_vx.legend()
plt.show()
fig.savefig(f"../Animaciones y graficas/Histogramas de 4 velocidades y {N} partículas.png")
