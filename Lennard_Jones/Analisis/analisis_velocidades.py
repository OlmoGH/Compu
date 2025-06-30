import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import maxwell

# === Parámetros de simulación ===
N, L, pasos, mod_v = np.loadtxt("Parametros.txt", usecols=1, dtype=int)

archivo_pos = f"Datos/Lennard_Jones_{N}_{pasos}_{L}.csv"
archivo_vel = f"Datos/Velocidades_{N}_{pasos}_{L}.csv"

# === Cargar datos ===
datos_pos = np.genfromtxt(archivo_pos, delimiter=',', skip_header=1)
datos_vel = np.genfromtxt(archivo_vel, delimiter=',', skip_header=1)

print(f"Número de partículas: {N}")

# === Análisis ===
# Se descartan primeros 'inicio' pasos
inicio = 10000
fin = pasos
vx = datos_vel[inicio:fin, 0::2]
vy = datos_vel[inicio:fin, 1::2]

v_cuadrado = vx**2 + vy**2
velocidades = np.sqrt(v_cuadrado)
v_flat = velocidades.ravel()
vx_flat = vx.ravel()
d = 2  # sistema 2D
T = np.mean(v_cuadrado) / d

print(f"|v| = {mod_v}")
print(f"T = {T:.3f}")

# === Gráfica ===
fig, (ax_v, ax_vx) = plt.subplots(ncols=2, figsize=(10, 5))

ax_v.hist(v_flat, bins=100, color='skyblue', edgecolor='black', density=True, label='Histograma v')
ax_vx.hist(vx_flat, bins=100, color='teal', edgecolor='black', density=True, label=r'Histograma $v_x$')

v_lin = np.linspace(0, v_flat.max(), 500)
vx_lin = np.linspace(vx_flat.min(), vx_flat.max(), 500)
ax_v.plot(v_lin, v_lin / T * np.exp(-v_lin ** 2 / (2 * T)), 'r-', label='Maxwell-Boltzmann')
ax_vx.plot(vx_lin, 1 / np.sqrt(2 * np.pi * T) * np.exp(-vx_lin ** 2 / (2 * T)), 'r-', label='Gauss')

ax_v.set_xlabel("Módulo de la velocidad")
ax_vx.set_xlabel("Velocidad horizontal")
ax_v.set_ylabel("Densidad de probabilidad")
ax_vx.set_ylabel("Densidad de probabilidad")
fig.suptitle("Distribución de velocidades")
ax_v.legend()
ax_vx.legend()
plt.grid(True)
fig.savefig("Distribuciones de probabilidad.pdf")
