import numpy as np
import matplotlib.pyplot as plt

# === Parámetros de simulación ===
N, L, pasos, mod_v = np.loadtxt("../Simulacion/Datos/Parametros.txt", usecols=1, dtype=int, converters=float)

archivo_vel = f"../Simulacion/Datos/Velocidades/Velocidades_{N}_{pasos}_{L}_{mod_v}.csv"

# === Cargar datos ===
datos_vel = np.genfromtxt(archivo_vel, delimiter=',', skip_header=1)

print(f"Número de partículas: {N}")

# === Análisis ===
# Se descartan primeros 'inicio' pasos
inicio = 10000
print(f"tamaño = {np.size(datos_vel)}")
if datos_vel.size == 2 * N:
    vx = datos_vel[0::2]
    vy = datos_vel[1::2]
else:
    vx = datos_vel[inicio:, 0::2]
    vy = datos_vel[inicio:, 1::2]

v_cuadrado = vx**2 + vy**2
velocidades = np.sqrt(v_cuadrado)
v_flat = velocidades.ravel()
vx_flat = vx.ravel()
d = 2  # sistema 2D
T = np.mean(v_cuadrado) / d

print(f"|v| = {mod_v}")
print(f"T = {T:.3f}")

# === Distribuciones teóricas ===
v_lin = np.linspace(0, v_flat.max(), 500)
vx_lin = np.linspace(vx_flat.min(), vx_flat.max(), 500)

dist_v = v_lin / T * np.exp(-v_lin ** 2 / (2 * T))
dist_vx = 1 / np.sqrt(2 * np.pi * T) * np.exp(-(vx_lin) ** 2 / (2 * T))

# === Gráfica ===
fig, (ax_v, ax_vx) = plt.subplots(ncols=2, figsize=(10, 5))
ax_v.grid(True)
ax_vx.grid(True)

ax_v.hist(v_flat, bins=100, color='skyblue', edgecolor='black', density=True, label='Histograma v')
ax_vx.hist(vx_flat, bins=100, color='teal', edgecolor='black', density=True, label=r'Histograma $v_x$')


ax_v.plot(v_lin, dist_v, 'r-', label='Maxwell-Boltzmann')
ax_vx.plot(vx_lin, dist_vx, 'r-', label='Gauss')

ax_v.set_xlabel("Módulo de la velocidad")
ax_vx.set_xlabel("Velocidad horizontal")
ax_v.set_ylabel("Densidad de probabilidad")
ax_vx.set_ylabel("Densidad de probabilidad")
fig.suptitle(fr"Distribución de velocidades para {N} partículas y t $\in$ [20, {int(pasos * 0.002)}]")
ax_v.legend()
ax_vx.legend()
plt.show()
# fig.savefig(f"../Animaciones y graficas/Distribuciones de probabilidad {pasos} pasos y {N} partículas.png")
