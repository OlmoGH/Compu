import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker

viridis = mpl.colormaps['viridis'].resampled(5)

# === Parámetros de simulación ===
N, L, pasos, mod_v = np.loadtxt("../Simulacion/Datos/Parametros.txt", usecols=1, dtype=int, converters=float)

rangos = np.array([0, 1, 2, 3, 4])
datos = np.empty((5, pasos, 2 * N))

for i in rangos:
    archivo_vel = f"../Simulacion/Datos/Velocidades/Velocidades_{N}_{pasos}_{L}_{mod_v}_{i}.csv"
    curr_data = np.genfromtxt(archivo_vel, delimiter=',', skip_header=1)
    datos[i] = curr_data

print(datos)
# === Análisis ===
# Se descartan primeros 'inicio' pasos
inicio = 10000
vx = datos[:, inicio:, 0::2].reshape(5, -1)
vy = datos[:, inicio:, 1::2].reshape(5, -1)

v_cuadrado = (vx * vx + vy * vy).mean(axis=1)

media = rangos
var =  np.full(5, v_cuadrado[0] / 2)

# === Distribuciones teóricas ===
vx_lin = np.linspace(vx.min(), vx.max(), 500)

dist_vx = np.array([1 / np.sqrt(2 * np.pi * var[i]) * np.exp(-(vx_lin - media[i]) ** 2 / (2 * var[i])) for i in range(5)])

# === Gráfica ===
fig, ax = plt.subplots()

ax.grid(True)
for i in range(5): 
    ax.plot(vx_lin, dist_vx[i], color=viridis(i), zorder=i+4, linewidth=2)
    ax.hist(vx[i], bins=30, color = viridis(i), density=True, alpha=0.8, zorder=i+3, rwidth=0.8, label=fr"$v_{{xinicial}} = {i}$")

ax.set_xlabel("Velocidad horizontal")
ax.set_ylabel("Densidad de probabilidad")
ax.set_xlim([-3, 7])
fig.suptitle(fr"Distribuciones de velocidades en el eje x para {N} partículas " + "\n" + rf" t $\in$ [20, {int(pasos * 0.002)}] y diferentes impulsos iniciales")
ax.legend()

ax.xaxis.set_major_locator(ticker.MultipleLocator(1))
plt.show()
fig.savefig(f"../Animaciones y graficas/Distribuciones eje x {N} partículas diferentes impulsos.png")
