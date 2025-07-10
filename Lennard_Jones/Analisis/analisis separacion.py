import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import pandas as pd

parametros_simulacion = np.loadtxt("C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Parametros.txt", usecols=1)
N = int(parametros_simulacion[0])
Lx = int(parametros_simulacion[1])
Ly = int(parametros_simulacion[2])
pasos = int(parametros_simulacion[3])
v = int(parametros_simulacion[4])

nombre = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Variacion_posicion/Separacion_{N}_{pasos}_{Lx}_{v}.csv"
nombre_velocidades = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Velocidades/Velocidades_{N}_{pasos}_{Lx}_{v}_0.csv"
datos = pd.read_csv(nombre, delimiter=',', header=0).to_numpy()
datos_vel = pd.read_csv(nombre_velocidades, delimiter=',', header=0).to_numpy()

velocidades_utiles = datos_vel[600000:900000]
vx = velocidades_utiles[0::2]
vy = velocidades_utiles[1::2]
v_cuadrado = vx ** 2 + vy ** 2
T = np.mean(v_cuadrado) / 2
N_datos = np.size(v_cuadrado.ravel())
dT = np.std(v_cuadrado) / (2 * np.sqrt(N_datos))
print(T)
print(dT)


incremento = 300000
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(datos, linewidth=1)
ax.set_xlim([0, len(datos)])
y_max = ax.get_ylim()[1]
for i, position in enumerate(range(incremento, 10 * incremento, incremento)):
    ax.axvline(position, color='r', linestyle='-')

fluido = Rectangle([2100000, 0], width=900000, height=y_max, color='r', alpha=0.4, label="Estado fluido")
solido = Rectangle([0, 0], width=2100000, height=y_max, color='b', alpha=0.4, label="Estado sólido")
ax.add_patch(solido)
ax.add_patch(fluido)
ax.set_title("Separación cuadrática media entre dos partículas con el tiempo \n" + r" $\langle (r_i(t)-r_j(t))^2 \rangle$ ")
ax.legend()

#plt.savefig("../Animaciones y graficas/Separacion.png")