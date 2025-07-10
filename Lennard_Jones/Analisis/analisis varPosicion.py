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

nombre = f"C:/Users/olmov/Desktop/compu/Compu/Lennard_Jones/Simulacion/Datos/Variacion_posicion/Var_posicion_{N}_{pasos}_{Lx}_{v}.csv"
datos = pd.read_csv(nombre, delimiter=',', header=0).to_numpy()

incrementos = [200000, 300000, 350000, 450000]
fig, ax = plt.subplots(figsize=(8, 5))
ax.set_xlim([0, len(datos)])
ax.plot(datos, linewidth=1)
for position in incrementos:
    ax.axvline(position, color='r', linestyle='-')

y_max = ax.get_ylim()[1]
y_min = ax.get_ylim()[0]
print(y_max)

fluido = Rectangle([incrementos[2], y_min], width=500000-incrementos[2], height=y_max - y_min, color='r', alpha=0.4, label="Fase fluida")
solido = Rectangle([0, y_min], width=incrementos[2], height=y_max - y_min, color='b', alpha=0.4, label="Fase sólida")
ax.add_patch(solido)
ax.add_patch(fluido)
ax.legend()
ax.set_title("Variación cuadrática media de la posición con el tiempo \n" + r" $\langle (r(t)-r_0)^2 \rangle$ ")
plt.savefig("../Animaciones y graficas/VarPosicion")