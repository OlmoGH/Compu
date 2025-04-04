import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

cuerpos = 6
datos = {}

# Leer los archivos de datos
for i in range(cuerpos):
    file_name = f"cuerpo_{i}.txt"
    datos[i] = pd.read_csv(file_name, delimiter='\t')

# Configurar la figura
plt.figure(figsize=(8, 8))
for i in range(1, cuerpos):  # Omitimos el Sol si está en el origen
    plt.plot(datos[i]["POSICION_X"], datos[i]["POSICION_Y"], label=f'Cuerpo {i}')

# Añadimos el Sol en el origen
plt.plot(0, 0, 'yo', label='Sol')

# Mejoras visuales
plt.xlabel("x (UA)")
plt.ylabel("y (UA)")
plt.title("Órbitas de los planetas")
plt.grid(True)
plt.gca().set_aspect('equal')  # Para que no se distorsione la órbita
plt.legend()

# Ajuste de los límites (automático o manual)
#plt.xlim([-2, 2])
#plt.ylim([-2, 2])

plt.show()

