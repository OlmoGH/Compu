import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-0.99, 0.99, 1000)
y = 1 / (np.pi * np.sqrt(1 - x ** 2))

data = np.random.uniform(0, 2 * np.pi, 1000000)
data = np.cos(data)

plt.hist(data, 50, density=True, alpha=0.8, edgecolor='k')
plt.plot(x, y, 'r-')
plt.grid(True)
plt.title(r"Distribución de probabilidad de $v_x$ y $v_y$ para $|v| = 1$")

plt.ylim([0, 2])
plt.xlabel(r"$v_i$")
plt.ylabel("PDF")
plt.savefig("../Animaciones y graficas/PDF arcoseno.png")