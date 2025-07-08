import numpy as np
import matplotlib.pyplot as plt

# Definir rango alrededor del mínimo (1.122σ)
x = np.linspace(0.9691155, 3.0, 1000)  # Centrado en ~1.122σ

# Potencial de Lennard-Jones (con σ=1)
y = 4 * (np.power(x, -12) - np.power(x, -6))

# Graficar
plt.figure(figsize=(10, 10))
plt.plot(x, y, 'b-', linewidth=2)
plt.title('Potencial de Lennard-Jones', fontsize=25)
plt.xlabel('r/σ', fontsize=18)
plt.ylabel('Energía [ε]', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)

# Marcar el mínimo teórico
min_x = 2**(1/6)  # ≈1.122
min_y = -1.0
plt.axvline(x=min_x, color='r', linestyle='--', alpha=0.5)
plt.axhline(y=min_y, color='r', linestyle='--', alpha=0.5)
plt.axhline(y=0, color='g', linestyle='--', alpha=0.8, label="E = 0")
plt.plot(min_x, min_y, 'ro', label=f'Mínimo: ({min_x:.3f}σ, -ε)')

plt.legend(fontsize=15)
plt.ylim(-1.5, 1.0)  # Ajustar para ver el mínimo
plt.xlim(0.5, 3)  # Ajustar para ver el mínimo
plt.axis('equal')
plt.tick_params('both', labelsize=15)
plt.savefig("../Animaciones y graficas/Potencial.png")
plt.show()