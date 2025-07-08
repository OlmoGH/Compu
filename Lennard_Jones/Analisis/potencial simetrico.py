import numpy as np
import matplotlib.pyplot as plt

# Definir rango alrededor del mínimo (1.122σ)
x = np.linspace(-10, 10, 10000)  # Centrado en ~1.122σ

# Potencial de Lennard-Jones (con σ=1)
x_despl = x - 1
def potencial(x):
    return  4 * (np.power(x, -12) - np.power(x, -6))

a = 2 ** (1/6)
y = potencial(x-a) + potencial(x+a)

y_aprox = -2 + 36 * 2 ** (2/3) * np.pow(x, 2)

# Graficar
plt.figure(figsize=(8, 8))
plt.plot(x, y_aprox, color='gray', linestyle='--', linewidth=2, label="Aproximación armónica")
plt.plot(x, y, 'b-', linewidth=2)
plt.title('Potencial simétrico', fontsize=25)
plt.xlabel('r/σ', fontsize=18)
plt.ylabel('Energía [ε]', fontsize=18)
plt.grid(True, linestyle='--', alpha=0.7)

# Marcar el mínimo teórico
min_x = 0  # ≈1.122
min_y = -2.0
plt.axvline(x=min_x, color='r', linestyle='--', alpha=0.5)
plt.axhline(y=min_y, color='r', linestyle='--', alpha=0.5)
plt.axhline(y=0, color='g', linestyle='--', alpha=0.8, label="E = 0")
plt.plot(min_x, min_y, 'ro', label=f'Mínimo: (0, -2ε)')

plt.legend(fontsize=15)
plt.ylim(-2.5, 1)  # Ajustar para ver el mínimo
plt.xlim(-.5, .5)  # Ajustar para ver el mínimo
plt.tick_params('both', labelsize=15)
plt.savefig("../Animaciones y graficas/Potencial simetrico.png")
plt.show()