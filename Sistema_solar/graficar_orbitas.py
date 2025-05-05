
import matplotlib.pyplot as plt
import os

def leer_orbita(nombre_archivo):
    x, y = [], []
    with open(nombre_archivo, 'r') as f:
        for linea in f:
            datos = linea.strip().split()
            if len(datos) >= 2:
                x.append(float(datos[0]))
                y.append(float(datos[1]))
    return x, y

def graficar_orbitas():
    colores = ['b', 'g', 'r', 'c', 'm', 'y', 'orange', 'purple']
    nombres = ['Mercurio', 'Venus', 'Tierra', 'Marte', 'Júpiter', 'Saturno', 'Urano', 'Neptuno']

    plt.figure(figsize=(8, 8))
    for i in range(1, 9):
        archivo = f'planet_{i}.txt'
        if os.path.exists(archivo):
            x, y = leer_orbita(archivo)
            plt.plot(x, y, label=nombres[i - 1], color=colores[i - 1])
        else:
            print(f"Archivo no encontrado: {archivo}")

    plt.plot(0, 0, 'yo', label='Sol')  # Sol en el centro
    plt.xlabel("X (UA)")
    plt.ylabel("Y (UA)")
    plt.title("Órbitas planetarias alrededor del Sol")
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    graficar_orbitas()
