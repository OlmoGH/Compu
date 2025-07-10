import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl

viridis = mpl.colormaps['viridis'].resampled(9)

datos = pd.read_csv("../Simulacion/Datos/Presiones2.txt", sep='\t', header=0, names=['N', 'L', 'Temp', 'p'])

T_range = [0, 2, 4, 8, 12, 17.98, 24, 31, 40]
L_range = [20, 30, 40, 50, 60, 70, 80, 90, 100]
N_range = [20, 30, 40, 50, 60, 70, 80, 90, 100]

Temperatura = (T_range, datos.Temp, 'Temp')
Numero = (N_range, datos.N, 'N')
Lado = (L_range, datos.L, 'L')

variable_dependiente = Temperatura
magitud_principal = Numero
magitud_secundaria = Lado

main_range = magitud_principal[0]
sec_range = magitud_secundaria[0]

main_variable = magitud_principal[1]
sec_variable = magitud_secundaria[1]

char = variable_dependiente[2]
indicador = magitud_principal[2]
mag_fija = magitud_secundaria[2]
valor_fijo = sec_range[5]

fig, ax = plt.subplots(figsize=(8, 5))

if (variable_dependiente[2] == 'L'):   
    for i, main_index in enumerate(main_range):

        mask = (np.abs(sec_variable - valor_fijo) < 1.2) & (datos.Temp < 1000) & (np.abs(main_variable - main_index) < 1.2)

        presiones_filtradas = datos.loc[mask, 'p'].to_numpy()

        variables_filtradas = datos.loc[mask, char].to_numpy()

        indices = np.argsort(variables_filtradas)
        x = 1 / variables_filtradas[indices] ** 2
        y = presiones_filtradas[indices]

        ax.plot(x, y, '-', color=viridis(i), label=indicador + f" = {main_index}")
        ax.set_xlabel(r"$1/L^2$", fontsize=12)
        ax.set_ylabel("Presión", fontsize=12)
        ax.set_title(fr"Presión frente a $1/L^2$ para {mag_fija} = {valor_fijo}", fontsize=12)
        ax.legend()
else:
    for i, main_index in enumerate(main_range):

        mask = (np.abs(sec_variable - valor_fijo) < 1.2) & (datos.Temp < 1000) & (np.abs(main_variable - main_index) < 1.2)

        presiones_filtradas = datos.loc[mask, 'p'].to_numpy()

        variables_filtradas = datos.loc[mask, char].to_numpy()

        indices = np.argsort(variables_filtradas)
        x = variables_filtradas[indices]
        y = presiones_filtradas[indices]

        ax.plot(x, y, '-', color=viridis(i), label=indicador + f" = {main_index}")
        ax.set_xlabel(fr"${char}$", fontsize=12)
        ax.set_ylabel("Presión", fontsize=12)
        ax.set_title(f"Presión frente a {char} para {mag_fija} = {valor_fijo}", fontsize=12)
        ax.legend()


temps = datos.loc[np.abs(datos.Temp - 40) < 5, 'Temp']

print(temps.mean())

plt.savefig(f"../Animaciones y graficas/Presión frente a {char} para {mag_fija} = {valor_fijo}.png")
plt.show()