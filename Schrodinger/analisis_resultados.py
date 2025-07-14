import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json

# Cargar el archivo JSON
with open('results.json', 'r') as f:
    data = json.load(open('results.json', 'r'))

# Convertir a DataFrame
df = pd.json_normalize(data['results'])
print(df)

# Verificar las columnas disponibles
print("Columnas disponibles:", df.columns.tolist())

# Extraer parámetros del nombre del benchmark
df[['opt', 'threads', 'N', 'steps']] = df['command'].str.extract(
    r'O(\d)-Threads(\d+)-N(\d+)-Steps(\d+)'
)

df['opt'] = df['opt'].astype(int)       # Nivel de optimización
df['threads'] = df['threads'].astype(int)  # Número de hilos
df['N'] = df['N'].astype(int)           # Tamaño del problema
df['steps'] = df['steps'].astype(int)   # Pasos de ejecución

# Convertir a numérico
df = df.astype({'opt': int, 'threads': int, 'N': int, 'steps': int})

# Tiempo frente a optimizacion
mask = (df['steps'] == df['steps'].median()) & (df['threads'] == 4)
print(df[mask])
for i in df['N'].unique():
    mask = (df['steps'] == df['steps'].median()) & (df['threads'] == 4)
    mask = mask & (df['N'] == i)
    plt.errorbar(x=df[mask]['opt'], y=df[mask]['mean'], yerr=df[mask]['stddev'] ,label=f"N: {i}")

plt.legend(loc='upper left')
plt.title("Tiempo de ejecucion frente a optimización")
plt.xlabel("Nivel de optimización")
plt.ylabel("Tiempo de ejecucion (s)")
plt.show()

#Tiempo frente a hilos
mask = (df['steps'] == df['steps'].median()) & (df['opt'] == 2)
print(df[mask])
for i in df['N'].unique():
    mask = (df['steps'] == df['steps'].median()) & (df['opt'] == 2)
    mask = mask & (df['N'] == i)
    plt.errorbar(x=df[mask]['threads'], y=df[mask]['mean'], yerr=df[mask]['stddev'] ,label=f"N: {i}")

plt.legend(loc='upper left')
plt.title("Tiempo de ejecucion frente a hilos")
plt.xlabel("Número de hilos")
plt.ylabel("Tiempo de ejecucion (s)")
plt.show()