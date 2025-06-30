import numpy as np
from numba import jit
import time

inicio = time.time()
@jit(nopython=True)
def multiply_by_gravitational_constant_withjit(matrix):
    g = 6.6743e-10
    result = np.empty_like(matrix)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            result[i,j] = matrix[i,j] * g
    return result



def multiply_by_gravitational_constant(matrix):
    g = 6.6743e-10
    result = np.empty_like(matrix)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            result[i,j] = matrix[i,j] * g
    return result

inicio3 = time.time()
matriz = np.random.rand(50000, 50000)
fin3 = time.time()

inicio2 = time.time()
matrizg1 = multiply_by_gravitational_constant(matriz)
fin2 = time.time()

inicio1 = time.time()
matrizg2 = multiply_by_gravitational_constant_withjit(matriz)
fin1 = time.time()

tiempo1 = fin1 - inicio1
tiempo2 = fin2 - inicio2
fin = time.time()
print(f"El tiempo que tarda en ejecutarse con el compilador es {tiempo1} s y sin el compilador es {tiempo2} s. El timepo que tarda en crear la matriz es {fin3-inicio3} s. El programa se ejecuta en {fin - inicio} s")