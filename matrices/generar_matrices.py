import numpy as np
from itertools import product
import os

def leer_matrices_desde_archivo(nombre_archivo):
    with open(nombre_archivo, 'r') as f:
        bloques = f.read().strip().split('\n\n')
    matrices = []
    for bloque in bloques:
        filas = [list(map(int, linea.strip().split())) for linea in bloque.strip().split('\n')]
        matrices.append(np.array(filas, dtype=int))
    return matrices

def escribir_matrices_a_archivo(matrices, nombre_archivo):
    with open(nombre_archivo, 'w') as f:
        for matriz in matrices:
            for fila in matriz:
                f.write(' '.join(map(str, fila)) + '\n')
            f.write('\n')  # Separador entre matrices

def es_valida(nueva_fila, matriz, n):
    for i in range(n):
        if nueva_fila[i] < max(matriz[i][n-1], i + 1) or nueva_fila[i] >= n+1:
            return False
    for j in range(1, n):
        if nueva_fila[j] < nueva_fila[j-1]:
            return False
    return True

def construir_matriz_extendida(matriz, nueva_fila, n):
    nueva_matriz = np.zeros((n+1, n+1), dtype=int)
    nueva_matriz[:n, :n] = matriz
    for i in range(n):
        nueva_matriz[i, n] = nueva_fila[i]
        nueva_matriz[n, i] = nueva_fila[i]
    for i in range(n+1):
        for j in range(i):
            nueva_matriz[i, j] = nueva_matriz[j, i] + 1
    return nueva_matriz

def extender_matriz(matriz):
    n = matriz.shape[0]
    posibilidades = []
    rangos = [range(max(matriz[i][n-1], i + 1), n+1) for i in range(n)]
    for nueva_fila in product(*rangos):
        if es_valida(nueva_fila, matriz, n):
            nueva_matriz = construir_matriz_extendida(matriz, nueva_fila, n)
            posibilidades.append(nueva_matriz)
    return posibilidades

def generar_siguiente_nivel(nombre_archivo_entrada, nombre_archivo_salida):
    matrices_n = leer_matrices_desde_archivo(nombre_archivo_entrada)
    todas_matrices = []
    for matriz in matrices_n:
        nuevas = extender_matriz(matriz)
        todas_matrices.extend(nuevas)
    escribir_matrices_a_archivo(todas_matrices, nombre_archivo_salida)
    print(f"Generadas {len(todas_matrices)} matrices en '{nombre_archivo_salida}'")
    
    
    
generar_siguiente_nivel("matrices_8.txt", "matrices_9.txt")



