from itertools import permutations
import numpy as np

def generar_permutaciones_reducidas(n):
    k = n*(n-1)//2
    valores = list(range(1, k))  # omitimos el k (mayor valor)
    return permutations(valores), k

def permutacion_a_matriz_con_max(permutacion, n, max_valor):
    matriz = np.zeros((n, n), dtype=int)
    idx = 0
    posiciones = []
    for i in range(n):
        for j in range(i+1, n):
            posiciones.append((i, j))
    posiciones.remove((0, n-1))  # dejamos el mayor para la esquina

    matriz[0][n-1] = max_valor
    matriz[n-1][0] = max_valor

    for (i, j) in posiciones:
        valor = permutacion[idx]
        matriz[i][j] = valor
        matriz[j][i] = valor
        idx += 1

    return matriz
