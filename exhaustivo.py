from itertools import permutations
from rpy2 import robjects
from rpy2.robjects import r, numpy2ri
import numpy as np

numpy2ri.activate()

def generar_permutaciones_reducidas(n):
    """
    Genera permutaciones de los valores excepto el máximo, que se colocará fijo en la esquina.
    """
    k = n*(n-1)//2
    valores = list(range(1, k))  # omitimos el k (mayor valor)
    return permutations(valores), k  # k = valor más grande

def permutacion_a_matriz_con_max(permutacion, n, max_valor):
    """
    Convierte una permutación en matriz simétrica nxn, 
    colocando el valor máximo fijo en la esquina superior derecha.
    """
    matriz = np.zeros((n, n), dtype=int)
    idx = 0
    # definimos las posiciones sobre la diagonal
    posiciones = []
    for i in range(n):
        for j in range(i+1, n):
            posiciones.append( (i,j) )
    # decidimos que la última posición sobre la diagonal (0,n-1) tendrá el valor máximo
    # eliminamos esa posición de las posiciones a permutar
    posiciones.remove( (0,n-1) )

    # colocamos el max_valor en la esquina
    matriz[0][n-1] = max_valor
    matriz[n-1][0] = max_valor

    # rellenamos el resto
    for (i,j) in posiciones:
        valor = permutacion[idx]
        matriz[i][j] = valor
        matriz[j][i] = valor
        idx += 1

    return matriz

def es_robinson(matriz):
    r('library(seriation)')
    r_matriz = robjects.r.matrix(robjects.IntVector(matriz.flatten()), nrow=matriz.shape[0], byrow=True)
    resultado = r['is.robinson'](r_matriz)
    return bool(resultado[0])

def analizar_matrices_reducidas(n):
    permutaciones, max_valor = generar_permutaciones_reducidas(n)
    total = 0
    robinson_count = 0
    for perm in permutaciones:
        matriz = permutacion_a_matriz_con_max(perm, n, max_valor)
        if es_robinson(matriz):
            robinson_count += 1
        total += 1
        if total % 1000 == 0:
            print(f"Procesadas: {total} matrices")

    print(f"De {total} matrices generadas, {robinson_count} son Robinson")

n = 4
analizar_matrices_reducidas(n)
