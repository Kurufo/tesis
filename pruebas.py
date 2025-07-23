import utils
import numpy as np
import time
import pandas as pd
import re
from scipy.optimize import linprog

# Arguments
from argparse import ArgumentParser

#Solve LP
from pulp import * # solve LP

# Cargar la matriz C desde el archivo CSV
C = pd.read_csv('C_construido.csv', index_col=0)


# Obtener las filas de tipo 2
filas_tipo_2 = C[C['type_constraint'] == 2].index
print("Filas de tipo 2 presentes en la matriz C:", list(filas_tipo_2))

l2=(C.shape[1]-1)/2

# Función para resolver el sistema de inequaciones
def resolver_sistema(C_matrix):
    type_constraint = C_matrix['type_constraint'].astype(int).to_numpy()
    C_vars = C_matrix.drop(columns=['type_constraint']).to_numpy()

    # Crear la función objetivo
    c = np.ones(C_vars.shape[1])

    # Crear las restricciones de desigualdad
    A_ineq_1 = -C_vars[type_constraint == 1]
    b_ineq_1 = -np.ones(A_ineq_1.shape[0])*l2

    A_ineq_2 = -C_vars[type_constraint == 2]
    b_ineq_2 = -np.ones(A_ineq_2.shape[0])*l2

    A_ineq = np.vstack([A_ineq_1, A_ineq_2])
    b_ineq = np.concatenate([b_ineq_1, b_ineq_2])

    # Resolver el sistema
    return linprog(c, A_ub=A_ineq, b_ub=b_ineq, method='highs')

# Resolver el sistema para la matriz C original
resultado_original = resolver_sistema(C)
if resultado_original.success:
    solucion_original = resultado_original.x
    print("Solución para la matriz C original:", solucion_original)
else:
    print("No se pudo encontrar una solución para la matriz C original.")

# Almacenar soluciones de submatrices
soluciones_submatrices = {}

# Iterar sobre cada fila de tipo 2 para crear y resolver submatrices
for fila in filas_tipo_2:
    # Crear una submatriz eliminando la fila actual
    submatriz = C.drop(index=fila)

    # Resolver el sistema para la submatriz
    resultado_submatriz = resolver_sistema(submatriz)

    # Verificar si la solución es exitosa
    if resultado_submatriz.success:
        soluciones_submatrices[fila] = resultado_submatriz.x
        print(f"Solución para la submatriz sin la fila {fila}:", resultado_submatriz.x)
    else:
        soluciones_submatrices[fila] = None
        print(f"No se pudo encontrar una solución para la submatriz sin la fila {fila}.")

# Comparar soluciones
print("\nComparación de soluciones:")
print("Solución original:", solucion_original)
for fila, solucion in soluciones_submatrices.items():
    if solucion is not None:
        diferencia = solucion - solucion_original
        print(f"Diferencia en la solución al eliminar la fila {fila}:", diferencia)
    else:
        print(f"Submatriz sin la fila {fila} no tiene solución.")

# Obtener el vector de soluciones
soluciones = resultado_original.x

# Dividir en dos mitades
n = len(soluciones) // 2
primera_mitad = soluciones[:n]
segunda_mitad = soluciones[n:]

# Calcular el acumulado y sumar el índice en cada mitad
resultado_acumulado = np.concatenate([
    [sum(primera_mitad[:i+1]) + i for i in range(n)],
    [sum(segunda_mitad[:i+1]) + i for i in range(n)]
])

print(resultado_acumulado)
print(primera_mitad)
print(segunda_mitad)


# 1. Cargar la matriz B desde el archivo CSV
B = pd.read_csv("B.csv", index_col=0)  # index_col=0 asume que la primera columna es un índice y no parte de los datos

# Convertir la matriz B a un array de numpy para cálculos de álgebra lineal
B_matrix = B.to_numpy()


# 2. Calcular el producto punto del vector acumulado con la matriz B
producto_punto = np.dot(B_matrix,resultado_acumulado)
print("Producto punto:", producto_punto)

# 3. Verificar si todos los valores del producto punto son positivos
if np.all(producto_punto > 0):
    print("Todos los valores del producto punto son positivos.")
else:
   print("No todos los valores del producto punto son positivos.")
   
# Print solution from Aux
legs = (resultado_acumulado[:n]-resultado_acumulado[n:])/2
distance2 = (resultado_acumulado[:n]+resultado_acumulado[n:])/2  
distance = resultado_acumulado[:n]-legs  
print(f'distances : {distance}')
#print(f'distances2: {distance2}')
print(f'legs      : {legs}')
print()   

resultado_original = resolver_sistema(C)
