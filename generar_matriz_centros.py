import numpy as np

def generar_matriz(n, seed=None):
    if seed is not None:
        np.random.seed(seed)

    M = np.zeros((n, n), dtype=int)

    # Diagonal principal
    for i in range(n):
        M[i, i] = i

    # Diagonal justo encima de la principal
    for i in range(n - 1):
        M[i, i + 1] = i

    # Parte superior derecha (más allá de la diagonal justo encima)
    for d in range(2, n):  # distancia desde la diagonal principal
        for i in range(n - d):
            j = i + d
            a = M[i, j - 1]
            b = M[i + 1, j]
            low = min(a, b)
            high = max(a, b)
            M[i, j] = np.random.randint(low, high + 1)

    # Parte inferior izquierda (simétrica + 1)
    for i in range(n):
        for j in range(i):
            M[i, j] = M[j, i] + 1

    return M
