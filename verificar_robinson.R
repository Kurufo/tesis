library('seriation')

# Definir los valores de la matriz en un vector
t_valores <- c(0, 3, 2, 4, 6,
               3, 0, 2, 4, 6,
               2, 2, 0, 1, 5,
               4, 4, 1, 0, 3,
               6, 6, 5, 3, 0)

#(0 2 3
# 2 0 2
# 3 2 0)

# Crear la matriz de 5x5 llenando por filas (byrow = TRUE)
t_matriz <- matrix(t_valores, nrow = 5, byrow = TRUE)

# Mostrar la matriz
print(t_matriz)

# Verificar si es Robinson de disimilitud (anti=TRUE, pre=TRUE)
is.robinson(t_matriz, anti = TRUE, pre = TRUE)

# Definir los valores de la matriz en un vector
valores <- c(0, 2, 4, 6,
             2, 0, 1, 5,
             4, 1, 0, 3,
             6, 5, 3, 0)

# Crear la matriz de 4x4 llenando por filas (byrow = TRUE)
matriz <- matrix(valores, nrow = 4, byrow = TRUE)

# Mostrar la matriz
print(matriz)

# Verificar si es Robinson de disimilitud (anti=TRUE, pre=TRUE)
is.robinson(matriz, anti = TRUE, pre = TRUE)

# Definir los valores de la matriz en un vector
valores <- c(0, 4, 4, 5, 5,
             4, 0, 3, 4, 5,
             4, 3, 0, 2, 5,
             5, 4, 2, 0, 1,
             5, 5, 5, 1, 0)

# Crear la matriz de 4x4 llenando por filas (byrow = TRUE)
matriz <- matrix(valores, nrow = 5, byrow = TRUE)

# Mostrar la matriz
print(matriz)

# Verificar si es Robinson de disimilitud (anti=TRUE, pre=TRUE)
is.robinson(matriz, anti = TRUE, pre = TRUE)