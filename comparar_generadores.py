# cosas.py
import numpy as np
import time
import utils
import subprocess
import sys
import csv
import statistics
from pulp import *
from generar_matriz_centros import generar_matriz

############################
# Argumentos
############################
if len(sys.argv) < 4:
    print("Uso: python cosas.py <tamaño_matriz> <n_iteraciones> <-r|-p>")
    sys.exit(1)

n_size = int(sys.argv[1])
n_iter = int(sys.argv[2])
mode = sys.argv[3]  # -r o -p

if mode not in ["-r", "-p"]:
    print("El tercer argumento debe ser -r (usar R) o -p (usar Python).")
    sys.exit(1)

############################
# Resultados
############################
resultados = []

############################
# Loop de iteraciones
############################
for rep in range(n_iter):
    print(f"\n=== Iteración {rep+1}/{n_iter}, tamaño {n_size}, modo {mode} ===")

    # Generar matriz
    if mode == "-r":
        subprocess.call(["/usr/bin/Rscript", "create_Robinson.r", str(n_size)])
        input_matrix = np.genfromtxt("random_robinson_matrix.csv", delimiter=',')
        metodo = "R"
        # Compute distance matrix
        distance_matrix = input_matrix.max() - input_matrix

        # Compute left/right centers
        max_closer = utils.compute_max_closer(distance_matrix)
        
    else:
        max_closer = generar_matriz(n_size)
        metodo = "Python"

    

    ##########################
    # Aux problem
    ##########################
    tic_aux = time.time()
    aux_problem = LpProblem("Aux_problem", LpMinimize)
    aux_distance = LpVariable.dicts("aux_distance", range(n_size), cat='Continuous')
    aux_leg = LpVariable.dicts("aux_leg", range(n_size), cat='Continuous')
    aux_problem += 0  # objetivo dummy

    two_points = [(i, j) for i in range(n_size) for j in range(n_size) if i != j]

    for (i, j) in two_points:
        if i < j:
            if (i == j-1) or ((max_closer[i, j-1] != max_closer[i, j]) and ((i == 0) or (max_closer[i-1, j] != max_closer[i, j]))):
                aux_problem += aux_distance[i] + aux_distance[j] - 2*aux_distance[max_closer[i, j]] - aux_leg[i] + aux_leg[j] >= 1
        elif i > j:
            if (i == j+1) or ((max_closer[i, j+1] != max_closer[i, j]) and ((i == (n_size-1)) or (max_closer[i+1, j] != max_closer[i, j]))):
                aux_problem += -aux_distance[i] - aux_distance[j] + 2*aux_distance[max_closer[i, j]] - aux_leg[i] + aux_leg[j] >= 1

    aux_status = aux_problem.solve(PULP_CBC_CMD(msg=0))
    toc_aux = time.time()
    aux_constraints = len(aux_problem.constraints)

    ##########################
    # Main problem
    ##########################
    tic = time.time()
    problem = LpProblem("Main_problem", LpMinimize)

    distance = LpVariable.dicts("distance", range(n_size), cat='Continuous', lowBound=0)
    leg = LpVariable.dicts("leg", range(n_size), cat='Continuous', lowBound=0)

    problem += lpSum([0*leg[i] for i in range(n_size)])

    for (i, j) in two_points:
        if (i < j-1):
            problem += distance[i] + distance[j] - 2*distance[max_closer[i, j]] - leg[i] + leg[j] >= 1
        elif (i > j+1):
            problem += -distance[i] - distance[j] + 2*distance[max_closer[i, j]] - leg[i] + leg[j] >= 1
        elif (i == j-1):
            problem += distance[i] + distance[j] - 2*distance[max_closer[i, j]] - leg[i] + leg[j] >= 1
        elif (i == j+1):
            problem += -distance[i] - distance[j] + 2*distance[max_closer[i, j]] - leg[i] + leg[j] >= 1

    main_status = problem.solve(PULP_CBC_CMD(msg=0))
    toc = time.time()
    main_constraints = len(problem.constraints)

    print(f"Aux problem: {LpStatus[aux_status]}, {toc_aux - tic_aux:.4f} s, {aux_constraints} restricciones")
    print(f"Main problem: {LpStatus[main_status]}, {toc - tic:.4f} s, {main_constraints} restricciones")

    # Guardar resultados en memoria
    resultados.append({
        "iteracion": rep + 1,
        "tamano": n_size,
        "metodo": metodo,
        "aux_status": LpStatus[aux_status],
        "main_status": LpStatus[main_status],
        "aux_time": toc_aux - tic_aux,
        "main_time": toc - tic,
        "aux_constraints": aux_constraints,
        "main_constraints": main_constraints
    })

############################
# Guardar resultados detallados a CSV
############################
detallado_filename = f"resultados_{metodo}_{n_size}.csv"
with open(detallado_filename, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=resultados[0].keys())
    writer.writeheader()
    writer.writerows(resultados)

############################
# Estadísticas globales
############################
aux_statuses = [r["aux_status"] for r in resultados]
main_statuses = [r["main_status"] for r in resultados]
aux_times = [r["aux_time"] for r in resultados]
main_times = [r["main_time"] for r in resultados]
aux_constraints_list = [r["aux_constraints"] for r in resultados]
main_constraints_list = [r["main_constraints"] for r in resultados]

summary = {
    "tamano": n_size,
    "metodo": metodo,
    "n_iter": n_iter,
    "aux_optimal_%": aux_statuses.count("Optimal") / n_iter * 100,
    "main_optimal_%": main_statuses.count("Optimal") / n_iter * 100,
    "aux_mean_time": statistics.mean(aux_times),
    "main_mean_time": statistics.mean(main_times),
    "aux_constraints_mean": statistics.mean(aux_constraints_list),
    "aux_constraints_std": statistics.pstdev(aux_constraints_list),
    "main_constraints_mean": statistics.mean(main_constraints_list),
    "main_constraints_std": statistics.pstdev(main_constraints_list)
}

# Guardar resumen
resumen_filename = f"resumen_{metodo}_{n_size}.csv"
with open(resumen_filename, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=summary.keys())
    writer.writeheader()
    writer.writerow(summary)

print("\n===== Estadísticas globales =====")
for k, v in summary.items():
    print(f"{k}: {v}")
print(f"\nResultados guardados en {detallado_filename}")
print(f"Resumen guardado en {resumen_filename}")
