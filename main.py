import subprocess
from permutaciones import generar_permutaciones_reducidas, permutacion_a_matriz_con_max
import utils
import numpy as np
import time
import pandas as pd
import re

#Solve LP
from pulp import * # solve LP

from scipy.optimize import linprog

def es_robinson_subprocess(matriz):
    matriz_str = ",".join(map(str, matriz.flatten()))
    n = matriz.shape[0]
    result = subprocess.run(
        [r"C:\Program Files\R\R-4.2.1\bin\Rscript.exe", 'check_robinson.R', matriz_str, str(n)],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print("Error en R:", result.stderr)
        return None
    output = result.stdout.strip()
    return output == "TRUE"

def analizar_matrices_reducidas(n):
    n_exitos=0
    n_fallas=0
    permutaciones, max_valor = generar_permutaciones_reducidas(n)
    total = 0
    robinson_count = 0
    for perm in permutaciones:
        matriz = permutacion_a_matriz_con_max(perm, n, max_valor)
        if es_robinson_subprocess(matriz):
            robinson_count += 1
            print("Cantidad de Robinsons: ", robinson_count)
            input_matrix = matriz

            distance_matrix=input_matrix 

            max_closer=utils.compute_max_closer(distance_matrix)
                
            #######################
            # Problem description #
            #######################
            n = len(input_matrix)
            ROWS = COLS = range(n)
            tol=0.11
            
            # Node labels
            nodelabel = ["S"+str(i) for i in range(n)]
            
            # Create the n-tuples  created, with the row and column index
            four_points = [ # order i<j<k<l 
                (i,j,k,l) 
                for l in range(3,n)
                for k in range(2,l)
                for j in range(1,k)
                for i in range(j)
            ]
            
            three_points = [
                (i,j,k) 
                for k in range(n)
                for j in range(n)
                for i in range(n)
            ]
            
            # Pairs
            two_points = [ 
                (i,j) 
                for j in range(n)
                for i in range(n)
            ]
            
            #### AUX PROBLEM #################################################################################################################################################################################
            # Time aux
            #Aux problem                                                                                                                                                                                   #
            aux_problem = LpProblem("A_has_a_positive_linear_combination", LpMinimize)                                                                                                                     #
            aux_distance = LpVariable.dicts("aux_distance", ROWS, cat='Continuous')                                                                                                                        #
            aux_leg      = LpVariable.dicts("aux_leg", ROWS, cat='Continuous')                                                                                                                             #
            aux_problem += 0 # feasibility                                                                                                                                                                 #
            # Constraints for each pair                                                                                                                                                                    #
            for (i,j) in two_points:                                                                                                                                                                       #
                if (i<j-1) :                                                                                                                                                                               #
                    aux_problem += LpConstraint( aux_distance[i] + aux_distance[j] - 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j))  #
                elif (i>j+1) :                                                                                                                                                                             #
                    aux_problem += LpConstraint( -aux_distance[i] - aux_distance[j] + 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j)) #
                elif (i==j-1) :                                                                                                                                                                            #
                    aux_problem += LpConstraint( aux_distance[i] + aux_distance[j] - 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j))  #
                elif (i==j+1) :                                                                                                                                                                            #
                    aux_problem += LpConstraint( -aux_distance[i] - aux_distance[j] + 2*aux_distance[max_closer[i,j]] - aux_leg[i] + aux_leg[j], sense=LpConstraintGE, rhs=1 , name='A'+str(i)+'_'+str(j)) #
         
            # Solve the aux_problem                                                                                                                                                                  #
            aux_problem.solve(PULP_CBC_CMD(msg=0))
            
            # Create Matrix A
            # create DataFrame
            list_constraints =  list(aux_problem.constraints)
            list_variables = list()
            for var in aux_problem.variables():
                list_variables.append(var.name)
            # remove dummy variable
            list_variables=list_variables[1:]
            
            # create empty matrix A
            A = pd.DataFrame(np.zeros(shape=(len(list_constraints),len(list_variables))), columns=list_variables, index=list_constraints)
            # problem dictionary
            dict = aux_problem.to_dict()
            # fill matrix A
            for constraint in dict['constraints']:
                for coefficient in constraint['coefficients']:
                    A.loc[constraint['name'], coefficient['name']]=coefficient['value']

            
            # Create matrix B
            l2 = int(len(list_variables)/2) # two parts
            columns_B = ['d+l_'+str(i) for i in range(l2)]+['d-l_'+str(i) for i in range(l2)]
            B = pd.DataFrame(np.zeros(shape=(len(list_constraints),len(list_variables))), columns=columns_B, index=list_constraints)
            
            for i in range(l2):
                B.iloc[:,i] = (A.iloc[:,i]+A.iloc[:,i+l2])/2
                B.iloc[:,i+l2] = (A.iloc[:,i]-A.iloc[:,i+l2])/2

            
            
            # Order rows of B 
            # add type constraint
            
            type_constraint = []
            type_constraint_order = []
            groups_size =[0,0,0]
            
            for (i,j) in two_points:
                if i != j :
                    type = 0
                    order = i
                    if (i<j) :                                                                                                                                                                               #
                        if max_closer[i,j] == i:
                            type = 0
                            order = i
                        else :
                            type = 1
                            order = j
                    elif (i>j) :                                                                                                                                                                         
                        if max_closer[i,j] == i:
                            type = 0
                            order = i
                        else :
                            type = 2
                            order = max_closer[i,j]
                    type_constraint.append(int(type))
                    type_constraint_order.append(int(order))
                    groups_size[type] = groups_size[type]+1
            
            # Add type and order to sort the rows
            B['type_constraint']=type_constraint
            B['type_constraint_order']=type_constraint_order
            
            restricciones=B[B['type_constraint']>0]
            restricciones=restricciones.sort_values(by="type_constraint")

            
            C = pd.DataFrame(np.zeros(shape=(restricciones.shape[0],len(list_variables))), columns=columns_B, index=restricciones.index)
            flag=0
            maxi_cent_der=0
            maxi_cent_izq=0
            max_q=0
            min_p=l2

            
            # Iterar sobre los índices de las filas en restricciones
            for idx in restricciones.index:
                # Extraer p y q del índice utilizando la expresión regular
                match = re.match(r"A(\d+)_(\d+)", idx)
                if match:
                    p, q = int(match.group(1)), int(match.group(2))
                    max_c = int(max_closer[p, q])  # Obtener el valor de max_closer para (p, q)
                   # print("p: ",p, ", q: ",q, ", max-closer: ", max_c)
                    # Calcular los rangos de columnas en base a p y q
                    




                    if p < q:
                        maxi_cent_izq=max(maxi_cent_izq,max_c)
                        # Si p < q, poner 1s desde p+1 hasta max_c y -1s desde max_c+1 hasta q
                        C.iloc[C.index.get_loc(idx), max_c+1:q+1] = 1
                        C.iloc[C.index.get_loc(idx), l2+p+1:l2+max_c+1] = -1
        #                 C.loc[idx, max_c+1:q+1] = 1  # np.arange funciona como índice para asignación
        #                 C.loc[idx, l2+p+1:l2+max_c+1] = -1
                    else:
                        maxi_cent_der=max(maxi_cent_der, max_c)
                        # Si p > q, poner -1s desde p+1 hasta max_c y 1s desde max_c+1 hasta q
                        C.iloc[C.index.get_loc(idx), max_c+1:p+1] = -1
                        C.iloc[C.index.get_loc(idx), l2+q+1:l2+max_c+1] = 1
        #                 C.loc[idx, max_c+1:p+1] = -1  # np.arange funciona como índice para asignación
        #                 C.loc[idx, l2+q+1:l2+max_c+1] = 1

                        
            C.astype(int).to_csv('C.csv')
            
            # Mostrar el resultado final de C
            #print("Matriz C")
            #print(C)
            
            coeff = np.array([i for i in range(l2)]+[i for i in range(l2)])
            type_constraint=[int(x) for x in list(restricciones['type_constraint'])]
            restricciones.drop(columns=['type_constraint', 'type_constraint_order'],inplace=True)
            restricciones.to_numpy()
            result = np.dot(restricciones,coeff)
            restricciones=pd.DataFrame(restricciones, columns=columns_B, index= C.index)
            restricciones['type_constraint']=type_constraint
            
            
            restricciones['suma_res']=result
            
            matriz_distancias=distance_matrix
            matriz_distancias=pd.DataFrame(matriz_distancias)
                    
            # Seleccionamos solo las columnas de variables en C
            C['type_constraint']=type_constraint
            C_vars = C.drop(columns=['type_constraint']).to_numpy()
            
            # Crear la función objetivo
            c = np.ones(C_vars.shape[1])
            
            # Crear las restricciones de desigualdad
            # Filtrar las filas con type_constraint == 1 (queremos suma > 0)
            A_ineq_1 = -C_vars[C['type_constraint'] == 1]  # Negamos para convertir en <=
            b_ineq_1 = -np.ones(A_ineq_1.shape[0])*l2  # Suma > 0 se convierte en -suma < -1
            
            # Filtrar las filas con type_constraint == 2 (queremos suma < 0)
            A_ineq_2 = -C_vars[C['type_constraint'] == 2]  # Mantenemos para suma < 0
            b_ineq_2 = -np.ones(A_ineq_2.shape[0])*l2  # Suma < 0 se representa como suma <= 1
            
            # Concatenamos las restricciones de desigualdad
            A_ineq = np.vstack([A_ineq_1, A_ineq_2])
            b_ineq = np.concatenate([b_ineq_1, b_ineq_2])
            C.drop(columns=['type_constraint'])
            # Resolver el sistema
            resulta = linprog(c, A_ub=A_ineq, b_ub=b_ineq, method='highs')
            
            # Verificamos el resultado
            if resulta.success:
                n_exitos+=1
            else:
                n_fallas+=1



        total += 1
        if total % 1000 == 0:
            print(f"Procesadas: {total} matrices")

    print(f"De {total} matrices, {robinson_count} son Robinson")
    print(f"De {robinson_count} espacios Robinson, en {n_exitos} se encuentra un dibujo valido")
n = 5
analizar_matrices_reducidas(n)
