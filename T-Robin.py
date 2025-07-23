import numpy as np

def min_distance_hyperedges(distance_matrix, labels):
    n = len(distance_matrix)
    hyperedges = []

    for i in range(n):
        row = distance_matrix[i]
        min_dist = np.min([row[j] for j in range(n) if j != i])
        edge = {labels[i]}
        for j in range(n):
            if i != j and distance_matrix[i][j] == min_dist:
                edge.add(labels[j])
        if len(edge) > 1 and edge not in hyperedges:
            hyperedges.append(edge)

    return hyperedges

def find_compatible_order_tree_and_hyperedges(distance_matrix, labels=None):
    X = labels if labels else [str(i) for i in range(len(distance_matrix))]
    E = min_distance_hyperedges(distance_matrix, X)

    X = list(X)
    k = len(X)
    X_k = X.copy()
    compatible_order = []
    tree_edges = []

    while k > 1:
        psi = {x: {y: 1 for y in X} for x in X}
        for x in X:
            psi[x][x] = 0

        for C in E:
            if len(set(X_k) & C) > 1:
                for x in C:
                    for y in X:
                        if y not in C:
                            psi[x][y] = 0

        found = False
        for x in X_k:
            for y in X_k:
                if psi[x][y] == 1:
                    x_k = x
                    compatible_order.append(x_k)

                    # Buscar padre real de x_k
                    for C in E:
                        if x_k in C:
                            candidates = list((C & set(X_k)) - {x_k})
                            if candidates:
                                parent = candidates[0]
                                tree_edges.append((parent, x_k))
                                break  # conectamos solo a uno

                    X_k.remove(x_k)
                    k -= 1
                    found = True
                    break
            if found:
                break

        if not found:
            print("No compatible order found.")
            return None, None, None

    compatible_order.append(X_k[0])
    return compatible_order, tree_edges, E



#A = np.array([
#    [0, 7,10, 5],
#    [7, 0, 7, 6],
#    [10,7, 0, 3],
#    [5, 6, 3, 0]
#])

A = np.array([
    [0, 5, 6, 1],
    [5, 0, 7, 2],
    [6, 7, 0, 3],
    [1, 2, 3, 0]
])
labels = ["a", "b", "c", "d"]

order, tree, hyperedges = find_compatible_order_tree_and_hyperedges(A, labels)

print("Compatible order:", order)
print("\nTree edges (parent → child):")
for p, c in tree:
    print(f"{p} → {c}")

print("\nHyperedges (sets):")
for h in hyperedges:
    print(h)
