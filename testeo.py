import numpy as np
from itertools import combinations

def find_compatible_order_with_full_hyperedges(distance_matrix, labels=None):
    """
    Parameters:
        distance_matrix (np.ndarray): symmetric distance matrix.
        labels (list): optional, list of labels for elements.

    Returns:
        order (list): compatible linear order.
        tree_edges (list of tuple): edges of the underlying tree.
        hyperedges (list of sets): list of hyperedges (as sets of labels or indices).
    """
    n = len(distance_matrix)
    if labels is None:
        labels = list(range(n))
    index_map = {i: labels[i] for i in range(n)}
    X = list(range(n))
    X_k = X.copy()
    k = len(X)
    order = [None] * k
    tree_edges = []
    hyperedges = []

    # Generate all hyperedges: for each node, group with nodes at same minimal distance (can be multiple sets)
    for i in range(n):
        dists = distance_matrix[i]
        sorted_indices = np.argsort(dists)
        prev_dist = -1
        for j in sorted_indices[1:]:  # skip self
            if prev_dist != -1 and not np.isclose(dists[j], prev_dist):
                break  # only group same distance levels
            prev_dist = dists[j]
            group = {i, j}
            for k2 in sorted_indices[1:]:
                if k2 != j and np.isclose(dists[k2], dists[j]):
                    group.add(k2)
            if len(group) > 1:
                group_frozenset = frozenset(group)
                if group_frozenset not in hyperedges:
                    hyperedges.append(group_frozenset)
            break  # only smallest distance per node group

    # Main loop
    while k > 1:
        # Initialize ψ: matrix of allowed eliminations
        psi = np.ones((n, n), dtype=bool)
        np.fill_diagonal(psi, 0)

        # Apply all hyperedges whose intersection with current X_k is > 1
        active_hyperedges = [C for C in hyperedges if len(set(C) & set(X_k)) > 1]
        for C in active_hyperedges:
            C = set(C)
            for x in C:
                for y in range(n):
                    if y not in C:
                        psi[x, y] = False

        # Find a pair (x,y) such that ψ_x(y) = 1
        found = False
        for x in X_k:
            for y in X_k:
                if x != y and psi[x, y]:
                    x_k = x
                    found = True
                    break
            if found:
                break

        if not found:
            print("No compatible order found.")
            return None, None, None

        # Append x_k to order
        order[k - 1] = index_map[x_k]

        # Find tree edge: connect to nearest node remaining in X_k
        remaining = list(set(X_k) - {x_k})
        if remaining:
            distances = [(distance_matrix[x_k][r], r) for r in remaining]
            min_dist, closest = min(distances)
            tree_edges.append((index_map[x_k], index_map[closest]))

        # Update X_k and k
        X_k.remove(x_k)
        k -= 1

    # Add the last remaining element
    order[0] = index_map[X_k[0]]

    # Convert hyperedges to label sets
    hyperedges_labeled = [set(index_map[i] for i in C) for C in hyperedges]

    return order, tree_edges, hyperedges_labeled

# Sample data for testing
A = np.array([
#     a   b   c   d   e
    [ 0,  1, 1, 10, 10],  # a
    [ 1,  0,  2, 10, 10],  # b
    [1, 2,  0,  3, 10],   # c
    [10,10, 3,  0,  10],    # d
    [10,10,10, 10,  0],     # e
])
labels = ["a", "b", "c", "d", "e"]

order, tree, hyperedges = find_compatible_order_with_full_hyperedges(A, labels)

print("Compatible order:", order)
print("\nTree edges (parent → child):")
for p, c in tree:
    print(f"{p} → {c}")

print("\nHyperedges (sets):")
for h in hyperedges:
    print(h)
