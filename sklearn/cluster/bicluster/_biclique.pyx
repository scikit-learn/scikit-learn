# encoding: utf-8
# cython: cdivision=True
# cython: boundscheck=False
# cython: wraparound=False

from collections import defaultdict

import numpy as np

cimport numpy as np
cimport cython

np.import_array()

from libc.stdlib cimport malloc, free


cdef find_pivot(set nodes, set new_candidates, dict neighbors, long lim):
    cdef long pivot = -1
    cdef long pivot_degree = -1
    cdef long n
    cdef long degree
    for n in nodes:
        degree = len(new_candidates & neighbors[n])
        if degree > pivot_degree:
            pivot = n
            pivot_degree = degree
        if degree == lim:
            break
    return pivot, pivot_degree


cdef precompute_neighbors(char[:, :] X):
    cdef long n_rows
    cdef long n_cols
    n_rows = X.shape[0]
    n_cols = X.shape[1]

    cdef dict neighbors
    cdef long row
    cdef long col
    neighbors = {}
    for row in range(n_rows):
        neighbors[row] = set(np.nonzero(X[row])[0] + n_rows)
    for col in range(n_cols):
        neighbors[col + n_rows] = set(np.nonzero(X[:, col])[0])

    second_neighbors = defaultdict(set)
    for row in range(n_rows):
        for col in neighbors[row]:
            second_neighbors[row].update(neighbors[col])
    for col in range(n_cols):
        col = col + n_rows
        for row in neighbors[col]:
            second_neighbors[col].update(neighbors[row])

    cdef dict all_neighbors
    cdef long node
    all_neighbors = {}
    for node in neighbors:
        all_neighbors[node] = neighbors[node] | second_neighbors[node]
        all_neighbors[node].discard(node)
    return all_neighbors


def find_bicliques(char[:, :] X):
    """Find all bicliques in a bipartite graph formed from array X.

    The graph has m+n nodes, where X has shape (m, n). If X[i, j]
    is nonzero, there is an edge between i and j.

    The bicliques are enumerated by a modification of the
    Bron-Kerbosch algorithm. This function is based on networkx's
    implementation.

    """
    cdef dict neighbors = precompute_neighbors(X)
    cdef dict degrees = {n: len(nbrs) for n, nbrs in neighbors.iteritems()}
    cdef long pivot = max(degrees, key=degrees.get)

    cdef set candidates = set(neighbors)
    cdef set small_candidates = set(candidates - neighbors[pivot])
    cdef set done = set()
    cdef list stack = []
    cdef list biclique_so_far = []

    cdef set candidates_new
    cdef set done_ndw

    cdef long node
    cdef list result
    cdef long n_new_candidates
    cdef long pivot_done
    cdef long degree_done
    cdef long pivot_cand
    cdef long degree_cand

    while small_candidates or stack:
        try:
            # any nodes left to check?
            node = small_candidates.pop()
        except KeyError:
            # go to next in stack
            candidates, done, small_candidates = stack.pop()
            biclique_so_far.pop()
            continue

        # add next node to biclique
        biclique_so_far.append(node)
        candidates.remove(node)
        done.add(node)
        nn = neighbors[node]
        new_candidates = candidates & nn
        new_done = done & nn

        # check if we have more to search
        if not new_candidates:
            if not new_done:
                result = biclique_so_far[:]
                if len(result) > 1:
                    # found a biclique
                    yield result
            biclique_so_far.pop()
            continue

        n_new_candidates = len(new_candidates)
        pivot_done, degree_done = find_pivot(new_done,
                                             new_candidates,
                                             neighbors,
                                             n_new_candidates)
        if degree_done == n_new_candidates:
            # shortcut: this part of tree already searched
            biclique_so_far.pop()
            continue

        pivot_cand, degree_cand = find_pivot(new_candidates,
                                             new_candidates,
                                             neighbors,
                                             n_new_candidates - 1)
        if degree_cand > degree_done:
            pivot = pivot_cand
        else:
            pivot = pivot_done

        # save search tree for later backout
        stack.append((candidates, done, small_candidates))
        candidates = new_candidates
        done = new_done
        small_candidates = new_candidates - neighbors[pivot]
