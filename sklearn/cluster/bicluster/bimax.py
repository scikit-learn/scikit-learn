"""Implements the BiMax biclustering algorithm.

Authors : Kemal Eren
License: BSD 3 clause

"""
from abc import ABCMeta
from collections import defaultdict

import numpy as np

from sklearn.base import BaseEstimator, BiclusterMixin
from sklearn.externals import six

from .utils import get_indicators


def _precompute_neighbors(X):
    n_rows, n_cols = X.shape

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

    all_neighbors = {}
    for node in neighbors:
        all_neighbors[node] = neighbors[node] | second_neighbors[node]
        all_neighbors[node].discard(node)
    return all_neighbors


def _find_pivot(nodes, degrees, lim):
    pivot = -1
    pivot_degree = -1
    for n in nodes:
        degree = degrees[n]
        if degree > pivot_degree:
            pivot = n
            pivot_degree = degree
        if degree == lim:
            break
    return pivot, pivot_degree


def _find_bicliques(X):
    """Find all bicliques in a bipartite graph formed from array X.

    The graph has m+n nodes, where X has shape (m, n). If X[i, j]
    is nonzero, there is an edge between i and j.

    The bicliques are enumerated by a modification of the
    Bron-Kerbosch algorithm. This function is based on networkx's
    implementation.

    """
    neighbors = _precompute_neighbors(X)
    degrees = {n: len(nbrs) for n, nbrs in neighbors.iteritems()}
    pivot = max(degrees, key=degrees.get)

    candidates = set(neighbors)
    small_candidates = set(candidates - neighbors[pivot])
    done = set()
    stack = []
    biclique_so_far = []

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
        pivot_done, degree_done = _find_pivot(new_done, degrees,
                                              n_new_candidates)
        if degree_done == n_new_candidates:
            # shortcut: this part of tree already searched
            biclique_so_far.pop()
            continue
        pivot_cand, degree_cand = _find_pivot(new_candidates, degrees,
                                              n_new_candidates - 1)
        if degree_cand > degree_done:
            pivot = pivot_cand
        else:
            pivot = pivot_done

        # save search tree for later backout
        stack.append((candidates, done, small_candidates))
        candidates = new_candidates
        done = new_done
        small_candidates = candidates - neighbors[pivot]


class BiMax(six.with_metaclass(ABCMeta, BaseEstimator,
                               BiclusterMixin)):
    """Method to find all maximal biclusters in a boolean array."""

    def __init__(self):
        pass

    def fit(self, X):
        """Creates a biclustering for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        """
        # TODO: check X
        result = list(_find_bicliques(X))

        all_rows = []
        all_cols = []
        n_rows = X.shape[0]
        for nodes in result:
            rows = list(n for n in nodes if n < n_rows)
            cols = list(n - n_rows for n in nodes if n >= n_rows)
            if not rows or not cols:
                continue
            rows, cols = np.array(rows), np.array(cols)
            row_idx, col_idx = get_indicators(rows, cols, X.shape)
            all_rows.append(row_idx)
            all_cols.append(col_idx)
        self.rows_ = np.vstack(all_rows)
        self.columns_ = np.vstack(all_cols)
