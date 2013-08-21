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


def _find_first_pivot(all_neighbors):
    """find first pivot (highest degree)"""
    max_degree = -1
    result = {}
    degrees = {}
    pivot_neighbors = set()  # handle empty graph
    for n, nbrs in all_neighbors.iteritems():
        degree = len(nbrs)
        degrees[n] = degree
        if degree > max_degree:
            result[n] = pivot_neighbors = nbrs
            max_degree = degree
        else:
            result[n] = nbrs
    return result, pivot_neighbors, degrees, max_degree


def _update_pivot(to_search, pivot_neighbors, new_candidates,
                  neighbors, all_degrees, lim):
    max_degree = -1
    n_candidates = len(new_candidates)
    for n in to_search:
        degree = n_candidates + all_degrees[n]
        if degree > max_degree:
            pivot_neighbors = new_candidates & neighbors[n]
            max_degree = degree
            if max_degree == lim:
                break
    return max_degree, pivot_neighbors


class BiMax(six.with_metaclass(ABCMeta, BaseEstimator,
                               BiclusterMixin)):
    """Method to find all maximal biclusters in a boolean array."""

    def __init__(self, random_state=None):
        self.random_state = random_state

    def _find_bicliques(self, X):
        """Find all bicliques in a bipartite graph formed from array X.

        The graph has m+n nodes, where X has shape (m, n). If X[i, j]
        is nonzero, there is an edge between i and j.

        The bicliques are enumerated by a modification of the
        Bron-Kerbosch algorithm. This function is based on networkx's
        implementation.

        """
        neighbors = _precompute_neighbors(X)
        neighbors, pivot_neighbors, degrees, max_degree = \
            _find_first_pivot(neighbors)
        candidates = set(neighbors)
        small_candidates = set(candidates - pivot_neighbors)
        done = set()
        stack = []
        biclique_so_far = []

        while small_candidates or stack:
            try:
                # any nodes left to check?
                n = small_candidates.pop()
            except KeyError:
                # go to next in stack
                candidates, done, small_candidates = stack.pop()
                biclique_so_far.pop()
                continue

            # add next node to biclique
            biclique_so_far.append(n)
            candidates.remove(n)
            done.add(n)
            nn = neighbors[n]
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

            n_candidates = len(new_candidates)

            # find pivot node; look in done nodes first
            max_degree_done, pivot_done_neighbors = \
                _update_pivot(new_done, pivot_neighbors,
                              new_candidates, neighbors, degrees,
                              n_candidates)

            # shortcut: this part of tree already searched
            if max_degree_done == n_candidates:
                biclique_so_far.pop()
                continue

            # still finding pivot node; look in candidates nodes
            max_degree, pivot_neighbors = \
                _update_pivot(new_candidates, pivot_neighbors,
                              new_candidates, neighbors,
                              degrees,
                              n_candidates - 1)

            # choose best pivot from ``candidates`` and ``done``
            if max_degree_done > max_degree:
                pivot_neighbors = pivot_done_neighbors

            # save search tree for later backout
            stack.append((candidates, done, small_candidates))
            candidates = new_candidates
            done = new_done
            small_candidates = candidates - pivot_neighbors

    def fit(self, X):
        """Creates a biclustering for X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)

        """
        # TODO: check X
        result = list(self._find_bicliques(X))

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
