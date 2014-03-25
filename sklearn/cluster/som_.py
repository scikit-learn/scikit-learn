"""
 Self-organizing map
"""

# Authors: Sebastien Campion <sebastien.campion@inria.fr>
#          naught101 <naught101@gmail.com>
# License: BSD

from __future__ import division
import numpy as np
from ..base import BaseEstimator


def _unserialise_coordinate(serial, spacing):
    coord = []
    for i in range(len(spacing)):
        coord.append(int(np.floor(serial/spacing[i])))
        serial = int(serial % spacing[i])

    return(coord)


def _serialise_coordinate(coord, spacing):
    return(np.dot(coord, spacing))


def _generate_adjacency_matrix(grid_dimensions):
    """Generate an adjacency matrix for nodes of an n-dimensional rectangular
    grid with dimensions given by grid_dimensions.
    """
    n_centers = np.prod(grid_dimensions)

    # initialise to 0 for self-distances, infinity non-connections
    adjacency = np.empty((n_centers, n_centers), dtype=float)
    adjacency.fill(np.inf)
    np.fill_diagonal(adjacency, 0)

    spacing = [int(np.prod(grid_dimensions[(k+1):]))
               for k in range(len(grid_dimensions))]

    for i in range(n_centers):
        coord = _unserialise_coordinate(i, spacing)
        neighbours = []
        for d in range(len(grid_dimensions)):
            if coord[d] > 0:
                down = coord[:]
                down[d] -= 1
                neighbours.append(down)
            if coord[d] < (grid_dimensions[d] - 1):
                up = coord[:]
                up[d] += 1
                neighbours.append(up)

        for neighbour in neighbours:
            serial = _serialise_coordinate(neighbour, spacing)
            adjacency[i, serial] = 1

    return(adjacency)


def _get_minimum_distances(adjacency):
    """Finds the shortest path between pairs of graph nodes, given an adjacency
    graph.

    Based on the `Floyd-Warshall algorithm
    <https://en.wikipedia.org/wiki/Floyd-Warshall_algorithm>`_.
    """
    assert (len(adjacency.shape) == 2 and adjacency.shape[0] == adjacency.shape[1]), \
        "adjacency isn't a square matrix"
    n_nodes = adjacency.shape[0]
    distance = adjacency.copy()
    for k in range(n_nodes):
        for i in range(n_nodes):
            for j in range(n_nodes):
                if distance[i, j] > distance[i, k] + distance[k, j]:
                    distance[i, j] = distance[i, k] + distance[k, j]

    return(distance)


class SelfOrganizingMap(BaseEstimator):
    """Self-Organizing Map

    Parameters
    ----------
    adjacency : tuple of integers, or ndarray, default: (4, 4)
        Form of the SOM grid to use. If a tuple of integers is passed,
        a orthotopic grid topology will be generated with those dimensions.
        If an ndarray is passed, it should be an adjacency matrix of the
        SOM grid, of dimension [n_centers, n_centers].

    n_iterations : int
        Number of iterations of the SOM algorithm to run

    learning_rate : float
        Learning rate (alpha in Kohonen [1990])

    init : 'random' or ndarray
        Method for initialization, defaults to 'random':

        'random' : use randomly chosen cluster centers.

        ndarray : an array of initial cluster centers [n_centers, n_features].


    Attributes
    ----------
    centers_ : array, [(x, y), n_features]
        Coordinates of centers and value

    labels_ :
        Labels of each point

    Notes
    ------
    References :
    - Kohonen, T., 1990. The Self-Organizing Map. Proceedings of the IEEE,
      78(9), pp.1464-1480. doi://10.1109/5.58325
    - Kohonen, T., 2013. Essentials of the self-organizing map. Neural
      Networks, 37, pp.52-65. doi://10.1016/j.neunet.2012.09.018
    """

    def __init__(self, adjacency=(4, 4), init='random', n_iterations=64,
                 learning_rate=1, callback=None):
        if isinstance(adjacency, int):
            adjacency = (adjacency,)

        if isinstance(adjacency, tuple):
            n_centers = np.prod(adjacency)
            if isinstance(init, np.ndarray) and (n_centers != init.shape[0]):
                raise ValueError("'init' contains %d centers, but 'adjacency' specifies %d clusters"
                                 % (init.shape[0], np.prod(adjacency)))

            adjacency = _generate_adjacency_matrix(adjacency)

        self.adjacency_matrix = adjacency
        self.distance_matrix = _get_minimum_distances(self.adjacency_matrix)
        self.graph_diameter = self.distance_matrix.max()
        self.n_centers = n_centers
        self.init = init
        self.n_iterations = n_iterations
        self.learning_rate = learning_rate
        self.callback = callback

    def fit(self, X):
        """Perform Self Organising Map clustering from features.

        Given an sample of X, we randomly choose one of them for each
        iteration.
        A good ratio, nb X = 2 or 3 x nbiter

        Parameters
        ----------
        X : ndarray [n_samples, n_features]
            Sample data array.

        """
        assert isinstance(X, np.ndarray), 'X is not an array!'
        self.centers_ = None
        self.dim = X.shape[-1]

        # init centers_
        if self.init == 'random':
            self.centers_ = np.random.rand(self.n_centers, self.dim)
        elif isinstance(self.init, np.ndarray):
            assert self.init.shape[-1] == self.dim
            self.centers_ = self.init

        # iteration loop
        iteration = 0
        # This can have duplicates. Would it make more sense to use
        # np.random.permutation(X)?
        indices = np.random.random_integers(0, len(X)-1, self.n_iterations)
        # TODO: this *was* based on the length of a square grid, now it's the
        # maximum diameter of the SOM topology. Is this OK?
        # See Kohonen (2013, p56)
        l = self.n_iterations / float(self.graph_diameter)
        for i in indices:
            lr = self.learning_rate * np.exp(-iteration / l)
            self._learn_x(X[i], lr, iteration)
            iteration += 1
            if self.callback is not None:
                self.callback(self, iteration)

        # assign labels
        self.labels_ = [self.best_matching_center(x) for x in X]
        return self

    def _learn_x(self, x, lr, iteration):
        winner = self.best_matching_center(x)
        radius = self.radius_of_the_neighborhood(iteration)
        updatable = self.centers_in_radius(winner, radius)
        # See Kohonen (2013, p56)
        neighborhood = np.exp(-np.sum((self.centers_[winner] - self.centers_[updatable])**2, axis=1)/(2*radius**2))
        self.centers_[updatable] = self.centers_[updatable] + \
            lr * np.asmatrix(neighborhood).T * np.asmatrix(x - self.centers_[winner])

    def best_matching_center(self, x):
        assert x.shape == self.centers_[1].shape
        distances = np.sum((x - self.centers_)**2, axis=1)
        return(distances.argmin())

    def centers_in_radius(self, winner, radius):
        return(np.where(self.distance_matrix[winner] < radius))

    def radius_of_the_neighborhood(self, iteration):
        # TODO: see above. This should initially cover about half the grid.
        l = self.n_iterations / self.graph_diameter
        return self.n_centers * np.exp(-iteration / l)
