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
        coord.append(int(np.floor(serial / spacing[i])))
        serial = int(serial % spacing[i])

    return(coord)


def _serialise_coordinate(coord, spacing):
    return(np.dot(coord, spacing))


def _generate_adjacency_matrix(grid_dimensions):
    """Generate an adjacency matrix of a rectangular grid.

    grid_dimensions : tuple of integers
    Specifies the dimensions of the grid for which to generate an adjacency
    matrix.
    """
    n_centers = np.prod(grid_dimensions)

    # initialise to 0 for self-distances, infinity non-connections
    adjacency = np.empty((n_centers, n_centers), dtype=float)
    adjacency.fill(np.inf)
    np.fill_diagonal(adjacency, 0)

    spacing = [int(np.prod(grid_dimensions[(k + 1):]))
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
    """Generate an all-pairs shortest path matrix from an adjacency matrix.

    adjacency : symmetric square ndarray
    Matrix with distances between pairs of graph nodes (usually 1s), or np.inf
    for no connection.

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
        Learning rate (initial alpha in Kohonen [1990])

    init : 'random' or ndarray
        Method for initialization, defaults to 'random':

        'random' : use randomly chosen cluster centers.

        ndarray : an array of initial cluster centers [n_centers, n_features].


    Attributes
    ----------
    cluster_centers_ : array, [(x, y), n_features]
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

        self.adjacency = adjacency
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
        # If adjacency is an int or tuple, we generate an adjacency matrix
        # of the specified dimensions. Otherwise if it's an adjacency matrix,
        # use that.
        if isinstance(self.adjacency, int):
            self.adjacency = (self.adjacency,)

        if isinstance(self.adjacency, tuple):
            n_centers = np.prod(self.adjacency)
            adjacency_matrix = _generate_adjacency_matrix(self.adjacency)
        else:
            assert isinstance(self.adjacency, np.ndarray), "'adjacency' is not an int, tuple or array!"
            adjacency_matrix = self.adjacency
            n_centers = adjacency_matrix.shape[0]

        if isinstance(self.init, np.ndarray) and (n_centers != self.init.shape[0]):
            raise ValueError("'init' contains %d centers, but 'adjacency' specifies %d clusters"
                             % (self.init.shape[0], n_centers))

        self.adjacency_matrix_ = adjacency_matrix
        self.distance_matrix_ = _get_minimum_distances(self.adjacency_matrix_)
        self.graph_diameter_ = self.distance_matrix_.max()
        self.n_centers_ = n_centers

        assert isinstance(X, np.ndarray), 'X is not an array!'
        self.cluster_centers_ = None
        self.dim_ = X.shape[-1]

        # init cluster_centers_
        if self.init == 'random':
            self.cluster_centers_ = np.random.rand(self.n_centers_, self.dim_)
        elif isinstance(self.init, np.ndarray):
            assert self.init.shape[-1] == self.dim_
            self.cluster_centers_ = self.init.copy()

        # iteration loop
        iteration = 0
        # This can have duplicates. Would it make more sense to use
        # np.random.permutation(X)?
        indices = np.random.random_integers(0, len(X) - 1, self.n_iterations)
        for i in indices:
            self._learn_x(X[i], iteration)
            iteration += 1
            if self.callback is not None:
                self.callback(self, iteration)

        # assign labels
        self.labels_ = np.array([self.best_matching_center(x) for x in X])
        return self

    def _learn_x(self, x, iteration):
        # alpha to suit Kohonen (2013, p56)
        # TODO: this *was* based on the length of a square grid, now it's the
        # maximum diameter of the SOM topology. Is this OK?
        l = self.n_iterations / float(self.graph_diameter_)
        self.alpha_ = self.learning_rate * np.exp(-iteration / l)
        winner = self.best_matching_center(x)
        radius = self.radius_of_the_neighborhood(iteration)
        updatable = self.cluster_centers_in_radius(winner, radius)
        distances = self.distance_matrix_[winner][updatable]
        # neighborhood function from Kohonen (2013, p56)
        neighborhood = self.alpha_ * np.exp(-distances / (2 * radius ** 2))
        self.cluster_centers_[updatable] = self.cluster_centers_[updatable] + \
            np.multiply(neighborhood, (x - self.cluster_centers_[updatable]).T).T

    def best_matching_center(self, x):
        assert x.shape == self.cluster_centers_[1].shape
        distances = np.sum((x - self.cluster_centers_) ** 2, axis=1)
        return(distances.argmin())

    def cluster_centers_in_radius(self, winner, radius):
        return(np.where(self.distance_matrix_[winner] < radius)[0])

    def radius_of_the_neighborhood(self, iteration):
        """Kohonen's sigma.

        Monotonically decreasing function of iteration. Kohonen (1990) doesn't
        specify a particular form.
        """
        # TODO: see above. This should initially cover about half the grid.
        l = self.n_iterations / self.graph_diameter_
        return self.n_centers_ * np.exp(-iteration / l)
