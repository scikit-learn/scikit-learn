"""
 Self-organizing map
"""

# Authors: Sebastien Campion <sebastien.campion@inria.fr>
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
    """Generate an adjacency matrix for nodes of an orthotopic grid with
    dimensions given by grid_dimensions
    """
    n_centres = np.prod(grid_dimensions)

    # initialise to 0 for self-distances, infinity non-connections
    adjacency = np.empty((n_centres, n_centres), dtype=float)
    adjacency.fill(np.inf)
    np.fill_diagonal(adjacency, 0)

    spacing = [int(np.prod(grid_dimensions[(k+1):])) for k in range(len(grid_dimensions))]

    for i in range(n_centres):
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
            adjacency[i,serial] = 1

    return(adjacency)

def _get_minimum_distances(adjacency):
    """Finds the shortest path between pairs of graph nodes, given an adjacency
    graph.

    Based on the [Floyd-Warshall algorithm](https://en.wikipedia.org/wiki/Floyd-Warshall_algorithm)
    """
    assert (len(adjacency.shape) == 2 and adjacency.shape[0] == adjacency.shape[1]), "adjacency isn't a square matrix"
    n_nodes = adjacency.shape[0]
    distance = adjacency.copy()
    for k in range(n_nodes):
        for i in range(n_nodes):
            for j in range(n_nodes):
                if distance[i,j] > distance[i,k] + distance[k,j]:
                    distance[i,j] = distance[i,k] + distance[k,j]

    return(distance)

class SelfOrganizingMap(BaseEstimator):
    """Self-Organizing Map

    Parameters
    ----------
    size : int
        Width and height of the square map as well as the number of
        centroids to generate. If init initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.

    n_iterations : int
        Number of iterations of the SOM algorithm to run

    learning_rate : float
        Learning rate (alpha in Kohonen [1990])

    init : {'random', 'matrix'}
        Method for initialization, defaults to 'random':

        'random' : use randomly chosen cluster centres.

        'matrix': interpret the size parameter as a size by M array
         of initial centres.


    Attributes
    ----------
    centres_ : array, [(x,y), n_features]
        Coordinates of centres and value

    labels_ :
        Labels of each point

    Notes
    ------
    Reference :
    Kohonen, T.; , "The self-organizing map,"
    Proceedings of the IEEE , vol.78, no.9, pp.1464-1480, Sep 1990
    """

    def __init__(self, size=16, init='random', n_iterations=64,
                 learning_rate=1, callback=None):
        self.size = size
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
        X = np.asanyarray(X)
        self.centres_ = None
        self.dim = X.shape[-1]

        # init centres_
        if self.init == 'random':
            self.centres_ = np.random.rand(self.size, self.size, self.dim)
        elif self.init == 'matrix':
            assert len(self.size.shape) == 3
            self.centres_ = self.size
            self.size = self.centres_.shape[0]

        # iteration loop
        iteration = 0
        indices = np.random.random_integers(0, len(X)-1, self.n_iterations)
        l = self.n_iterations / self.size
        for i in indices:
            lr = self.learning_rate * np.exp(-iteration / l)
            self._learn_x(X[i], lr, iteration)
            iteration += 1
            if self.callback != None:
                self.callback(self, iteration)

        # assign labels
        self.labels_ = [self.best_matching_centre(x) for x in X]
        return self

    def _learn_x(self, x, lr, iteration):
        winner = self.best_matching_centre(x)
        radius = self.radius_of_the_neighborhood(iteration)
        for n in self.centres_in_radius(winner, radius):
            nx, ny = n
            wt = self.centres_[nx][ny]
            dr = self.dist(winner, n, radius)
            self.centres_[nx][ny] = wt + dr * lr * (x - wt)

    def best_matching_centre(self, x):
        assert x.shape[0] == self.centres_.shape[-1]
        x = np.resize(x, self.centres_.shape)
        dists = np.sum((x - self.centres_)**2, axis=-1)
        min = dists.argmin()
        #w = np.unravel_index(min,dists.shape)
        return divmod(min, self.size)

    def dist(self, w, n, radius):
        wx, wy = w
        nx, ny = n
        d = (wx - nx)**2 + (wy - ny)**2
        # Official paper implementation : return np.exp(-d/2*radius**2)
        return np.exp(-d / radius)

    def centres_in_radius(self, winner, radius):
        wi, wj = winner
        x = y = np.arange(self.size)
        xx, yy = np.meshgrid(x, y)
        v = np.sqrt((xx - wi)**2 + (yy - wj)**2) < radius
        return np.c_[np.nonzero(v)]

    def radius_of_the_neighborhood(self, iteration):
        l = self.n_iterations / self.size
        return self.size * np.exp(-iteration / l)
