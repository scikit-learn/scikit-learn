"""
 Self-organizing map

 Reference : (to check)
 Kohonen, T.; , "The self-organizing map,"
 Proceedings of the IEEE , vol.78, no.9, pp.1464-1480, Sep 1990
"""
# Authors: Sebastien Campion <sebastien.campion@inria.fr>
# License: BSD
from __future__ import division
import math
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

    adjacency = np.zeros((n_centres, n_centres))

    spacing = [int(np.prod(grid_dimensions[(k+1):])) for k in range(len(grid_dimensions))]

    for i in range(n_centres):
        coord = _unserialise_coordinate(i, spacing)
        neighbours = []
        for d in range(len(grid_dimensions)):
            if coord[d] > 0:
                down = coord.copy()
                down[d] -= 1
                neighbours.append(down)
            if coord[d] < (grid_dimensions[d] - 1):
                up = coord.copy()
                up[d] += 1
                neighbours.append(up)

        for neighbour in neighbours:
            serial = _serialise_coordinate(neighbour, spacing)
            adjacency[i,serial] = 1

    return(adjacency)

class SelfOrganizingMap(BaseEstimator):
    """Self-Organizing Map

    Parameters
    ----------
    X : ndarray
        A M by N array of M observations in N dimensions or a length
        M array of N one-dimensional observations.

    size : int
        Width and height of the square map as well as the number of
        centroids to generate. If init initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.

    n_iterations : int
        Number of iterations of the som algrithm to run

    learning_rate : float
        Learning rate

    init : {'random', 'matrix'}
        Method for initialization, defaults to 'random':

        'random': randomly points choosed

        'matrix': interpret the size parameter as a size by M array
         of initial neurons.

    Methods
    -------
    fit(X):
        Compute SOM

    Attributes
    ----------
    neurons_: array, [(x,y), n_features]
        Coordinates of neurons and value

    labels_:
        Labels of each point

    Notes
    ------

    """

    def __init__(self, size=16, init='random', n_iterations=64,
                 learning_rate=1, callback=None):
        self.size = size
        self.init = init
        self.n_iterations = n_iterations
        self.learning_rate = learning_rate
        self.callback = callback

    def fit(self, X, **params):
        """Given an sample of X, we randomly choose one of them for each
        iteration.
        A good ratio, nb X = 2 or 3 x nbiter"""
        X = np.asanyarray(X)
        self._set_params(**params)
        self.neurons_ = None
        dim = X.shape[-1]

        # init neurons_
        if self.init == 'random':
            self.neurons_ = np.random.rand(self.size, self.size, dim)
        elif self.init == 'matrix':
            assert len(self.size.shape) == 3
            self.neurons_ = self.size
            self.size = self.neurons_.shape[0]

        # iteration loop
        iteration = 0 
        indices = np.random.random_integers(0, len(X)-1, self.n_iterations)
        for i in indices:
            l = self.n_iterations / self.size
            lr = self.learning_rate * math.exp(-iteration / l)
            self._learn_vector(X[i], lr, iteration)
            iteration += 1
            if self.callback != None:
                self.callback(self, iteration)

        # assign labels
        self.labels_ = [self.bmu(x) for x in X]
        return self

    def _learn_vector(self, vector, lr, iteration):
        winner = self.bmu(vector)
        radius = self.radius_of_the_neighbordhood(iteration)
        for n in self.neurons_in_radius(winner, radius):
            nx, ny = n
            wt = self.neurons_[nx][ny]
            dr = self.dist(winner, n, radius)
            self.neurons_[nx][ny] = wt + dr * lr * (vector - wt)

    def bmu(self, vector):
        """Best matching unit
        """
        assert vector.shape[0] == self.neurons_.shape[-1]
        vector = np.resize(vector, self.neurons_.shape)
        dists = np.sum((vector - self.neurons_)**2, axis=-1)
        min = dists.argmin()
        #w = np.unravel_index(min,dists.shape)
        return divmod(min, self.size)

    def dist(self, w, n, radius):
        wx, wy = w
        nx, ny = n
        d = (wx - nx)**2 + (wy - ny)**2
        # offcial paper implementation : return math.exp(-d/2*radius**2)
        return math.exp(-d / radius)

    def neurons_in_radius(self, winner, radius):
        wi, wj = winner
        x = y = np.arange(self.size)
        xx, yy = np.meshgrid(x, y)
        v = np.sqrt((xx - wi)**2 + (yy - wj)**2) < radius
        return np.c_[np.nonzero(v)]

    def radius_of_the_neighbordhood(self, iteration):
        l = self.n_iterations / self.size
        return self.size * math.exp(-iteration / l)
