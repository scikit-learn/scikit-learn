"""Self-Organizing Map"""

# Authors: Julio Faracco <jcfaracco@gmail.com>
# License: BSD 3 clause

import time
import warnings

import numpy as np

from ..base import TransformerMixin
from ..base import ClusterMixin
from ..base import BaseEstimator
from ..utils import check_random_state
from ..utils.validation import check_is_fitted
from ..utils.validation import _deprecate_positional_args
from ..metrics.pairwise import euclidean_distances

def asymptotic_decay(learning_rate, t, n_iter):
    """Decay function of the learning process

    Parameters
    ----------
    learning_rate : float
        Current learning rate.
    t : int
        Current iteration.
    n_iter : int
        Maximum number of iterations for the training.

    Returns
    -------
    learning_rate: float
        New learning rate updated based on iteration
    """
    return learning_rate * (1.0 - (t / n_iter))


class SOM(TransformerMixin, ClusterMixin, BaseEstimator):
    """Kohonen Self-Organizing Maps (SOM).

    Read more in the :ref:`User Guide <som>`.

    Parameters
    ----------

    width : int, default=3
        Width of Self-Organizing Map network topology.
        Obs: the default=3 is used to pass all scikit test cases.

    height : int, default=1
        Height of Self-Organizing Map network topology.
        Obs: the default=1 is used to pass all scikit test cases. With the
        width=3, SOM have a linear topology of (1, 3) shape.

    sigma : float, default=1.0
        Sigma value for neighborhood distance calculation.

    learning_rate : flat, default=0.5
        Learning rate training parameter adjusted all iterations.

    max_iter : int, default=3
        Maximum number of iterations of the SOM algorithm for a
        single run. It consider a N iterations repeatedly for a training
        sample.

    decay_function : func, default=asymptotic_decay()
        The learning rate decay of each interaction.

    neighborhood_function : {"gaussian", "mexican_hat", "bubble", "triangle"}, default="gaussian"
        Function which determines the rate of change of the neighborhood
        around the BMUs (winner neuron).

    topology : {"hexagonal", "rectangular"}, default="rectangular"
        SOM weights matrix topology. It can be rectangular and hexagonal.
        Other topologies are not implemented.

    activation_distance : {"euclidean", "cosine", "manhattan", "chebyshev"}, default="euclidean"
        Function that calculates the distance activation of a winner neuron.

    verbose : int, default=0
        Verbosity mode.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    Attributes
    ----------
    _activation_map : ndarray of shape (width, height)
        A map corresponding the activation neurons for a respective input.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point

    n_samples_ : int
        Number of samples to fit

    n_iter_ : int
        Number of iterations run.

    Examples
    --------

    >>> from sklearn.neural_network import SOM
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0]])
    >>> som = SOM(width=1, height=2, random_state=0).fit(X)
    >>> som.labels_
    array([1, 1, 1, 0, 0, 0])
    >>> som.predict([[0, 0], [12, 3]])
    array([1, 0])

    References
    ----------

    [1] Kohonen, T., "Self-Organizing Map", 2nd ed., Springer-Verlag,
        Berlin, 1995, pp. 113.

    [2] Kiviluoto, K., "Topology Preservation in Self-Organizing Maps",
        in the proceeding of International Conference on Neural
        Networks (ICNN), 1996, pp. 294-299.

    [3] Vettigli, Giuseppe. "MiniSom: minimalistic and NumPy-based
        implementation of the Self Organizing Map." (2013).
        https://github.com/JustGlowing/minisom/

    [4] Comitani, Federico. "SimpSOM" (2019).
        https://github.com/fcomitani/SimpSOM/
    """
    @_deprecate_positional_args
    def __init__(self, width=1, height=3, sigma=1.0, learning_rate=0.5,
                 max_iter=5, decay_function=asymptotic_decay,
                 neighborhood_function='gaussian', topology='rectangular',
                 activation_distance='euclidean', verbose=0, random_state=0):
        self.height, self.width = height, width
        self.sigma = sigma
        self.learning_rate = learning_rate
        self.max_iter = max_iter
        self.decay_function = decay_function
        self.topology = topology
        self.verbose = verbose
        self.random_state = random_state
        self.neighborhood_function = neighborhood_function
        self.activation_distance = activation_distance

    def _check_params(self):
        self.__neighborhood_function = ['gaussian', 'mexican_hat',
                                        'bubble', 'triangle']

        self.__topology = ['hexagonal', 'rectangular']

        self.__activation_distance = ['euclidean', 'cosine',
                                      'manhattan', 'chebyshev']

        if self.neighborhood_function not in self.__neighborhood_function:
            msg = '%s not supported only as neighborhood function'
            raise ValueError(msg % neighborhood_function)

        self._neighborhood_function = getattr(self, "_" + self.neighborhood_function)

        if self.activation_distance not in self.__activation_distance:
            msg = '%s not supported only as activation distance method'
            raise ValueError(msg % activation_distance)

        self._activation_distance = getattr(self, "_" + self.activation_distance)

        if self.topology not in self.__topology:
            msg = '%s not supported only hexagonal and rectangular available'
            raise ValueError(msg % self.topology)

        if self.neighborhood_function in ['triangle', 'bubble'] and \
           (np.divmod(self.sigma, 1)[1] != 0 or self.sigma < 1):
            warnings.warn('triangle or bubble requires a sigma >= 1')

        self._random_state = check_random_state(self.random_state)

        self._labels_ = None

    def _setup(self):
        """Setup components of SOM used in addition with weights
        """
        self._activation_map = np.zeros((self.width, self.height))

        self._weigths = np.zeros((self.width, self.height))

        self._n_x = np.arange(self.width)
        self._n_y = np.arange(self.height)

        self._xx, self._yy = np.meshgrid(self._n_x, self._n_y)
        self._xx = self._xx.astype(float)
        self._yy = self._yy.astype(float)

    def _randomize_weights(self, features):
        """Initializes the weights of the SOM with random values.
        """
        self._check_params()

        self._setup()

        r = check_random_state(self.random_state)

        self._weights = r.rand(self.width, self.height,
                               features)

    def _randomize_weights_from_data(self, X):
        """Initializes the weights of the SOM picking random samples from data.

        Parameters
        ----------
        X : np.array
            Data sample to randomize SOM weights array
        """
        X = self._check_test_data(X)

        self._check_params()

        self._setup()

        self._weights = np.zeros((self.width, self.height, len(X[-1])))

        it = np.nditer(self._activation_map, flags=['multi_index'])
        while not it.finished:
            rand_i = self._random_state.randint(len(X))
            self._weights[it.multi_index] = X[rand_i]
            it.iternext()

    def _distance_weights(self, X):
        """Returns a matrix d where d[i,j] is the euclidean distance between
        X[i] and the j-th weight.
        """
        input_data = np.array(X)
        weights_flat = self._weights.reshape(-1, self._weights.shape[2])
        input_data_sq = np.power(input_data, 2).sum(axis=1, keepdims=True)
        weights_flat_sq = np.power(weights_flat, 2).sum(axis=1, keepdims=True)
        cross_term = np.dot(input_data, weights_flat.T)
        dist = -2 * cross_term + input_data_sq + weights_flat_sq.T
        dist[dist < 0] = 0
        return np.sqrt(dist)

    def quantization(self, X):
        """Assigns a code book (weights vector of the winning neuron)
        to each sample in data.

        Parameters
        ----------
        X : np.array
            Data sample

        Returns
        -------
        weights : np.array
            Weight values of the winner neurons selected using BMU
        """
        winners_coords = self._get_bmus(X)

        return self._weights[np.unravel_index(winners_coords,
                                              self._weights.shape[:2])]

    def _quantization_error(self, X):
        """Returns the quantization error computed as the average
        distance between each input sample and its best matching unit.
        """
        q = self.quantization(X)

        return self._euclidean(X, q).mean()


    def _topological_error(self, X):
        """Returns the topographic error computed by finding
        the best-matching and second-best-matching neuron in the map
        for each input and then evaluating the positions.
        A sample for which these two nodes are not adjacent counts as
        an error. The topographic error is given by the
        the total number of errors divided by the total of samples.
        If the topographic error is 0, no error occurred.
        If 1, the topology was not preserved for any of the samples.

        Parameters
        ----------
        X : np.array
            Data sample to calculate topological error

        Returns
        -------
        mean : float
            Mean of all second BMUs that is not a neighbour of the BMU
            itself
        """
        X = self._check_test_data(X)

        if np.prod(self._activation_map.shape) == 1:
            warn('The topographic error is not defined for a 1-by-1 map.')
            return np.nan

        b2mu_idxs = np.argsort(self._distance_weights(X), axis=1)[:, :2]
        b2my_xy = np.unravel_index(b2mu_idxs, self._weights.shape[:2])
        if self.topology == 'hexagonal':
            new_b2my_xy = list()
            for array in np.array(b2my_xy):
                new_b2my_xy.append(np.apply_along_axis(self._transform_weights_to_hex, 1, array))
            b2my_xy = tuple(new_b2my_xy)
        b2mu_x, b2mu_y = b2my_xy[0], b2my_xy[1]
        dxdy = np.hstack([np.diff(b2mu_x), np.diff(b2mu_y)])
        distance = np.linalg.norm(dxdy, axis=1)
        return (distance >= np.sqrt(2)).mean()

    def _get_bmus(self, X):
        X = self._check_test_data(X)

        return np.argmin(self._distance_weights(X), axis=1)

    def get_coordinates(self):
        """Returns the position of the neurons on an euclidean
        plane that reflects the chosen topology in two meshgrids xx and yy.
        Neuron with map coordinates (1, 4) has coordinate (xx[1, 4], yy[1, 4])
        in the euclidean plane.
        Only useful if the topology chosen is not rectangular.
        """
        if self.topology == "hexagonal":
            if len(self._xx.T) != len(self._yy.T):
                raise ValueError("Meshgrid xx and yy have different lengths.")

            for i in range(0, len(self._xx.T)):
                for j in range(0, len(self._xx.T[i])):
                    coord = [self._xx.T[i][j], self._yy.T[i][j]]
                    hex_coord = self._transform_weights_to_hex(coord)
                    self._xx.T[i][j] = hex_coord[0]
                    self._yy.T[i][j] = hex_coord[1]
        return self._xx.T, self._yy.T

    def distance_map(self):
        """Returns the distance map of the weights.
        Each cell is the normalised sum of the distances between
        a neuron and its neighbours. Note that this method uses
        the euclidean distance."""
        um = np.zeros((self._weights.shape[0],
                       self._weights.shape[1],
                       8))  # 2 spots more for hexagonal topology

        if self.topology == 'hexagonal':
            ii = [[1, 1, 1, 0, -1, 0], [0, 1, 0, -1, -1, -1]]
            jj = [[1, 0, -1, -1, 0, 1], [1, 0, -1, -1, 0, 1]]
        else:
            ii = [[0, -1, -1, -1, 0, 1, 1, 1]]*2
            jj = [[-1, -1, 0, 1, 1, 1, 0, -1]]*2

        for x in range(self._weights.shape[0]):
            for y in range(self._weights.shape[1]):
                w_2 = self._weights[x, y]
                e = y % 2 == 0   # only used on hexagonal topology
                for k, (i, j) in enumerate(zip(ii[e], jj[e])):
                    if (x+i >= 0 and x+i < self._weights.shape[0] and
                            y+j >= 0 and y+j < self._weights.shape[1]):
                        w_1 = self._weights[x+i, y+j]
                        w_d = w_2 - w_1
                        um[x, y, k] = np.sqrt(np.dot(w_d, w_d.T))

        um = um.sum(axis=2)
        return um/um.max()

    def quality(self, X):
        """Calculate the quality of the SOM network based on topology
        and neuron activation versus neighborhood.

        Parameters
        ----------
        X : np.array
            Data sample to calculate the quantization error and 
            topology error

        Returns
        -------
        qe, te : tuple
            A tuple with both error rates
        """
        qe = self._quantization_error(X)
        te = self._topological_error(X)
        return qe, te

    def set_weights(self, weights):
        """Set initial weights to SOM network

        Parameters
        ----------
        weights : np.array
            A predefined weights (maybe based on previous training).
        """
        self._weights = weights

    def get_weights(self):
        """Get initial weights to SOM network

        Returns
        -------
        weights : np.array
            The weights of SOM network. Even if it is defined, trained or
            set only.
        """
        if not hasattr(self, "_weights"):
            return np.empty((self.width, self.height))

        return self._weights

    def _transform_weights_to_hex(self, point):
        """Convert Cartesian coordinates to hexagonal tiling coordinates.

        Parameters
        ----------
        x : float
            Position along the x-axis of Cartesian coordinates.

        y : float
            Position along the y-axis of Cartesian coordinates.

        Returns
        -------
        array : np.array
            A 2d array containing the coordinates in the new space based
            on SOM topology.
        """
        if len(point) != 2:
            return np.nan

        x, y = point[0], point[1]

        newy = y * 2 / np.sqrt(3) * 3 / 4
        newx = x
        if y % 2:
            newx += 0.5
        return [newx, newy]

    ##########################################
    # Section of activation distance functions
    ##########################################
    def _cosine(self, x, w):
        num = (w * x).sum(axis=2)
        denum = np.multiply(np.linalg.norm(w, axis=2), np.linalg.norm(x))
        return 1 - num / (denum + 1e-8)

    def _euclidean(self, x, w):
        return np.linalg.norm(np.subtract(x, w), axis=-1)

    def _manhattan(self, x, w):
        return np.linalg.norm(np.subtract(x, w), ord=1, axis=-1)

    def _chebyshev(self, x, w):
        return np.max(np.subtract(x, w), axis=-1)

    ###################################
    # Section of neighborhood functions
    ###################################
    def _gaussian(self, c, sigma):
        """Returns a Gaussian centered in c."""
        d = 2 * (sigma ** 2)
        ax = np.exp(-np.power(self._xx-self._xx.T[c], 2)/d)
        ay = np.exp(-np.power(self._yy-self._yy.T[c], 2)/d)
        return (ax * ay).T  # the external product gives a matrix

    def _mexican_hat(self, c, sigma):
        """Mexican hat centered in c."""
        p = np.power(self._xx-self._xx.T[c], 2) + np.power(self._yy-self._yy.T[c], 2)
        d = 2*sigma*sigma
        return (np.exp(-p/d)*(1-2/d*p)).T

    def _bubble(self, c, sigma):
        """Constant function centered in c with spread sigma.
        sigma should be an odd value.
        """
        ax = np.logical_and(self._n_x > c[0] - sigma,
                            self._n_x < c[0] + sigma)
        ay = np.logical_and(self._n_y > c[1] - sigma,
                            self._n_y < c[1] + sigma)
        return np.outer(ax, ay) * 1.

    def _triangle(self, c, sigma):
        """Triangular function centered in c with spread sigma."""
        triangle_x = (-np.abs(c[0] - self._n_x)) + sigma
        triangle_y = (-np.abs(c[1] - self._n_y)) + sigma
        triangle_x[triangle_x < 0] = 0.
        triangle_y[triangle_y < 0] = 0.
        return np.outer(triangle_x, triangle_y)

    def _build_iteration_indexes(self, input_len, num_iterations, random=False):
        """Returns an iterable with the indexes of the samples
        to pick at each iteration of the training.
        If random_generator is not None, it must be an instalce
        of numpy.random.RandomState and it will be used
        to randomize the order of the samples."""
        iterations = np.arange(num_iterations) % input_len
        if random:
            self._random_state.shuffle(iterations)
        return iterations

    def _check_test_data(self, X):
        X = self._validate_data(X, accept_sparse='csr', reset=False,
                                dtype=[np.float64, np.float32],
                                order='C', accept_large_sparse=False)
        return X

    def _activate(self, x):
        """Updates matrix activation_map, in this matrix
           the element i,j is the response of the neuron i,j to x."""
        self._activation_map = self._activation_distance(x, self._weights)

    def _winner(self, x):
        """Computes the coordinates of the winning neuron for the sample x."""
        self._activate(x)
        return np.unravel_index(self._activation_map.argmin(),
                                self._activation_map.shape)

    def _update_weights(self, x, win, t, max_iteration):
        """Updates the weights of the neurons.

        Parameters
        ----------
        x : np.array
            Current pattern to learn.

        win : tuple
            Position of the winning neuron for x (array or tuple).

        t : int
            Iteration index

        max_iteration : int
            Maximum number of training itarations.
        """
        eta = self.decay_function(self.learning_rate, t, max_iteration)
        # sigma and learning rate decrease with the same rule
        sig = self.decay_function(self.sigma, t, max_iteration)
        # improves the performances
        g = self._neighborhood_function(win, sig) * eta
        # w_new = eta * neighborhood_function * (x-w)
        self._weights += np.einsum('ij, ijk->ijk', g, x - self._weights)

    def weights(self):
        return self._weights

    def fit(self, X, y=None, random_data=False, random_order=False):
        """Train SOM.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training instances to cluster. It must be noted that the data
            will be converted to C ordering, which will cause a memory
            copy if the given data is not C-contiguous.
            If a sparse matrix is passed, a copy will be made if it's not in
            CSR format.

        y : Ignored
            Not used, present here for API consistency by convention.

        random_data : bool, default=False
            Use random data to initialize weights

        random_order : bool, default=False
            Randomize input data over iteration to get a better random weights
            during training.

        Returns
        -------
        self
            Fitted estimator.
        """
        if self.max_iter < 1:
            raise ValueError("Number of iterations needs to be higher than 1")

        X = self._check_test_data(X)

        self.n_samples_ = X.shape[0]
        self.n_features_in_ = X.shape[1]

        # SOM weights are not initialized previously
        if random_data:
            self._randomize_weights_from_data(X)
        else:
            self._randomize_weights(self.n_features_in_)

        self._check_params()

        self.n_iter_ = len(X) * self.max_iter

        random_generator = None
        iterations = self._build_iteration_indexes(len(X), self.n_iter_,
                                                   random_order)

        verbose = self.verbose
        begin = time.time()
        total = 0
        for t, iteration in enumerate(iterations):
            self._update_weights(X[iteration], self._winner(X[iteration]),
                                 t, self.n_iter_)

            if verbose:
                end = time.time()
                qe, te = self.quality(X)
                print("[%s] Epoch %d Iteration %d, mean quantization error = %.2f,"
                      " topographic error = %.2f, time = %.2fs"
                      % (type(self).__name__, epoch, t, qe, te,
                          end - begin))
                total += end - begin
                begin = end

        self._labels_ = self._get_bmus(X)

        if verbose:
            print("Time total to fit: %.2fs" % total)

        return self

    def fit_transform(self, X, y=None, random_order=False):
        """Train SOM transform X to cluster-distance space.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training instances to cluster. It must be noted that the data
            will be converted to C ordering, which will cause a memory
            copy if the given data is not C-contiguous.
            If a sparse matrix is passed, a copy will be made if it's not in
            CSR format.

        y : Ignored
            Not used, present here for API consistency by convention.

        random_order : bool, default=False
            Randomize input data over iteration to get a better random weights
            during training.

        Returns
        -------
        self
            Fitted estimator.
        """
        return self.fit(X, random_order=random_order)._transform(X)

    @property
    def labels_(self):
        check_is_fitted(self)

        return self._labels_

    def transform(self, X):
        """Transform X to a cluster-distance space.

        In the new space, each dimension is the distance to the cluster
        centers. Note that even if X is sparse, the array returned by
        `transform` will typically be dense.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to transform.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_clusters)
            X transformed in the new space.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        return self._transform(X)

    def _transform(self, X):
        """Guts of transform method; no input validation."""
        return self._euclidean(X, self.quantization(X))

    def predict(self, X, sample_weight=None):
        """Predict the closest cluster each sample in X belongs to.

        In the vector quantization literature, `weights` is called
        the code book and each value returned by `predict` is the index of
        the closest code in the code book.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to predict.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        check_is_fitted(self)

        return self._get_bmus(X)
