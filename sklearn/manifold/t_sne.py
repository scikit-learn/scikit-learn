# Author: Alexander Fabisch  -- <afabisch@informatik.uni-bremen.de>
# Author: Christopher Moody <chrisemoody@gmail.com>
# Author: Nick Travers <nickt@squareup.com>
# License: BSD 3 clause (C) 2014

# This is the exact and Barnes-Hut t-SNE implementation. There are other
# modifications of the algorithm:
# * Fast Optimization for t-SNE:
#   http://cseweb.ucsd.edu/~lvdmaaten/workshops/nips2010/papers/vandermaaten.pdf

import numpy as np
from scipy import linalg
import scipy.sparse as sp
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from ..neighbors import BallTree
from ..base import BaseEstimator
from ..utils import check_array
from ..utils import check_random_state
from ..utils.extmath import _ravel
from ..decomposition import PCA
from ..metrics.pairwise import pairwise_distances
from . import _utils
from . import _barnes_hut_tsne
from ..utils.fixes import astype
from ..externals.six import string_types


MACHINE_EPSILON = np.finfo(np.double).eps


def _joint_probabilities(distances, desired_perplexity, verbose):
    """Compute joint probabilities p_ij from distances.

    Parameters
    ----------
    distances : array, shape (n_samples * (n_samples-1) / 2,)
        Distances of samples are stored as condensed matrices, i.e.
        we omit the diagonal and duplicate entries and store everything
        in a one-dimensional array.

    desired_perplexity : float
        Desired perplexity of the joint probability distributions.

    verbose : int
        Verbosity level.

    Returns
    -------
    P : array, shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.
    """
    # Compute conditional probabilities such that they approximately match
    # the desired perplexity
    distances = astype(distances, np.float32, copy=False)
    conditional_P = _utils._binary_search_perplexity(
        distances, None, desired_perplexity, verbose)
    P = conditional_P + conditional_P.T
    sum_P = np.maximum(np.sum(P), MACHINE_EPSILON)
    P = np.maximum(squareform(P) / sum_P, MACHINE_EPSILON)
    return P


def _joint_probabilities_nn(distances, neighbors, desired_perplexity, verbose):
    """Compute joint probabilities p_ij from distances using just nearest
    neighbors.

    This method is approximately equal to _joint_probabilities. The latter
    is O(N), but limiting the joint probability to nearest neighbors improves
    this substantially to O(uN).

    Parameters
    ----------
    distances : array, shape (n_samples * (n_samples-1) / 2,)
        Distances of samples are stored as condensed matrices, i.e.
        we omit the diagonal and duplicate entries and store everything
        in a one-dimensional array.

    desired_perplexity : float
        Desired perplexity of the joint probability distributions.

    verbose : int
        Verbosity level.

    Returns
    -------
    P : array, shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.
    """
    # Compute conditional probabilities such that they approximately match
    # the desired perplexity
    distances = astype(distances, np.float32, copy=False)
    neighbors = astype(neighbors, np.int64, copy=False)
    conditional_P = _utils._binary_search_perplexity(
        distances, neighbors, desired_perplexity, verbose)
    m = "All probabilities should be finite"
    assert np.all(np.isfinite(conditional_P)), m
    P = conditional_P + conditional_P.T
    sum_P = np.maximum(np.sum(P), MACHINE_EPSILON)
    P = np.maximum(squareform(P) / sum_P, MACHINE_EPSILON)
    assert np.all(np.abs(P) <= 1.0)
    return P


def _kl_divergence(params, P, degrees_of_freedom, n_samples, n_components,
                   skip_num_points=0):
    """t-SNE objective function: gradient of the KL divergence
    of p_ijs and q_ijs and the absolute error.

    Parameters
    ----------
    params : array, shape (n_params,)
        Unraveled embedding.

    P : array, shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.

    degrees_of_freedom : float
        Degrees of freedom of the Student's-t distribution.

    n_samples : int
        Number of samples.

    n_components : int
        Dimension of the embedded space.

    skip_num_points : int (optional, default:0)
        This does not compute the gradient for points with indices below
        `skip_num_points`. This is useful when computing transforms of new
        data where you'd like to keep the old data fixed.

    Returns
    -------
    kl_divergence : float
        Kullback-Leibler divergence of p_ij and q_ij.

    grad : array, shape (n_params,)
        Unraveled gradient of the Kullback-Leibler divergence with respect to
        the embedding.
    """
    X_embedded = params.reshape(n_samples, n_components)

    # Q is a heavy-tailed distribution: Student's t-distribution
    n = pdist(X_embedded, "sqeuclidean")
    n += 1.
    n /= degrees_of_freedom
    n **= (degrees_of_freedom + 1.0) / -2.0
    Q = np.maximum(n / (2.0 * np.sum(n)), MACHINE_EPSILON)

    # Optimization trick below: np.dot(x, y) is faster than
    # np.sum(x * y) because it calls BLAS

    # Objective: C (Kullback-Leibler divergence of P and Q)
    kl_divergence = 2.0 * np.dot(P, np.log(P / Q))

    # Gradient: dC/dY
    grad = np.ndarray((n_samples, n_components))
    PQd = squareform((P - Q) * n)
    for i in range(skip_num_points, n_samples):
        np.dot(_ravel(PQd[i]), X_embedded[i] - X_embedded, out=grad[i])
    grad = grad.ravel()
    c = 2.0 * (degrees_of_freedom + 1.0) / degrees_of_freedom
    grad *= c

    return kl_divergence, grad


def _kl_divergence_error(params, P, neighbors, degrees_of_freedom, n_samples,
                         n_components):
    """t-SNE objective function: the absolute error of the
    KL divergence of p_ijs and q_ijs.

    Parameters
    ----------
    params : array, shape (n_params,)
        Unraveled embedding.

    P : array, shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.

    neighbors : array (n_samples, K)
        The neighbors is not actually required to calculate the
        divergence, but is here to match the signature of the
        gradient function

    degrees_of_freedom : float
        Degrees of freedom of the Student's-t distribution.

    n_samples : int
        Number of samples.

    n_components : int
        Dimension of the embedded space.

    Returns
    -------
    kl_divergence : float
        Kullback-Leibler divergence of p_ij and q_ij.

    grad : array, shape (n_params,)
        Unraveled gradient of the Kullback-Leibler divergence with respect to
        the embedding.
    """
    X_embedded = params.reshape(n_samples, n_components)

    # Q is a heavy-tailed distribution: Student's t-distribution
    n = pdist(X_embedded, "sqeuclidean")
    n += 1.
    n /= degrees_of_freedom
    n **= (degrees_of_freedom + 1.0) / -2.0
    Q = np.maximum(n / (2.0 * np.sum(n)), MACHINE_EPSILON)

    # Optimization trick below: np.dot(x, y) is faster than
    # np.sum(x * y) because it calls BLAS

    # Objective: C (Kullback-Leibler divergence of P and Q)
    if len(P.shape) == 2:
        P = squareform(P)
    kl_divergence = 2.0 * np.dot(P, np.log(P / Q))

    return kl_divergence


def _kl_divergence_bh(params, P, neighbors, degrees_of_freedom, n_samples,
                      n_components, angle=0.5, skip_num_points=0,
                      verbose=False):
    """t-SNE objective function: KL divergence of p_ijs and q_ijs.

    Uses Barnes-Hut tree methods to calculate the gradient that
    runs in O(NlogN) instead of O(N^2)

    Parameters
    ----------
    params : array, shape (n_params,)
        Unraveled embedding.

    P : array, shape (n_samples * (n_samples-1) / 2,)
        Condensed joint probability matrix.

    neighbors: int64 array, shape (n_samples, K)
        Array with element [i, j] giving the index for the jth
        closest neighbor to point i.

    degrees_of_freedom : float
        Degrees of freedom of the Student's-t distribution.

    n_samples : int
        Number of samples.

    n_components : int
        Dimension of the embedded space.

    angle : float (default: 0.5)
        This is the trade-off between speed and accuracy for Barnes-Hut T-SNE.
        'angle' is the angular size (referred to as theta in [3]) of a distant
        node as measured from a point. If this size is below 'angle' then it is
        used as a summary node of all points contained within it.
        This method is not very sensitive to changes in this parameter
        in the range of 0.2 - 0.8. Angle less than 0.2 has quickly increasing
        computation time and angle greater 0.8 has quickly increasing error.

    skip_num_points : int (optional, default:0)
        This does not compute the gradient for points with indices below
        `skip_num_points`. This is useful when computing transforms of new
        data where you'd like to keep the old data fixed.

    verbose : int
        Verbosity level.

    Returns
    -------
    kl_divergence : float
        Kullback-Leibler divergence of p_ij and q_ij.

    grad : array, shape (n_params,)
        Unraveled gradient of the Kullback-Leibler divergence with respect to
        the embedding.
    """
    params = astype(params, np.float32, copy=False)
    X_embedded = params.reshape(n_samples, n_components)
    neighbors = astype(neighbors, np.int64, copy=False)
    if len(P.shape) == 1:
        sP = squareform(P).astype(np.float32)
    else:
        sP = P.astype(np.float32)

    grad = np.zeros(X_embedded.shape, dtype=np.float32)
    error = _barnes_hut_tsne.gradient(sP, X_embedded, neighbors,
                                      grad, angle, n_components, verbose,
                                      dof=degrees_of_freedom)
    c = 2.0 * (degrees_of_freedom + 1.0) / degrees_of_freedom
    grad = grad.ravel()
    grad *= c

    return error, grad


def _gradient_descent(objective, p0, it, n_iter, objective_error=None,
                      n_iter_check=1, n_iter_without_progress=50,
                      momentum=0.5, learning_rate=1000.0, min_gain=0.01,
                      min_grad_norm=1e-7, min_error_diff=1e-7, verbose=0,
                      args=None, kwargs=None):
    """Batch gradient descent with momentum and individual gains.

    Parameters
    ----------
    objective : function or callable
        Should return a tuple of cost and gradient for a given parameter
        vector. When expensive to compute, the cost can optionally
        be None and can be computed every n_iter_check steps using
        the objective_error function.

    p0 : array-like, shape (n_params,)
        Initial parameter vector.

    it : int
        Current number of iterations (this function will be called more than
        once during the optimization).

    n_iter : int
        Maximum number of gradient descent iterations.

    n_iter_check : int
        Number of iterations before evaluating the global error. If the error
        is sufficiently low, we abort the optimization.

    objective_error : function or callable
        Should return a tuple of cost and gradient for a given parameter
        vector.

    n_iter_without_progress : int, optional (default: 30)
        Maximum number of iterations without progress before we abort the
        optimization.

    momentum : float, within (0.0, 1.0), optional (default: 0.5)
        The momentum generates a weight for previous gradients that decays
        exponentially.

    learning_rate : float, optional (default: 1000.0)
        The learning rate should be extremely high for t-SNE! Values in the
        range [100.0, 1000.0] are common.

    min_gain : float, optional (default: 0.01)
        Minimum individual gain for each parameter.

    min_grad_norm : float, optional (default: 1e-7)
        If the gradient norm is below this threshold, the optimization will
        be aborted.

    min_error_diff : float, optional (default: 1e-7)
        If the absolute difference of two successive cost function values
        is below this threshold, the optimization will be aborted.

    verbose : int, optional (default: 0)
        Verbosity level.

    args : sequence
        Arguments to pass to objective function.

    kwargs : dict
        Keyword arguments to pass to objective function.

    Returns
    -------
    p : array, shape (n_params,)
        Optimum parameters.

    error : float
        Optimum.

    i : int
        Last iteration.
    """
    if args is None:
        args = []
    if kwargs is None:
        kwargs = {}

    p = p0.copy().ravel()
    update = np.zeros_like(p)
    gains = np.ones_like(p)
    error = np.finfo(np.float).max
    best_error = np.finfo(np.float).max
    best_iter = 0

    for i in range(it, n_iter):
        new_error, grad = objective(p, *args, **kwargs)
        grad_norm = linalg.norm(grad)

        inc = update * grad >= 0.0
        dec = np.invert(inc)
        gains[inc] += 0.05
        gains[dec] *= 0.95
        np.clip(gains, min_gain, np.inf)
        grad *= gains
        update = momentum * update - learning_rate * grad
        p += update

        if (i + 1) % n_iter_check == 0:
            if new_error is None:
                new_error = objective_error(p, *args)
            error_diff = np.abs(new_error - error)
            error = new_error

            if verbose >= 2:
                m = "[t-SNE] Iteration %d: error = %.7f, gradient norm = %.7f"
                print(m % (i + 1, error, grad_norm))

            if error < best_error:
                best_error = error
                best_iter = i
            elif i - best_iter > n_iter_without_progress:
                if verbose >= 2:
                    print("[t-SNE] Iteration %d: did not make any progress "
                          "during the last %d episodes. Finished."
                          % (i + 1, n_iter_without_progress))
                break
            if grad_norm <= min_grad_norm:
                if verbose >= 2:
                    print("[t-SNE] Iteration %d: gradient norm %f. Finished."
                          % (i + 1, grad_norm))
                break
            if error_diff <= min_error_diff:
                if verbose >= 2:
                    m = "[t-SNE] Iteration %d: error difference %f. Finished."
                    print(m % (i + 1, error_diff))
                break

        if new_error is not None:
            error = new_error

    return p, error, i


def trustworthiness(X, X_embedded, n_neighbors=5, precomputed=False):
    """Expresses to what extent the local structure is retained.

    The trustworthiness is within [0, 1]. It is defined as

    .. math::

        T(k) = 1 - \frac{2}{nk (2n - 3k - 1)} \sum^n_{i=1}
            \sum_{j \in U^{(k)}_i (r(i, j) - k)}

    where :math:`r(i, j)` is the rank of the embedded datapoint j
    according to the pairwise distances between the embedded datapoints,
    :math:`U^{(k)}_i` is the set of points that are in the k nearest
    neighbors in the embedded space but not in the original space.

    * "Neighborhood Preservation in Nonlinear Projection Methods: An
      Experimental Study"
      J. Venna, S. Kaski
    * "Learning a Parametric Embedding by Preserving Local Structure"
      L.J.P. van der Maaten

    Parameters
    ----------
    X : array, shape (n_samples, n_features) or (n_samples, n_samples)
        If the metric is 'precomputed' X must be a square distance
        matrix. Otherwise it contains a sample per row.

    X_embedded : array, shape (n_samples, n_components)
        Embedding of the training data in low-dimensional space.

    n_neighbors : int, optional (default: 5)
        Number of neighbors k that will be considered.

    precomputed : bool, optional (default: False)
        Set this flag if X is a precomputed square distance matrix.

    Returns
    -------
    trustworthiness : float
        Trustworthiness of the low-dimensional embedding.
    """
    if precomputed:
        dist_X = X
    else:
        dist_X = pairwise_distances(X, squared=True)
    dist_X_embedded = pairwise_distances(X_embedded, squared=True)
    ind_X = np.argsort(dist_X, axis=1)
    ind_X_embedded = np.argsort(dist_X_embedded, axis=1)[:, 1:n_neighbors + 1]

    n_samples = X.shape[0]
    t = 0.0
    ranks = np.zeros(n_neighbors)
    for i in range(n_samples):
        for j in range(n_neighbors):
            ranks[j] = np.where(ind_X[i] == ind_X_embedded[i, j])[0][0]
        ranks -= n_neighbors
        t += np.sum(ranks[ranks > 0])
    t = 1.0 - t * (2.0 / (n_samples * n_neighbors *
                          (2.0 * n_samples - 3.0 * n_neighbors - 1.0)))
    return t


class TSNE(BaseEstimator):
    """t-distributed Stochastic Neighbor Embedding.

    t-SNE [1] is a tool to visualize high-dimensional data. It converts
    similarities between data points to joint probabilities and tries
    to minimize the Kullback-Leibler divergence between the joint
    probabilities of the low-dimensional embedding and the
    high-dimensional data. t-SNE has a cost function that is not convex,
    i.e. with different initializations we can get different results.

    It is highly recommended to use another dimensionality reduction
    method (e.g. PCA for dense data or TruncatedSVD for sparse data)
    to reduce the number of dimensions to a reasonable amount (e.g. 50)
    if the number of features is very high. This will suppress some
    noise and speed up the computation of pairwise distances between
    samples. For more tips see Laurens van der Maaten's FAQ [2].

    Read more in the :ref:`User Guide <t_sne>`.

    Parameters
    ----------
    n_components : int, optional (default: 2)
        Dimension of the embedded space.

    perplexity : float, optional (default: 30)
        The perplexity is related to the number of nearest neighbors that
        is used in other manifold learning algorithms. Larger datasets
        usually require a larger perplexity. Consider selecting a value
        between 5 and 50. The choice is not extremely critical since t-SNE
        is quite insensitive to this parameter.

    early_exaggeration : float, optional (default: 4.0)
        Controls how tight natural clusters in the original space are in
        the embedded space and how much space will be between them. For
        larger values, the space between natural clusters will be larger
        in the embedded space. Again, the choice of this parameter is not
        very critical. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high.

    learning_rate : float, optional (default: 1000)
        The learning rate can be a critical parameter. It should be
        between 100 and 1000. If the cost function increases during initial
        optimization, the early exaggeration factor or the learning rate
        might be too high. If the cost function gets stuck in a bad local
        minimum increasing the learning rate helps sometimes.

    n_iter : int, optional (default: 1000)
        Maximum number of iterations for the optimization. Should be at
        least 200.

    n_iter_without_progress : int, optional (default: 30)
        Only used if method='exact'
        Maximum number of iterations without progress before we abort the
        optimization. If method='barnes_hut' this parameter is fixed to
        a value of 30 and cannot be changed.

        .. versionadded:: 0.17
           parameter *n_iter_without_progress* to control stopping criteria.

    min_grad_norm : float, optional (default: 1e-7)
        Only used if method='exact'
        If the gradient norm is below this threshold, the optimization will
        be aborted. If method='barnes_hut' this parameter is fixed to a value
        of 1e-3 and cannot be changed.

    metric : string or callable, optional
        The metric to use when calculating distance between instances in a
        feature array. If metric is a string, it must be one of the options
        allowed by scipy.spatial.distance.pdist for its metric parameter, or
        a metric listed in pairwise.PAIRWISE_DISTANCE_FUNCTIONS.
        If metric is "precomputed", X is assumed to be a distance matrix.
        Alternatively, if metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays from X as input and return a value indicating
        the distance between them. The default is "euclidean" which is
        interpreted as squared euclidean distance.

    init : string or numpy array, optional (default: "random")
        Initialization of embedding. Possible options are 'random', 'pca',
        and a numpy array of shape (n_samples, n_components).
        PCA initialization cannot be used with precomputed distances and is
        usually more globally stable than random initialization.

    verbose : int, optional (default: 0)
        Verbosity level.

    random_state : int or RandomState instance or None (default)
        Pseudo Random Number generator seed control. If None, use the
        numpy.random singleton. Note that different initializations
        might result in different local minima of the cost function.

    method : string (default: 'barnes_hut')
        By default the gradient calculation algorithm uses Barnes-Hut
        approximation running in O(NlogN) time. method='exact'
        will run on the slower, but exact, algorithm in O(N^2) time. The
        exact algorithm should be used when nearest-neighbor errors need
        to be better than 3%. However, the exact method cannot scale to
        millions of examples.

        .. versionadded:: 0.17
           Approximate optimization *method* via the Barnes-Hut.

    angle : float (default: 0.5)
        Only used if method='barnes_hut'
        This is the trade-off between speed and accuracy for Barnes-Hut T-SNE.
        'angle' is the angular size (referred to as theta in [3]) of a distant
        node as measured from a point. If this size is below 'angle' then it is
        used as a summary node of all points contained within it.
        This method is not very sensitive to changes in this parameter
        in the range of 0.2 - 0.8. Angle less than 0.2 has quickly increasing
        computation time and angle greater 0.8 has quickly increasing error.


    Attributes
    ----------
    embedding_ : array-like, shape (n_samples, n_components)
        Stores the embedding vectors.

    kl_divergence_ : float
        Kullback-Leibler divergence after optimization.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.manifold import TSNE
    >>> X = np.array([[0, 0, 0], [0, 1, 1], [1, 0, 1], [1, 1, 1]])
    >>> model = TSNE(n_components=2, random_state=0)
    >>> np.set_printoptions(suppress=True)
    >>> model.fit_transform(X) # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    array([[ 0.00017599,  0.00003993],
           [ 0.00009891,  0.00021913],
           [ 0.00018554, -0.00009357],
           [ 0.00009528, -0.00001407]])

    References
    ----------

    [1] van der Maaten, L.J.P.; Hinton, G.E. Visualizing High-Dimensional Data
        Using t-SNE. Journal of Machine Learning Research 9:2579-2605, 2008.

    [2] van der Maaten, L.J.P. t-Distributed Stochastic Neighbor Embedding
        http://homepage.tudelft.nl/19j49/t-SNE.html

    [3] L.J.P. van der Maaten. Accelerating t-SNE using Tree-Based Algorithms.
        Journal of Machine Learning Research 15(Oct):3221-3245, 2014.
        http://lvdmaaten.github.io/publications/papers/JMLR_2014.pdf
    """

    def __init__(self, n_components=2, perplexity=30.0,
                 early_exaggeration=4.0, learning_rate=1000.0, n_iter=1000,
                 n_iter_without_progress=30, min_grad_norm=1e-7,
                 metric="euclidean", init="random", verbose=0,
                 random_state=None, method='barnes_hut', angle=0.5):
        if not ((isinstance(init, string_types) and
                init in ["pca", "random"]) or
                isinstance(init, np.ndarray)):
            msg = "'init' must be 'pca', 'random', or a numpy array"
            raise ValueError(msg)
        self.n_components = n_components
        self.perplexity = perplexity
        self.early_exaggeration = early_exaggeration
        self.learning_rate = learning_rate
        self.n_iter = n_iter
        self.n_iter_without_progress = n_iter_without_progress
        self.min_grad_norm = min_grad_norm
        self.metric = metric
        self.init = init
        self.verbose = verbose
        self.random_state = random_state
        self.method = method
        self.angle = angle
        self.embedding_ = None

    def _fit(self, X, skip_num_points=0):
        """Fit the model using X as training data.

        Note that sparse arrays can only be handled by method='exact'.
        It is recommended that you convert your sparse array to dense
        (e.g. `X.toarray()`) if it fits in memory, or otherwise using a
        dimensionality reduction technique (e.g. TruncatedSVD).

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the metric is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row. Note that this
            when method='barnes_hut', X cannot be a sparse array and if need be
            will be converted to a 32 bit float array. Method='exact' allows
            sparse arrays and 64bit floating point inputs.

        skip_num_points : int (optional, default:0)
            This does not compute the gradient for points with indices below
            `skip_num_points`. This is useful when computing transforms of new
            data where you'd like to keep the old data fixed.
        """
        if self.method not in ['barnes_hut', 'exact']:
            raise ValueError("'method' must be 'barnes_hut' or 'exact'")
        if self.angle < 0.0 or self.angle > 1.0:
            raise ValueError("'angle' must be between 0.0 - 1.0")
        if self.method == 'barnes_hut' and sp.issparse(X):
            raise TypeError('A sparse matrix was passed, but dense '
                            'data is required for method="barnes_hut". Use '
                            'X.toarray() to convert to a dense numpy array if '
                            'the array is small enough for it to fit in '
                            'memory. Otherwise consider dimensionality '
                            'reduction techniques (e.g. TruncatedSVD)')
        else:
            X = check_array(X, accept_sparse=['csr', 'csc', 'coo'],
                            dtype=np.float64)
        random_state = check_random_state(self.random_state)

        if self.early_exaggeration < 1.0:
            raise ValueError("early_exaggeration must be at least 1, but is "
                             "%f" % self.early_exaggeration)

        if self.n_iter < 200:
            raise ValueError("n_iter should be at least 200")

        if self.metric == "precomputed":
            if isinstance(self.init, string_types) and self.init == 'pca':
                raise ValueError("The parameter init=\"pca\" cannot be used "
                                 "with metric=\"precomputed\".")
            if X.shape[0] != X.shape[1]:
                raise ValueError("X should be a square distance matrix")
            distances = X
        else:
            if self.verbose:
                print("[t-SNE] Computing pairwise distances...")

            if self.metric == "euclidean":
                distances = pairwise_distances(X, metric=self.metric,
                                               squared=True)
            else:
                distances = pairwise_distances(X, metric=self.metric)

        if not np.all(distances >= 0):
            raise ValueError("All distances should be positive, either "
                             "the metric or precomputed distances given "
                             "as X are not correct")

        # Degrees of freedom of the Student's t-distribution. The suggestion
        # degrees_of_freedom = n_components - 1 comes from
        # "Learning a Parametric Embedding by Preserving Local Structure"
        # Laurens van der Maaten, 2009.
        degrees_of_freedom = max(self.n_components - 1.0, 1)
        n_samples = X.shape[0]
        # the number of nearest neighbors to find
        k = min(n_samples - 1, int(3. * self.perplexity + 1))

        neighbors_nn = None
        if self.method == 'barnes_hut':
            if self.verbose:
                print("[t-SNE] Computing %i nearest neighbors..." % k)
            if self.metric == 'precomputed':
                # Use the precomputed distances to find
                # the k nearest neighbors and their distances
                neighbors_nn = np.argsort(distances, axis=1)[:, :k]
            else:
                # Find the nearest neighbors for every point
                bt = BallTree(X)
                # LvdM uses 3 * perplexity as the number of neighbors
                # And we add one to not count the data point itself
                # In the event that we have very small # of points
                # set the neighbors to n - 1
                distances_nn, neighbors_nn = bt.query(X, k=k + 1)
                neighbors_nn = neighbors_nn[:, 1:]
            P = _joint_probabilities_nn(distances, neighbors_nn,
                                        self.perplexity, self.verbose)
        else:
            P = _joint_probabilities(distances, self.perplexity, self.verbose)
        assert np.all(np.isfinite(P)), "All probabilities should be finite"
        assert np.all(P >= 0), "All probabilities should be zero or positive"
        assert np.all(P <= 1), ("All probabilities should be less "
                                "or then equal to one")

        if isinstance(self.init, np.ndarray):
            X_embedded = self.init
        elif self.init == 'pca':
            pca = PCA(n_components=self.n_components, svd_solver='randomized',
                      random_state=random_state)
            X_embedded = pca.fit_transform(X)
        elif self.init == 'random':
            X_embedded = None
        else:
            raise ValueError("Unsupported initialization scheme: %s"
                             % self.init)

        return self._tsne(P, degrees_of_freedom, n_samples, random_state,
                          X_embedded=X_embedded,
                          neighbors=neighbors_nn,
                          skip_num_points=skip_num_points)

    def _tsne(self, P, degrees_of_freedom, n_samples, random_state,
              X_embedded=None, neighbors=None, skip_num_points=0):
        """Runs t-SNE."""
        # t-SNE minimizes the Kullback-Leiber divergence of the Gaussians P
        # and the Student's t-distributions Q. The optimization algorithm that
        # we use is batch gradient descent with three stages:
        # * early exaggeration with momentum 0.5
        # * early exaggeration with momentum 0.8
        # * final optimization with momentum 0.8
        # The embedding is initialized with iid samples from Gaussians with
        # standard deviation 1e-4.

        if X_embedded is None:
            # Initialize embedding randomly
            X_embedded = 1e-4 * random_state.randn(n_samples,
                                                   self.n_components)
        params = X_embedded.ravel()

        opt_args = {"n_iter": 50, "momentum": 0.5, "it": 0,
                    "learning_rate": self.learning_rate,
                    "n_iter_without_progress": self.n_iter_without_progress,
                    "verbose": self.verbose, "n_iter_check": 25,
                    "kwargs": dict(skip_num_points=skip_num_points)}
        if self.method == 'barnes_hut':
            m = "Must provide an array of neighbors to use Barnes-Hut"
            assert neighbors is not None, m
            obj_func = _kl_divergence_bh
            objective_error = _kl_divergence_error
            sP = squareform(P).astype(np.float32)
            neighbors = neighbors.astype(np.int64)
            args = [sP, neighbors, degrees_of_freedom, n_samples,
                    self.n_components]
            opt_args['args'] = args
            opt_args['min_grad_norm'] = 1e-3
            opt_args['n_iter_without_progress'] = 30
            # Don't always calculate the cost since that calculation
            # can be nearly as expensive as the gradient
            opt_args['objective_error'] = objective_error
            opt_args['kwargs']['angle'] = self.angle
            opt_args['kwargs']['verbose'] = self.verbose
        else:
            obj_func = _kl_divergence
            opt_args['args'] = [P, degrees_of_freedom, n_samples,
                                self.n_components]
            opt_args['min_error_diff'] = 0.0
            opt_args['min_grad_norm'] = self.min_grad_norm

        # Early exaggeration
        P *= self.early_exaggeration

        params, kl_divergence, it = _gradient_descent(obj_func, params,
                                                      **opt_args)
        opt_args['n_iter'] = 100
        opt_args['momentum'] = 0.8
        opt_args['it'] = it + 1
        params, kl_divergence, it = _gradient_descent(obj_func, params,
                                                      **opt_args)
        if self.verbose:
            print("[t-SNE] KL divergence after %d iterations with early "
                  "exaggeration: %f" % (it + 1, kl_divergence))
        # Save the final number of iterations
        self.n_iter_final = it

        # Final optimization
        P /= self.early_exaggeration
        opt_args['n_iter'] = self.n_iter
        opt_args['it'] = it + 1
        params, error, it = _gradient_descent(obj_func, params, **opt_args)

        if self.verbose:
            print("[t-SNE] Error after %d iterations: %f"
                  % (it + 1, kl_divergence))

        X_embedded = params.reshape(n_samples, self.n_components)
        self.kl_divergence_ = kl_divergence

        return X_embedded

    def fit_transform(self, X, y=None):
        """Fit X into an embedded space and return that transformed
        output.

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the metric is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row.

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Embedding of the training data in low-dimensional space.
        """
        embedding = self._fit(X)
        self.embedding_ = embedding
        return self.embedding_

    def fit(self, X, y=None):
        """Fit X into an embedded space.

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the metric is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row. If the method
            is 'exact', X may be a sparse matrix of type 'csr', 'csc'
            or 'coo'.
        """
        self.fit_transform(X)
        return self
