import numpy as np
from scipy import linalg
from scipy.optimize import fmin_l_bfgs_b
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform
from ..base import BaseEstimator, TransformerMixin
from ..utils import check_random_state
from ..neighbors import NearestNeighbors
from ..neighbors.base import _get_weights
from . import _binary_search


EPSILON = 1e-16


def _gradient_descent(objective, p0, n_iter, momentum=0.5,
                      learning_rate=100.0, min_gain=0.01,
                      min_grad_norm=1e-6, min_error_diff=1e-5, verbose=0):
    """Batch gradient descent with momentum and individual gains.

    Parameters
    ----------
    objective : function or callable
        Should return a tuple of cost and gradient for a given parameter
        vector.

    p0 : array-like, shape (n_params,)
        Initial parameter vector.

    n_iter : int
        Maximum number of gradient descent iterations.

    momentum : float, within (0.0, 1.0), optional (default: 0.5)
        The momentum generates a weight for previous gradients that decays
        exponentially.

    learning_rate : float, optional (default: 100.0)
        The learning rate should be extremely high for t-SNE! Values in the
        range [100.0, 500.0] are common.

    min_gain : float, optional (default: 0.01)
        Minimum individual gain for each parameter.

    min_grad_norm : float, optional (default: 1e-6)
        If the gradient norm is below this threshold, the optimization will
        be aborted.

    min_error_diff : float, optional (default: 1e-5)
        If the absolute difference of two successive cost function values
        is below this threshold, the optimization will be aborted.

    verbose : int, optional (default: 0)
        Verbosity level.
    """
    p = p0.copy().ravel()
    update = np.zeros_like(p)
    gains = np.ones_like(p)
    error = np.finfo(np.float).max
    for i in range(n_iter):
        new_error, grad = objective(p)
        error_diff = np.abs(new_error - error)
        error = new_error
        grad_norm = linalg.norm(grad)

        if min_grad_norm >= grad_norm:
            if verbose:
                print("Iteration %d: gradient norm %f. Finished."
                      % (i + 1, grad_norm))
            break
        if min_error_diff >= error_diff:
            if verbose:
                print("Iteration %d: error difference %f. Finished."
                      % (i + 1, error_diff))
            break

        inc = update * grad >= 0.0
        dec = np.invert(inc)
        gains[inc] += 0.05
        gains[dec] *= 0.95;
        np.clip(gains, min_gain, np.inf)
        grad *= gains
        update = momentum * update - learning_rate * grad
        p += update

        if verbose and (i+1) % 10 == 0:
            print("Iteration %d: error = %.6f, gradient norm = %.6f"
                  % (i + 1, error, grad_norm))

    return p, error


def trustworthiness(X, X_embedded, n_neighbors=5, precomputed=False):
    """Expresses to what extent the local structure is retained.

    The trustworthiness is within [0, 1]. It is defined as

    .. math::

        T(k) = 1 - \frac{2}{nk (2n - 3k - 1)} \sum^n_{i=1} \sum_{j \in U^{(k)}_i (r(i, j) - k)}

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
        If the distance is 'precomputed' X must be a square distance
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
        dist_X = squareform(pdist(X, "sqeuclidean"))
    dist_X_embedded = squareform(pdist(X_embedded, "sqeuclidean"))
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


class TSNE(BaseEstimator, TransformerMixin):
    """t-distributed Stochastic Neighbor Embedding.

    t-SNE is a tool to visualize high-dimensional data. It converts
    similarities between data points to joint probabilities and tries
    to minimize the Kullback-Leibler divergence between the joint
    probabilities of the low-dimensional embedding and the
    high-dimensional data.

    It is highly recommended to use another dimensionality reduction
    method (e.g. PCA) to reduce the number of dimensions to a reasonable
    amount (e.g. 50) if the number of features is very high. This will
    often improve the visualization.

    Usually t-SNE does not generalize, i.e. it can only compute the
    embedding of the training data. However, we use a heuristic to
    embed unseen data: first we determine the n_neighbors nearest
    neighbors from the training set of the test sample in the original
    space and then we compute a distance-weighted average in the
    embedded space to obtain the transformed test sample. Note that
    this does not work if the distance matrix is precomputed.

    .. topic:: References:

        * `"t-Distributed Stochastic Neighbor Embedding"
        <http://homepage.tudelft.nl/19j49/t-SNE.html>`_
        L.J.P. van der Maaten
        * `"Visualizing High-Dimensional Data Using t-SNE"
        <http://jmlr.csail.mit.edu/papers/volume9/vandermaaten08a/vandermaaten08a.pdf>`_
        L.J.P. van der Maaten, G.E. Hinton, 2008

    Parameters
    ----------
    n_components : int, optional (default: 2)
        Number of components to keep.

    perplexity : float, optional (default: 30)
        The perplexity is related to the number of nearest neighbors that is
        used in other manifold learning algorithms. Larger datasets usually
        require a larger perplexity. Consider selcting a value between 5 and
        50.

    optimizer : string
        One of 'gradient_descent' and 'lbfgs'.

    n_iter : int, optional (default: 1000)
        Maximum number of iterations for the optimizer.

    distances : string, optional (default: sqeuclidean)
        A distance metric that is defined in scipy.spatial.distance or
        'precomputed'.

    n_neighbors : int, optional, (default: 3)
        Number of neighbors that will be used to transform the data to the
        embedded space.

    fit_inverse_transform : bool, optional, (default: False)
        Learn the inverse transform for non-precomputed distances.

    verbose : int, optional (default: 0)
        Verbosity level.

    random_state : int or RandomState instance or None (default)
        Pseudo Random Number generator seed control. If None, use the
        numpy.random singleton.
    """
    def __init__(self, n_components=2, perplexity=30.0,
                 optimizer="gradient_descent", n_iter=1000,
                 distances="sqeuclidean", n_neighbors=3,
                 fit_inverse_transform=False, verbose=0, random_state=None):
        if fit_inverse_transform and distances == "precomputed":
            raise ValueError("Cannot fit_inverse_transform with a precomputed "
                             "distance matrix.")
        self.n_components = n_components
        self.perplexity = perplexity
        self.optimizer = optimizer
        self.n_iter = n_iter
        self.distances = distances
        self.n_neighbors = n_neighbors
        self.fit_inverse_transform = fit_inverse_transform
        self.verbose = verbose
        self.random_state = random_state

    def fit(self, X, y=None):
        """Fit the model using X as training data.

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the distance is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row.

        Returns
        -------
        self : TSNE
            This object.
        """
        if not self.optimizer in ["gradient_descent", "lbfgs"]:
            raise ValueError("Unknown optimizer '%s'" % self.optimizer)

        if self.distances == "precomputed":
            if X.shape[0] != X.shape[1]:
                raise ValueError("X should be a square distance matrix")
            dist = X
        else:
            if self.verbose:
                print("Computing pairwise distances...")
            dist = squareform(pdist(X, self.distances))

        self.embedding_ = self._tsne(dist)

        if self.distances != "precomputed":
            self.knn_ = NearestNeighbors(n_neighbors=self.n_neighbors)
            self.knn_.fit(X)

        if self.fit_inverse_transform:
            self.knn_embedded_ = NearestNeighbors(n_neighbors=self.n_neighbors)
            self.knn_embedded_.fit(self.embedding_)

        return self

    def _tsne(self, dist):
        """Runs t-SNE."""
        n_samples = dist.shape[0]
        # Degrees of freedom of the Student's t-distribution. The suggestion
        # alpha = n_components - 1 comes from "Learning a Parametric Embedding
        # by Preserving Local Structure" Laurens van der Maaten, 2009.
        alpha = self.n_components - 1.0
        random_state = check_random_state(self.random_state)

        # Compute conditional probabilities such that they approximately match
        # the desired perplexity
        conditional_P = _binary_search._binary_search_perplexity(dist,
            self.perplexity, self.verbose)
        P = conditional_P + conditional_P.T
        # Normalization and "early exaggeration" (* 4, see below for an
        # explanation)
        sum_P = np.maximum(np.sum(P), EPSILON)
        P = np.maximum(4.0 * squareform(P) / sum_P, EPSILON)
        # Note that P and Q are stored as condensed matrices, i.e. we omit
        # the diagonal and duplicate entries and store everything in a
        # one-dimensional array.

        # Initialize embedding
        Y = random_state.randn(n_samples, self.n_components)

        # Objective function: KL divergence of pij and qij
        def objective(p):
            Y = p.reshape(n_samples, self.n_components)

            # Q is a heavy-tailed distribution: Student's t-distribution
            n = ((1.0 + pdist(Y, "sqeuclidean") / alpha)
                 ** ((alpha + 1.0) / -2.0))
            Q = np.maximum(n / (2.0 * np.sum(n)), EPSILON)

            # Gradient: dC/dY
            grad = np.ndarray((n_samples, self.n_components))
            PQd = squareform((P - Q) * n)
            c = 2.0 * (alpha + 1.0) / alpha
            for i in range(n_samples):
                grad[i] = c * np.sum(PQd[i].reshape(-1, 1) * (Y[i] - Y), axis=0)
            grad = grad.ravel()

            # Objective: C (Kullback-Leibler divergence of P and Q)
            kl_divergence = 2.0 * np.sum(P * np.log(P / Q))

            return kl_divergence, grad

        if self.optimizer == "lbfgs":
            optimizer = fmin_l_bfgs_b
            pre_opt_args = dict(disp=self.verbose, maxfun=10)
            opt_args = dict(disp=self.verbose, maxfun=self.n_iter)
        elif self.optimizer == "gradient_descent":
            optimizer = _gradient_descent
            pre_opt_args = dict(verbose=self.verbose, momentum=0.5,
                                learning_rate=100.0, n_iter=100)
            opt_args = dict(verbose=self.verbose, momentum=0.9,
                            learning_rate=100.0, n_iter=self.n_iter)

        # Early exaggeration (from "Visualizing Data using t-SNE" by
        # Laurens van der Maaten and Geoffrey Hinton):
        # "A less obvious way to improve the optimization, which we call
        # 'early exaggeration', is to multiply all of the pij's by, for
        # example, 4, in the initial stages of the optimization. This means
        # that almost all of the qij's, which still add up to 1, are much
        # too small to model their corresponding pij's. As a result, the
        # optimization is encouraged to focus on modeling the large pij's
        # by fairly large qij's. The effect is that the natural clusters in
        # the data tend to form tight widely separated clusters in the map.
        # This creates a lot of relatively empty space in the map, which
        # makes it much easier for the clusters to move around relative to
        # one another in order to find a good global organization."
        p = Y.ravel()
        r = optimizer(objective, p, **pre_opt_args)
        p, error = r[:2]
        if self.verbose:
            print("Error after initial optimization: %f" % error)

        P /= 4.0  # Undo exaggeration
        r = optimizer(objective, p, **opt_args)
        p, error = r[:2]
        if self.verbose:
            print("Error after final optimization: %f" % error)

        Y = p.reshape(n_samples, self.n_components)

        return Y

    def score(self, X, y=None, n_neighbors=5):
        """Compute trustworthiness.

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the distance is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row.

        n_neighbors : int, optional (default: 5)
            Number of neighbors k that will be considered.

        Returns
        -------
        trustworthiness : float
            Trustworthiness of the low-dimensional embedding.
        """
        return trustworthiness(X, self.transform(X), n_neighbors=n_neighbors,
                               precomputed=self.distances == "precomputed")

    def inverse_transform(self, X):
        """Transform X back to original space.

        Parameters
        ----------
        X : array, shape (n_samples, n_components)
            Embedded data.

        Returns
        -------
        X_new : array, shape (n_samples, n_features)
            Data in the original space.
        """
        if not self.fit_inverse_transform:
            raise ValueError("Inverse transform was not fitted!")

        neigh_dist, neigh_ind = self.knn_embedded_.kneighbors(X)
        neigh_dist = np.array([np.maximum(d, EPSILON) for d in neigh_dist])

        weights = _get_weights(neigh_dist, "distance")

        if weights is None:
            X_original = np.array([np.mean(self.X[ind], axis=0)
                                   for ind in neigh_ind])
        else:
            X_original = np.array([(np.average(self.X[ind], axis=0,
                                               weights=weights[i]))
                                   for (i, ind) in enumerate(neigh_ind)])

        return X_original

    def transform(self, X):
        """Transform X to the embedded space.

        Parameters
        ----------
        X : array, shape (n_samples, n_features) or (n_samples, n_samples)
            If the distance is 'precomputed' X must be a square distance
            matrix. Otherwise it contains a sample per row.

        Returns
        -------
        X_new : array, shape (n_samples, n_components)
            Embedding of the training data in low-dimensional space.
        """
        if self.distances == "precomputed":
            if X is self.X:
                return self.embedding_
            else:
                raise ValueError("Can only transform trainin data when "
                                 "distances are precomputed.")
        neigh_dist, neigh_ind = self.knn_.kneighbors(X)
        neigh_dist = np.array([np.maximum(d, EPSILON) for d in neigh_dist])

        weights = _get_weights(neigh_dist, "distance")

        if weights is None:
            X_embedded = np.array([np.mean(self.embedding_[ind], axis=0)
                                   for ind in neigh_ind])
        else:
            X_embedded = np.array([(np.average(self.embedding_[ind], axis=0,
                                               weights=weights[i]))
                                   for (i, ind) in enumerate(neigh_ind)])

        return X_embedded
