"""Fuzzy C-Means (FCM) clustering."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
from numbers import Integral, Real

import numpy as np

from sklearn.base import (
    BaseEstimator,
    ClusterMixin,
    TransformerMixin,
    _fit_context,
)
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.utils import check_array, check_random_state
from sklearn.utils._param_validation import Interval, StrOptions, validate_params
from sklearn.utils.validation import (
    _check_feature_names_in,
    _is_arraylike_not_scalar,
    check_is_fitted,
    validate_data,
)


def _dist_to_membership(D, m):
    """Convert a distance matrix to a fuzzy membership matrix.

    Parameters
    ----------
    D : ndarray of shape (n_samples, n_clusters)
        The distance matrix.

    m : float
        Fuzziness exponent (m > 1).

    Returns
    -------
    U : ndarray of shape (n_samples, n_clusters)
        The membership matrix.
    """
    p = 2.0 / (m - 1.0)

    # Handle zero distances to avoid DivisionByZero / NaNs
    zero_mask = D == 0
    zero_rows = np.any(zero_mask, axis=1)

    U = np.zeros_like(D)

    if np.any(zero_rows):
        # Assign membership equally among the centers with distance 0, and 0 to others
        num_zeros = np.sum(zero_mask, axis=1, keepdims=True)
        num_zeros_safe = np.where(num_zeros == 0, 1, num_zeros)
        U[zero_rows] = zero_mask[zero_rows] / num_zeros_safe[zero_rows]

        # For rows with non-zero distance, compute the standard formula
        non_zero_rows = ~zero_rows
        if np.any(non_zero_rows):
            D_sub = D[non_zero_rows]
            inv_D = 1.0 / D_sub
            inv_D_p = np.power(inv_D, p)
            sum_inv_D_p = np.sum(inv_D_p, axis=1, keepdims=True)
            U[non_zero_rows] = inv_D_p / sum_inv_D_p
    else:
        inv_D = 1.0 / D
        inv_D_p = np.power(inv_D, p)
        sum_inv_D_p = np.sum(inv_D_p, axis=1, keepdims=True)
        U = inv_D_p / sum_inv_D_p

    return U


@validate_params(
    {
        "X": ["array-like"],
        "n_clusters": [Interval(Integral, 1, None, closed="left")],
        "m": [Interval(Real, 1.0, None, closed="neither")],
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "tol": [Interval(Real, 0.0, None, closed="left")],
        "init": [
            StrOptions({"random_membership", "k-means++", "random_centers"}),
            callable,
            "array-like",
        ],
        "random_state": ["random_state"],
    },
    prefer_skip_nested_validation=True,
)
def fuzzy_c_means(
    X,
    n_clusters,
    *,
    m=2.0,
    max_iter=300,
    tol=1e-4,
    init="random_membership",
    random_state=None,
):
    """Perform Fuzzy C-Means (FCM) clustering.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        The data to cluster.

    n_clusters : int
        The number of clusters to form.

    m : float, default=2.0
        Fuzziness parameter (m > 1). Controls the degree of fuzziness of the clusters.

    max_iter : int, default=300
        Maximum number of iterations.

    tol : float, default=1e-4
        Tolerance for convergence (maximum change in membership matrix U).

    init : {'random_membership', 'k-means++', 'random_centers'}, callable or \
            array-like, default='random_membership'
        Method for initialization:
        - `'random_membership'`: initialize the fuzzy membership matrix U randomly.
        - `'k-means++'`: initialize cluster centers using k-means++ seeding.
        - `'random_centers'`: choose random points from X as initial centers.
        - array-like: shape (n_clusters, n_features) giving initial cluster centers.
        - callable: must accept X, n_clusters, random_state and return cluster centers.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation.

    Returns
    -------
    cluster_centers : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers.

    u : ndarray of shape (n_samples, n_clusters)
        Fuzzy membership matrix.

    labels : ndarray of shape (n_samples,)
        Hard labels (the cluster index with the maximum membership value).

    n_iter : int
        Number of iterations run.
    """
    X = check_array(X, accept_sparse=False, dtype=[np.float64, np.float32])
    n_samples, n_features = X.shape

    if n_samples < n_clusters:
        raise ValueError(f"n_samples={n_samples} should be >= n_clusters={n_clusters}.")

    random_state = check_random_state(random_state)

    # Initialization
    if isinstance(init, str) and init == "random_membership":
        U = random_state.uniform(size=(n_samples, n_clusters))
        U /= np.sum(U, axis=1, keepdims=True)
        # Compute initial centers from random membership
        Um = np.power(U, m)
        col_sums = np.sum(Um, axis=0, keepdims=True).T
        col_sums = np.where(col_sums == 0, 1e-15, col_sums)
        centers = (Um.T @ X) / col_sums
    else:
        if isinstance(init, str) and init == "k-means++":
            from sklearn.cluster._kmeans import _kmeans_plusplus

            x_squared_norms = np.sum(X**2, axis=1)
            sample_weight = np.ones(n_samples, dtype=X.dtype)
            centers, _ = _kmeans_plusplus(
                X, n_clusters, x_squared_norms, sample_weight, random_state
            )
        elif isinstance(init, str) and init == "random_centers":
            seeds = random_state.choice(n_samples, size=n_clusters, replace=False)
            centers = X[seeds]
        elif _is_arraylike_not_scalar(init):
            centers = np.array(init, dtype=X.dtype)
            if centers.shape != (n_clusters, n_features):
                raise ValueError(
                    f"The shape of initial centers {centers.shape} does not match "
                    f"the required shape {(n_clusters, n_features)}."
                )
        elif callable(init):
            centers = init(X, n_clusters, random_state=random_state)
            centers = check_array(centers, dtype=X.dtype)
            if centers.shape != (n_clusters, n_features):
                raise ValueError(
                    f"Callable init must return shape {(n_clusters, n_features)}, "
                    f"got {centers.shape}."
                )
        else:
            raise ValueError(f"Invalid init parameter: {init}")

        D = euclidean_distances(X, centers)
        U = _dist_to_membership(D, m)

    n_iter = 0
    for i in range(max_iter):
        U_old = U.copy()

        # 1. Update Centroids
        Um = np.power(U, m)
        col_sums = np.sum(Um, axis=0, keepdims=True).T
        col_sums = np.where(col_sums == 0, 1e-15, col_sums)
        centers = (Um.T @ X) / col_sums

        # 2. Update Membership
        D = euclidean_distances(X, centers)
        U = _dist_to_membership(D, m)

        n_iter += 1

        # Check convergence
        change = np.max(np.abs(U - U_old))
        if change < tol:
            break

    # Calculate hard labels from final membership matrix
    labels = np.argmax(U, axis=1)

    return centers, U, labels, n_iter


class FuzzyCMeans(TransformerMixin, ClusterMixin, BaseEstimator):
    """Fuzzy C-Means clustering.

    Read more in the :ref:`User Guide <fuzzy_c_means>`.

    Parameters
    ----------
    n_clusters : int, default=8
        The number of clusters to form as well as the number of centroids to generate.

    m : float, default=2.0
        Fuzziness parameter (m > 1). Controls the degree of fuzziness of the clusters.
        Higher values produce fuzzier, more overlapping clusters.

    max_iter : int, default=300
        Maximum number of iterations of the fuzzy c-means algorithm for a single run.

    tol : float, default=1e-4
        Relative tolerance with regards to the Frobenius norm of the difference
        in the cluster membership matrix U to declare convergence.

    init : {'random_membership', 'k-means++', 'random_centers'}, callable or \
            array-like, default='random_membership'
        Method for initialization:
        - `'random_membership'`: initialize the fuzzy membership matrix U randomly.
        - `'k-means++'`: initialize cluster centers using k-means++ seeding.
        - `'random_centers'`: choose random points from X as initial centers.
        - array-like: shape (n_clusters, n_features) giving initial cluster centers.
        - callable: must accept X, n_clusters, random_state and return cluster centers.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for initialization. Pass an int for
        reproducible outputs.

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers.

    u_ : ndarray of shape (n_samples, n_clusters)
        The final fuzzy membership matrix.

    labels_ : ndarray of shape (n_samples,)
        Hard labels of each point (the cluster with the maximum membership).

    n_iter_ : int
        Number of iterations run.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (n_features_in_,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    KMeans : The classic implementation of the clustering method based on the
        Lloyd's algorithm.

    Examples
    --------
    >>> from sklearn.cluster import FuzzyCMeans
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0]])
    >>> fcm = FuzzyCMeans(n_clusters=2, random_state=0).fit(X)
    >>> fcm.labels_
    array([0, 0, 0, 1, 1, 1])
    >>> fcm.predict_proba([[0, 0], [12, 3]])
    array([[0.98484848, 0.01515152],
           [0.00684932, 0.99315068]])
    >>> fcm.cluster_centers_
    array([[ 1.,  2.],
           [10.,  2.]])
    """

    _parameter_constraints: dict = {
        "n_clusters": [Interval(Integral, 1, None, closed="left")],
        "m": [Interval(Real, 1.0, None, closed="neither")],
        "max_iter": [Interval(Integral, 1, None, closed="left")],
        "tol": [Interval(Real, 0.0, None, closed="left")],
        "init": [
            StrOptions({"random_membership", "k-means++", "random_centers"}),
            callable,
            "array-like",
        ],
        "random_state": ["random_state"],
    }

    def __init__(
        self,
        n_clusters=8,
        *,
        m=2.0,
        max_iter=300,
        tol=1e-4,
        init="random_membership",
        random_state=None,
    ):
        self.n_clusters = n_clusters
        self.m = m
        self.max_iter = max_iter
        self.tol = tol
        self.init = init
        self.random_state = random_state

    @_fit_context(prefer_skip_nested_validation=True)
    def fit(self, X, y=None):
        """Compute Fuzzy C-Means clustering.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training instances to cluster.

        y : Ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        # Validate data
        X = validate_data(
            self,
            X,
            accept_sparse=False,
            dtype=[np.float64, np.float32],
            ensure_min_samples=self.n_clusters,
        )

        self.cluster_centers_, self.u_, self.labels_, self.n_iter_ = fuzzy_c_means(
            X,
            n_clusters=self.n_clusters,
            m=self.m,
            max_iter=self.max_iter,
            tol=self.tol,
            init=self.init,
            random_state=self.random_state,
        )
        return self

    def fit_predict(self, X, y=None):
        """Compute cluster centers and predict hard labels for each sample.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            New data to cluster.

        y : Ignored
            Not used, present here for API consistency.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        return self.fit(X).labels_

    def predict(self, X):
        """Predict the hard cluster label for each sample in X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            New data to predict.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample is closest to.
        """
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse=False, reset=False)

        # Compute distances and membership probabilities
        D = euclidean_distances(X, self.cluster_centers_)
        U = _dist_to_membership(D, self.m)
        return np.argmax(U, axis=1)

    def predict_proba(self, X):
        """Predict fuzzy membership probabilities for each sample in X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            New data to predict.

        Returns
        -------
        U : ndarray of shape (n_samples, self.n_clusters)
            Fuzzy membership matrix where each row sums to 1.
        """
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse=False, reset=False)

        D = euclidean_distances(X, self.cluster_centers_)
        return _dist_to_membership(D, self.m)

    def transform(self, X):
        """Transform X to a cluster-distance space.

        In the new space, each dimension is the distance to the cluster centers.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            New data to transform.

        Returns
        -------
        X_new : ndarray of shape (n_samples, self.n_clusters)
            X transformed in the new distance space.
        """
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse=False, reset=False)
        return euclidean_distances(X, self.cluster_centers_)

    def score(self, X, y=None):
        """Opposite of the value of FCM objective function on X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            New data.

        y : Ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        score : float
            Opposite of the FCM objective value on X.
        """
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse=False, reset=False)

        D = euclidean_distances(X, self.cluster_centers_)
        U = _dist_to_membership(D, self.m)
        Um = np.power(U, self.m)

        # Calculate FCM objective sum_{i,j} U_ij^m * d_ij^2
        obj_value = np.sum(Um * (D**2))
        return -obj_value

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Only used to validate feature names with the names seen in ``fit``.
            The output names are always cluster-distance names and are
            unaffected by this value.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Names of the form ``["fuzzycmeans0", ..., "fuzzycmeans{n_clusters-1}"]``,
            matching the ``n_clusters`` columns produced by :meth:`transform`.
        """
        check_is_fitted(self)
        _check_feature_names_in(self, input_features, generate_names=False)

        return np.asarray(
            [f"fuzzycmeans{i}" for i in range(self.n_clusters)],
            dtype=object,
        )

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        # Fuzzy C-Means in this implementation does not support sparse matrices directly
        tags.input_tags.sparse = False
        return tags
