"""Bisecting K-means clustering."""
# Author: Michal Krawczyk <mkrwczyk.1@gmail.com>

import warnings

import numpy as np
import scipy.sparse as sp
from threadpoolctl import threadpool_limits

from ..exceptions import ConvergenceWarning
from ..exceptions import EfficiencyWarning

from ._kmeans import KMeans
from ._kmeans import _kmeans_single_elkan
from ._kmeans import _kmeans_single_lloyd

from ._k_means_common import _inertia_dense
from ._k_means_common import _inertia_sparse

from ._k_means_lloyd import lloyd_iter_chunked_dense
from ._k_means_lloyd import lloyd_iter_chunked_sparse

from ..utils.extmath import row_norms
from ..utils._openmp_helpers import _openmp_effective_n_threads

from ..utils.validation import check_is_fitted
from ..utils.validation import _check_sample_weight
from ..utils.validation import check_random_state


def _check_labels(X, sample_weight, x_squared_norms, centers, n_threads=1):
    """Compute the labels of the given samples and centers.

    Parameters
    ----------
    X : {ndarray, sparse matrix} of shape (n_samples, n_features)
        The input samples to assign to the labels. If sparse matrix, must
        be in CSR format.

    sample_weight : ndarray of shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : ndarray of shape (n_samples,)
        Precomputed squared euclidean norm of each data point, to speed up
        computations.

    centers : ndarray of shape (n_clusters, n_features)
        The cluster centers.

    n_threads : int, default=1
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        The resulting assignment (Labels of each point).
    """
    n_samples = X.shape[0]
    n_clusters = centers.shape[0]

    labels = np.full(n_samples, -1, dtype=np.int32)
    weight_in_clusters = np.zeros(n_clusters, dtype=centers.dtype)
    center_shift = np.zeros_like(weight_in_clusters)

    if sp.issparse(X):
        _labels = lloyd_iter_chunked_sparse
    else:
        _labels = lloyd_iter_chunked_dense

    _labels(
        X,
        sample_weight,
        x_squared_norms,
        centers,
        centers,
        weight_in_clusters,
        labels,
        center_shift,
        n_threads,
        update_centers=False,
    )

    return labels


def _check_labels_threadpool_limit(
    X, sample_weight, x_squared_norms, centers, n_threads=1
):
    """Same as _check_labels but in a threadpool_limits context."""
    with threadpool_limits(limits=1, user_api="blas"):
        labels = _check_labels(X, sample_weight, x_squared_norms, centers, n_threads)

    return labels


class BisectKMeans(KMeans):
    """Bisecting K-Means clustering.

    Read more in the :ref:`User Guide <bisect_k_means>`.

    .. versionadded:: 1.1

    Parameters
    ----------
    n_clusters : int, default=8
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'k-means++', 'random'} or callable, default='k-means++'
        Method for initialization:

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose `n_clusters` observations (rows) at random from data
        for the initial centroids.

        If a callable is passed, it should take arguments X, n_clusters and a
        random state and return an initialization.

    n_init : int, default=10
        Number of time the inner k-means algorithm will be run with different
        centroid seeds in each bisection.
        That will result producing for each bisection best output of n_init
        consecutive runs in terms of inertia.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    max_iter : int, default=30
        Maximum number of iterations of the inner k-means algorithm at each
        bisection.

    verbose : int, default=0
        Verbosity mode.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations  to declare
        convergence. Used in inner k-means algorithm at each bisection to pick
        best possible clusters.

    copy_x : bool, default=True
        When pre-computing distances it is more numerically accurate to center
        the data first. If copy_x is True (default), then the original data is
        not modified. If False, the original data is modified, and put back
        before the function returns, but small numerical differences may be
        introduced by subtracting and then adding the data mean. Note that if
        the original data is not C-contiguous, a copy will be made even if
        copy_x is False. If the original data is sparse, but not in CSR format,
        a copy will be made even if copy_x is False.

    algorithm : {"auto", "full", "elkan"}, default="auto"
        K-means algorithm to use. The classical EM-style algorithm is "full".
        The "elkan" variation is more efficient on data with well-defined
        clusters, by using the triangle inequality. However it's more memory
        intensive due to the allocation of an extra array of shape
        (n_samples, n_clusters).
        For now "auto" (kept for backward compatibiliy) chooses "elkan" but it
        might change in the future for a better heuristic.

    bisect_strategy : {"biggest_sse", "largest_cluster"}, default="biggest_sse"
        Defines how should bisection by performed:

        - "biggest_sse" means that Bisect K-Means will always check
            all calculated cluster for cluster with biggest SSE
            (Sum of squared errors) and bisect it. This approach concentrates on
            precision, but may be costly in terms of execution time (especially for
            larger ammount of data points).

        - "largest_cluster" - Bisect K-Means will always split cluster with
            largest amount of points assigned to it from all clusters
            previously calculated. That should work faster than picking by SSE
            ('biggest_sse') and may produce similar results in most cases.

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point.

    inertia_ : float
        Sum of squared distances of samples to their closest cluster center,
        weighted by the sample weights if provided.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    KMeans : Original implementation of K-Means algorithm.

    Notes
    -----
    Bisection cannot be performed if n_cluster < 2.
    Despite that in case when n_cluster == 1 - centroid will be created and points
    will be assigned to it.

    Also it might be inefficient when n_cluster is equal to 2, due to unnecassary
    calculations for that case.

    Examples
    --------
    >>> from sklearn.cluster import BisectKMeans
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0],
    ...               [10, 6], [10, 8], [10, 10]])
    >>> bisect_means = BisectKMeans(n_clusters=3, random_state=0).fit(X)
    >>> bisect_means.labels_
    array([0, 0, 0, 1, 1, 1, 2, 2, 2], dtype=int32)
    >>> bisect_means.predict([[0, 0], [12, 3]])
    array([0, 1], dtype=int32)
    >>> bisect_means.cluster_centers_
    array([[ 1.,  2.],
           [10.,  2.],
           [10.,  8.]])
    """

    def __init__(
        self,
        n_clusters=8,
        init="k-means++",
        n_init=10,
        random_state=None,
        max_iter=30,
        verbose=0,
        tol=1e-4,
        copy_x=True,
        algorithm="auto",
        bisect_strategy="biggest_sse",
    ):

        super().__init__(
            n_clusters=n_clusters,
            init=init,
            max_iter=max_iter,
            verbose=verbose,
            random_state=random_state,
            tol=tol,
            n_init=n_init,
            copy_x=copy_x,
            algorithm=algorithm,
        )

        self.bisect_strategy = bisect_strategy

    def _compute_bisect_errors(self, X, centers, labels, sample_weight):
        """
        Calculate the sum of squared errors (inertia) and group them by label (center).

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        centers : ndarray of shape (n_clusters, n_features)
            The cluster centers.

        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X.

        Returns
        -------
        errors_by_label : dict
            dictionary containing squared error of each point by label
            as ndarray.
        """
        errors_by_label = {}

        _inertia = _inertia_sparse if sp.issparse(X) else _inertia_dense

        for value in range(centers.shape[0]):
            indexes = labels == value

            data = X[indexes]
            weights = sample_weight[indexes]
            center = centers[value][np.newaxis, :]
            label = np.zeros(data.shape[0], dtype=np.intc)

            errors_by_label[value] = _inertia(
                data, weights, center, label, self._n_threads
            )

        return errors_by_label

    def _check_params(self, X):
        super()._check_params(X)

        # bisect_strategy
        if self.bisect_strategy not in ["biggest_sse", "largest_cluster"]:
            raise ValueError(
                "Bisect Strategy must be 'biggest_sse', "
                "or 'largest_cluster' "
                f"got {self.bisect_strategy} instead."
            )

        # Regular K-Means should do less computations when there are only
        # less than 3 clusters
        if self.n_clusters < 3:
            warnings.warn(
                "BisectKMeans might be inefficient for n_cluster "
                "smaller than 3  "
                "- Use Normal KMeans from sklearn.cluster instead.",
                EfficiencyWarning,
            )

        if X.shape[0] <= 1:
            raise ValueError(
                "BisectKMeans needs more than one sample to perform bisection."
            )

        if hasattr(self.init, "__array__"):
            raise ValueError("BisectKMeans does not support init as array.")

    def _bisect(self, X, sample_weight, random_state):
        """Bisection of data
        Attempts to get best bisection of data by performing regular K-Means
        for different pairs of centroids

        .. note:: Number of attempts is specified by self.n_init

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

            Training instances to cluster.

            .. note:: The data will be converted to C ordering,
                which will cause a memory copy
                if the given data is not C-contiguous.

        sample_weight : array-like of shape (n_samples,)
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        random_state : int, RandomState instance
            Determines random number generation for centroid initialization.

        Returns
        -------
        self
            Fitted estimator.
        """
        init = self.init
        x_squared_norms = row_norms(X, squared=True)

        best_inertia = None

        for i in range(self.n_init):
            centers_init = self._init_centroids(
                X, x_squared_norms, init, random_state, n_centroids=2
            )

            labels, inertia, centers, _ = self._kmeans_single(
                X,
                sample_weight,
                centers_init,
                max_iter=self.max_iter,
                verbose=self.verbose,
                tol=self.tol,
                x_squared_norms=x_squared_norms,
                n_threads=self._n_threads,
            )

            # allow small tolerance on the inertia to accommodate for
            # non-deterministic rounding errors due to parallel computation
            if best_inertia is None or inertia < best_inertia * (1 - 1e-6):
                best_labels = labels
                best_centers = centers
                best_inertia = inertia

        distinct_clusters = len(set(best_labels))
        if distinct_clusters != 2:
            warnings.warn(
                "Number of distinct clusters ({}) found smaller than "
                "n_clusters ({}). Possibly due to duplicate points "
                "in X.".format(distinct_clusters, 2),
                ConvergenceWarning,
                stacklevel=2,
            )

        return best_centers, best_labels

    def fit(self, X, y=None, sample_weight=None):
        """Compute bisecting k-means clustering.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

            Training instances to cluster.

            .. note:: The data will be converted to C ordering,
                which will cause a memory copy
                if the given data is not C-contiguous.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        self
            Fitted estimator.
        """
        X = self._validate_data(
            X,
            accept_sparse="csr",
            dtype=[np.float64, np.float32],
            order="C",
            copy=self.copy_x,
            accept_large_sparse=False,
        )

        self._check_params(X)
        random_state = check_random_state(self.random_state)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)
        self._n_threads = _openmp_effective_n_threads()

        if self._algorithm == "full":
            self._kmeans_single = _kmeans_single_lloyd
            self._check_mkl_vcomp(X, X.shape[0])
        else:
            self._kmeans_single = _kmeans_single_elkan

        _inertia = _inertia_sparse if sp.issparse(X) else _inertia_dense

        if self.verbose:
            print("Running Bisecting K-Means with parameters:")
            print(f"-> number of clusters: {self.n_clusters}")
            print(f"-> number of centroid initializations: {self.n_init}")
            print(f"-> relative tolerance: {self.tol:.4e}")
            print(f"-> bisect strategy: {self.bisect_strategy} \n")

        # Subtract of mean of X for more accurate distance computations
        if not sp.issparse(X):
            self._X_mean = X.mean(axis=0)
            X -= self._X_mean

        # Only assign to created centroid when n_clusters == 1
        if self.n_clusters == 1:
            x_squared_norms = row_norms(X, squared=True)

            clusters = self._init_centroids(
                X, x_squared_norms, self.init, random_state, n_centroids=1
            )
            warnings.warn(
                "Bisection won't be performed - needs at least two clusters to run."
            )

            self.labels_ = np.zeros(X.shape[0], dtype=np.intc)

        else:
            # Run proper bisection to gather
            # self.cluster_centers_ and self.labels_
            clusters, self.labels_ = self._run_bisect_kmeans(
                X, sample_weight, random_state
            )

        # Restore original data
        if not sp.issparse(X):
            X += self._X_mean
            clusters += self._X_mean

        self.cluster_centers_ = clusters

        self.inertia_ = _inertia(
            X, sample_weight, self.cluster_centers_, self.labels_, self._n_threads
        )

        return self

    def _run_bisect_kmeans(self, X, sample_weight, random_state):
        """Performs Bisecting K-Means, which splits always cluster depending
        on 'bisect_strategy' attribute:

        - "biggest sse": Picks cluster with biggest SSE (Sum of Squared Errors)
         from all calculated.

        - "largest_cluster": Picks always cluster with largest number of
         points assigned from all calculated. That method will perform faster
          than picking by SSE methods, while producing similar results.

        .. note:: All of passed parameters must be pre-calculated

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training instances to cluster.

        sample_weight : array-like of shape (n_samples,)
            The weights for each observation in X.

        random_state : int, RandomState instance
            Determines random number generation for centroid initialization.

        Returns
        ----------
        cluster_centers_ : ndarray of shape (n_clusters, n_features)
            Coordinates of cluster centers.

        labels_ : ndarray of shape (n_samples,)
            Labels of each point.
        """
        strategy = self.bisect_strategy

        # Dictionary to imitate tree view of created centers.
        # Used to keep hierarchical ordering.
        tree_dict = {-1: {"children": None, "center": None}}

        # Leaves Dictionary used for clustering.
        # Stores information which data points are assigned to given cluster
        # along with error or amount of assigned points, depending on
        # specified 'bisect_strategy'
        leaves_dict = {
            -1: {"samples": np.ones(X.shape[0], dtype=bool), "error_or_size": None}
        }

        # ID of biggest center stored in centers_dict
        parent_id = -1

        last_center_id = 0

        for n_iter in range(self.n_clusters - 1):
            picked_samples = leaves_dict[parent_id]["samples"]

            # Pick data and weights to bisect
            picked_data = X[picked_samples]
            picked_weights = sample_weight[picked_samples]

            # Perform Bisection
            _centers, _labels = self._bisect(picked_data, picked_weights, random_state)

            if self.verbose:
                print(f"Centroid Found: {_centers[0]}")
                print(f"Centroid Found: {_centers[1]}")

            if strategy == "biggest_sse":
                # "biggest_sse" uses computed SSE (Sum of Squared Errors)
                # to pick cluster with biggest SSE Error
                metrics_values = self._compute_bisect_errors(
                    picked_data, _centers, _labels, picked_weights
                )
            else:
                # strategy == 'largest_cluster'
                # 'largest_cluster' takes counts occurances of each label
                # to pick cluster with largest number of data points assigned
                metrics_values = np.bincount(_labels)

            # "Create Hierarchy":
            # Cluster with smaller metrics value (SSE or ammount of points)
            # will be at 'left side' and cluster with higher at 'right side'
            ordered_labels = (
                (0, 1) if metrics_values[0] <= metrics_values[1] else (1, 0)
            )

            # Assign calculated nested clusters to their root
            tree_dict[parent_id]["children"] = [last_center_id, last_center_id + 1]

            for ii, label in enumerate(ordered_labels):
                child_id = tree_dict[parent_id]["children"][ii]

                # Save Results on Tree
                tree_dict[child_id] = {"children": None, "center": _centers[label]}

                # Create Mask for samples for selecting proper data points
                samples_mask = leaves_dict[parent_id]["samples"].copy()
                samples_mask[picked_samples] = _labels == label

                # Save recently generated leaves
                leaves_dict[child_id] = {
                    "samples": samples_mask,
                    "error_or_size": metrics_values[label],
                }

            # Split cluster is no longer leaf
            del leaves_dict[parent_id]

            last_center_id += 2

            # Pick new 'biggest cluster' to split
            parent_id, _ = max(leaves_dict.items(), key=lambda x: x[1]["error_or_size"])

        # Delete Initial cluster
        del tree_dict[-1]

        # Initialize Labels
        labels = np.full(X.shape[0], -1, dtype=np.intc)

        centers = []

        for i, (leaf_id, leaf) in enumerate(leaves_dict.items()):
            # Save leaf cluster centers
            centers.append(tree_dict[leaf_id]["center"])

            # Assign labels to proper data points
            labels[leaf["samples"]] = i

        centers = np.asarray(centers)

        # Inner Tree will be later used at 'predict' method:
        self._inner_tree = tree_dict

        return centers, labels

    def predict(self, X, sample_weight=None):
        """Predict which cluster each sample in X belongs to.

        Prediction is made by going down the hierarchical tree
        in searching of closest leaf cluster.

        In the vector quantization literature, `cluster_centers_` is called
        the code book and each value returned by `predict` is the index of
        the closest code in the code book.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to predict.

        sample_weight : Ignored
            Not used, present here for API consistency by convention.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        x_squared_norms = row_norms(X, squared=True)

        # sample_weight is needed for labeling in _check_labels_threadpool_limit()
        # as argument, but not used there since its not updating centers.
        # In this case, it doesn't have to be shape of X - only ndarray with same dtype
        sample_weight = np.zeros(1, dtype=X.dtype)

        # With only one cluster all points have same label
        if self.n_clusters == 1:
            return np.zeros(X.shape[0], dtype=np.intc)

        # Get all clusters from inner tree as array
        all_centers = np.array([val["center"] for val in self._inner_tree.values()])

        # Centers stored in tree have subtracted mean for calculating centroids
        if not sp.issparse(X):
            all_centers += self._X_mean

        # Init labels for first two clusters
        centers_hierarchical = [0, 1]
        labels = _check_labels_threadpool_limit(
            X, sample_weight, x_squared_norms, all_centers[:2], self._n_threads
        )

        # Go down the tree with samples to assign them to proper leaves
        while len(centers_hierarchical) != self.n_clusters:
            nested_centers = []

            for idx in centers_hierarchical:
                center_ids = [idx]

                # If center has children - relabel its data to child centers
                if self._inner_tree[idx]["children"] is not None:
                    center_ids = self._inner_tree[idx]["children"]
                    picked_samples = labels == idx

                    new_labels = _check_labels_threadpool_limit(
                        X[picked_samples],
                        sample_weight,
                        x_squared_norms[picked_samples],
                        all_centers[center_ids[0] : (center_ids[1] + 1)],
                        self._n_threads,
                    )

                    # Label with id of center
                    new_labels[new_labels == 0] = center_ids[0]
                    new_labels[new_labels == 1] = center_ids[1]

                    # Move new assignment to predicted labels
                    labels[picked_samples] = new_labels

                nested_centers.extend(center_ids)
            centers_hierarchical = nested_centers

        # Get all ids of leaves for proper assignment
        leaves_ids = [
            x[0] for x in self._inner_tree.items() if x[1]["children"] is None
        ]

        # Copy of labels will be used to make sure that samples are picked correctly
        temp_labels = labels.copy()

        for i, center_id in enumerate(leaves_ids):
            # Assign labels to proper data points
            labels[temp_labels == center_id] = i

        return labels
