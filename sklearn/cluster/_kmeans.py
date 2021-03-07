"""K-means clustering."""

# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Thomas Rueckstiess <ruecksti@in.tum.de>
#          James Bergstra <james.bergstra@umontreal.ca>
#          Jan Schlueter <scikit-learn@jan-schlueter.de>
#          Nelle Varoquaux
#          Peter Prettenhofer <peter.prettenhofer@gmail.com>
#          Olivier Grisel <olivier.grisel@ensta.org>
#          Mathieu Blondel <mathieu@mblondel.org>
#          Robert Layton <robertlayton@gmail.com>
# License: BSD 3 clause

import warnings

import numpy as np
import scipy.sparse as sp
from threadpoolctl import threadpool_limits
from threadpoolctl import threadpool_info

from ..base import BaseEstimator, ClusterMixin, TransformerMixin
from ..metrics.pairwise import euclidean_distances
from ..utils.extmath import row_norms, stable_cumsum
from ..utils.sparsefuncs_fast import assign_rows_csr
from ..utils.sparsefuncs import mean_variance_axis
from ..utils.validation import _deprecate_positional_args
from ..utils import check_array
from ..utils import gen_batches
from ..utils import check_random_state
from ..utils import deprecated
from ..utils.validation import check_is_fitted, _check_sample_weight
from ..utils._openmp_helpers import _openmp_effective_n_threads
from ..exceptions import ConvergenceWarning
from ._k_means_fast import CHUNK_SIZE
from ._k_means_fast import _inertia_dense
from ._k_means_fast import _inertia_sparse
from ._k_means_fast import _mini_batch_update_csr
from ._k_means_lloyd import lloyd_iter_chunked_dense
from ._k_means_lloyd import lloyd_iter_chunked_sparse
from ._k_means_elkan import init_bounds_dense
from ._k_means_elkan import init_bounds_sparse
from ._k_means_elkan import elkan_iter_chunked_dense
from ._k_means_elkan import elkan_iter_chunked_sparse


###############################################################################
# Initialization heuristic


def _kmeans_plusplus(X, n_clusters, x_squared_norms,
                     random_state, n_local_trials=None):
    """Computational component for initialization of n_clusters by
    k-means++. Prior validation of data is assumed.

    Parameters
    ----------
    X : {ndarray, sparse matrix} of shape (n_samples, n_features)
        The data to pick seeds for.

    n_clusters : int
        The number of seeds to choose.

    x_squared_norms : ndarray of shape (n_samples,)
        Squared Euclidean norm of each data point.

    random_state : RandomState instance
        The generator used to initialize the centers.
        See :term:`Glossary <random_state>`.

    n_local_trials : int, default=None
        The number of seeding trials for each center (except the first),
        of which the one reducing inertia the most is greedily chosen.
        Set to None to make the number of trials depend logarithmically
        on the number of seeds (2+log(k)); this is the default.

    Returns
    -------
    centers : ndarray of shape (n_clusters, n_features)
        The inital centers for k-means.

    indices : ndarray of shape (n_clusters,)
        The index location of the chosen centers in the data array X. For a
        given index and center, X[index] = center.
    """
    n_samples, n_features = X.shape

    centers = np.empty((n_clusters, n_features), dtype=X.dtype)

    # Set the number of local seeding trials if none is given
    if n_local_trials is None:
        # This is what Arthur/Vassilvitskii tried, but did not report
        # specific results for other than mentioning in the conclusion
        # that it helped.
        n_local_trials = 2 + int(np.log(n_clusters))

    # Pick first center randomly and track index of point
    center_id = random_state.randint(n_samples)
    indices = np.full(n_clusters, -1, dtype=int)
    if sp.issparse(X):
        centers[0] = X[center_id].toarray()
    else:
        centers[0] = X[center_id]
    indices[0] = center_id

    # Initialize list of closest distances and calculate current potential
    closest_dist_sq = euclidean_distances(
        centers[0, np.newaxis], X, Y_norm_squared=x_squared_norms,
        squared=True)
    current_pot = closest_dist_sq.sum()

    # Pick the remaining n_clusters-1 points
    for c in range(1, n_clusters):
        # Choose center candidates by sampling with probability proportional
        # to the squared distance to the closest existing center
        rand_vals = random_state.random_sample(n_local_trials) * current_pot
        candidate_ids = np.searchsorted(stable_cumsum(closest_dist_sq),
                                        rand_vals)
        # XXX: numerical imprecision can result in a candidate_id out of range
        np.clip(candidate_ids, None, closest_dist_sq.size - 1,
                out=candidate_ids)

        # Compute distances to center candidates
        distance_to_candidates = euclidean_distances(
            X[candidate_ids], X, Y_norm_squared=x_squared_norms, squared=True)

        # update closest distances squared and potential for each candidate
        np.minimum(closest_dist_sq, distance_to_candidates,
                   out=distance_to_candidates)
        candidates_pot = distance_to_candidates.sum(axis=1)

        # Decide which candidate is the best
        best_candidate = np.argmin(candidates_pot)
        current_pot = candidates_pot[best_candidate]
        closest_dist_sq = distance_to_candidates[best_candidate]
        best_candidate = candidate_ids[best_candidate]

        # Permanently add best center candidate found in local tries
        if sp.issparse(X):
            centers[c] = X[best_candidate].toarray()
        else:
            centers[c] = X[best_candidate]
        indices[c] = best_candidate

    return centers, indices


###############################################################################
# K-means batch estimation by EM (expectation maximization)

def _tolerance(X, tol):
    """Return a tolerance which is independent of the dataset."""
    if tol == 0:
        return 0
    if sp.issparse(X):
        variances = mean_variance_axis(X, axis=0)[1]
    else:
        variances = np.var(X, axis=0)
    return np.mean(variances) * tol


@_deprecate_positional_args
def k_means(X, n_clusters, *, sample_weight=None, init='k-means++',
            precompute_distances='deprecated', n_init=10, max_iter=300,
            verbose=False, tol=1e-4, random_state=None, copy_x=True,
            n_jobs='deprecated', algorithm="auto", return_n_iter=False):
    """K-means clustering algorithm.

    Read more in the :ref:`User Guide <k_means>`.

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        The observations to cluster. It must be noted that the data
        will be converted to C ordering, which will cause a memory copy
        if the given data is not C-contiguous.

    n_clusters : int
        The number of clusters to form as well as the number of
        centroids to generate.

    sample_weight : array-like of shape (n_samples,), default=None
        The weights for each observation in X. If None, all observations
        are assigned equal weight.

    init : {'k-means++', 'random'}, callable or array-like of shape \
            (n_clusters, n_features), default='k-means++'
        Method for initialization:

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose `n_clusters` observations (rows) at random from data
        for the initial centroids.

        If an array is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

        If a callable is passed, it should take arguments X, n_clusters and a
        random state and return an initialization.

    precompute_distances : {'auto', True, False}
        Precompute distances (faster but takes more memory).

        'auto' : do not precompute distances if n_samples * n_clusters > 12
        million. This corresponds to about 100MB overhead per job using
        double precision.

        True : always precompute distances

        False : never precompute distances

        .. deprecated:: 0.23
            'precompute_distances' was deprecated in version 0.23 and will be
            removed in 1.0 (renaming of 0.25). It has no effect.

    n_init : int, default=10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm to run.

    verbose : bool, default=False
        Verbosity mode.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations to declare
        convergence.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    copy_x : bool, default=True
        When pre-computing distances it is more numerically accurate to center
        the data first. If copy_x is True (default), then the original data is
        not modified. If False, the original data is modified, and put back
        before the function returns, but small numerical differences may be
        introduced by subtracting and then adding the data mean. Note that if
        the original data is not C-contiguous, a copy will be made even if
        copy_x is False. If the original data is sparse, but not in CSR format,
        a copy will be made even if copy_x is False.

    n_jobs : int, default=None
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

        ``None`` or ``-1`` means using all processors.

        .. deprecated:: 0.23
            ``n_jobs`` was deprecated in version 0.23 and will be removed in
            1.0 (renaming of 0.25).

    algorithm : {"auto", "full", "elkan"}, default="auto"
        K-means algorithm to use. The classical EM-style algorithm is "full".
        The "elkan" variation is more efficient on data with well-defined
        clusters, by using the triangle inequality. However it's more memory
        intensive due to the allocation of an extra array of shape
        (n_samples, n_clusters).

        For now "auto" (kept for backward compatibility) chooses "elkan" but it
        might change in the future for a better heuristic.

    return_n_iter : bool, default=False
        Whether or not to return the number of iterations.

    Returns
    -------
    centroid : ndarray of shape (n_clusters, n_features)
        Centroids found at the last iteration of k-means.

    label : ndarray of shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia : float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).

    best_n_iter : int
        Number of iterations corresponding to the best results.
        Returned only if `return_n_iter` is set to True.
    """
    est = KMeans(
        n_clusters=n_clusters, init=init, n_init=n_init, max_iter=max_iter,
        verbose=verbose, precompute_distances=precompute_distances, tol=tol,
        random_state=random_state, copy_x=copy_x, n_jobs=n_jobs,
        algorithm=algorithm
    ).fit(X, sample_weight=sample_weight)
    if return_n_iter:
        return est.cluster_centers_, est.labels_, est.inertia_, est.n_iter_
    else:
        return est.cluster_centers_, est.labels_, est.inertia_


def _kmeans_single_elkan(X, sample_weight, centers_init, max_iter=300,
                         verbose=False, x_squared_norms=None, tol=1e-4,
                         n_threads=1):
    """A single run of k-means elkan, assumes preparation completed prior.

    Parameters
    ----------
    X : {ndarray, sparse matrix} of shape (n_samples, n_features)
        The observations to cluster. If sparse matrix, must be in CSR format.

    sample_weight : array-like of shape (n_samples,)
        The weights for each observation in X.

    centers_init : ndarray of shape (n_clusters, n_features)
        The initial centers.

    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm to run.

    verbose : bool, default=False
        Verbosity mode.

    x_squared_norms : array-like, default=None
        Precomputed x_squared_norms.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations to declare
        convergence.
        It's not advised to set `tol=0` since convergence might never be
        declared due to rounding errors. Use a very small number instead.

    n_threads : int, default=1
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

    Returns
    -------
    centroid : ndarray of shape (n_clusters, n_features)
        Centroids found at the last iteration of k-means.

    label : ndarray of shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia : float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).

    n_iter : int
        Number of iterations run.
    """
    n_samples = X.shape[0]
    n_clusters = centers_init.shape[0]

    # Buffers to avoid new allocations at each iteration.
    centers = centers_init
    centers_new = np.zeros_like(centers)
    weight_in_clusters = np.zeros(n_clusters, dtype=X.dtype)
    labels = np.full(n_samples, -1, dtype=np.int32)
    labels_old = labels.copy()
    center_half_distances = euclidean_distances(centers) / 2
    distance_next_center = np.partition(np.asarray(center_half_distances),
                                        kth=1, axis=0)[1]
    upper_bounds = np.zeros(n_samples, dtype=X.dtype)
    lower_bounds = np.zeros((n_samples, n_clusters), dtype=X.dtype)
    center_shift = np.zeros(n_clusters, dtype=X.dtype)

    if sp.issparse(X):
        init_bounds = init_bounds_sparse
        elkan_iter = elkan_iter_chunked_sparse
        _inertia = _inertia_sparse
    else:
        init_bounds = init_bounds_dense
        elkan_iter = elkan_iter_chunked_dense
        _inertia = _inertia_dense

    init_bounds(X, centers, center_half_distances,
                labels, upper_bounds, lower_bounds)

    strict_convergence = False

    for i in range(max_iter):
        elkan_iter(X, sample_weight, centers, centers_new,
                   weight_in_clusters, center_half_distances,
                   distance_next_center, upper_bounds, lower_bounds,
                   labels, center_shift, n_threads)

        # compute new pairwise distances between centers and closest other
        # center of each center for next iterations
        center_half_distances = euclidean_distances(centers_new) / 2
        distance_next_center = np.partition(
            np.asarray(center_half_distances), kth=1, axis=0)[1]

        if verbose:
            inertia = _inertia(X, sample_weight, centers, labels)
            print(f"Iteration {i}, inertia {inertia}")

        centers, centers_new = centers_new, centers

        if np.array_equal(labels, labels_old):
            # First check the labels for strict convergence.
            if verbose:
                print(f"Converged at iteration {i}: strict convergence.")
            strict_convergence = True
            break
        else:
            # No strict convergence, check for tol based convergence.
            center_shift_tot = (center_shift**2).sum()
            if center_shift_tot <= tol:
                if verbose:
                    print(f"Converged at iteration {i}: center shift "
                          f"{center_shift_tot} within tolerance {tol}.")
                break

        labels_old[:] = labels

    if not strict_convergence:
        # rerun E-step so that predicted labels match cluster centers
        elkan_iter(X, sample_weight, centers, centers, weight_in_clusters,
                   center_half_distances, distance_next_center,
                   upper_bounds, lower_bounds, labels, center_shift,
                   n_threads, update_centers=False)

    inertia = _inertia(X, sample_weight, centers, labels)

    return labels, inertia, centers, i + 1


def _kmeans_single_lloyd(X, sample_weight, centers_init, max_iter=300,
                         verbose=False, x_squared_norms=None, tol=1e-4,
                         n_threads=1):
    """A single run of k-means lloyd, assumes preparation completed prior.

    Parameters
    ----------
    X : {ndarray, sparse matrix} of shape (n_samples, n_features)
        The observations to cluster. If sparse matrix, must be in CSR format.

    sample_weight : ndarray of shape (n_samples,)
        The weights for each observation in X.

    centers_init : ndarray of shape (n_clusters, n_features)
        The initial centers.

    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm to run.

    verbose : bool, default=False
        Verbosity mode

    x_squared_norms : ndarray of shape (n_samples,), default=None
        Precomputed x_squared_norms.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations to declare
        convergence.
        It's not advised to set `tol=0` since convergence might never be
        declared due to rounding errors. Use a very small number instead.

    n_threads : int, default=1
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

    Returns
    -------
    centroid : ndarray of shape (n_clusters, n_features)
        Centroids found at the last iteration of k-means.

    label : ndarray of shape (n_samples,)
        label[i] is the code or index of the centroid the
        i'th observation is closest to.

    inertia : float
        The final value of the inertia criterion (sum of squared distances to
        the closest centroid for all observations in the training set).

    n_iter : int
        Number of iterations run.
    """
    n_clusters = centers_init.shape[0]

    # Buffers to avoid new allocations at each iteration.
    centers = centers_init
    centers_new = np.zeros_like(centers)
    labels = np.full(X.shape[0], -1, dtype=np.int32)
    labels_old = labels.copy()
    weight_in_clusters = np.zeros(n_clusters, dtype=X.dtype)
    center_shift = np.zeros(n_clusters, dtype=X.dtype)

    if sp.issparse(X):
        lloyd_iter = lloyd_iter_chunked_sparse
        _inertia = _inertia_sparse
    else:
        lloyd_iter = lloyd_iter_chunked_dense
        _inertia = _inertia_dense

    strict_convergence = False

    # Threadpoolctl context to limit the number of threads in second level of
    # nested parallelism (i.e. BLAS) to avoid oversubsciption.
    with threadpool_limits(limits=1, user_api="blas"):
        for i in range(max_iter):
            lloyd_iter(X, sample_weight, x_squared_norms, centers, centers_new,
                       weight_in_clusters, labels, center_shift, n_threads)

            if verbose:
                inertia = _inertia(X, sample_weight, centers, labels)
                print(f"Iteration {i}, inertia {inertia}.")

            centers, centers_new = centers_new, centers

            if np.array_equal(labels, labels_old):
                # First check the labels for strict convergence.
                if verbose:
                    print(f"Converged at iteration {i}: strict convergence.")
                strict_convergence = True
                break
            else:
                # No strict convergence, check for tol based convergence.
                center_shift_tot = (center_shift**2).sum()
                if center_shift_tot <= tol:
                    if verbose:
                        print(f"Converged at iteration {i}: center shift "
                              f"{center_shift_tot} within tolerance {tol}.")
                    break

            labels_old[:] = labels

        if not strict_convergence:
            # rerun E-step so that predicted labels match cluster centers
            lloyd_iter(X, sample_weight, x_squared_norms, centers, centers,
                       weight_in_clusters, labels, center_shift, n_threads,
                       update_centers=False)

    inertia = _inertia(X, sample_weight, centers, labels)

    return labels, inertia, centers, i + 1


def _labels_inertia(X, sample_weight, x_squared_norms, centers,
                    n_threads=None):
    """E step of the K-means EM algorithm.

    Compute the labels and the inertia of the given samples and centers.

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

    n_threads : int, default=None
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

    Returns
    -------
    labels : ndarray of shape (n_samples,)
        The resulting assignment.

    inertia : float
        Sum of squared distances of samples to their closest cluster center.
    """
    n_samples = X.shape[0]
    n_clusters = centers.shape[0]

    n_threads = _openmp_effective_n_threads(n_threads)

    labels = np.full(n_samples, -1, dtype=np.int32)
    weight_in_clusters = np.zeros(n_clusters, dtype=centers.dtype)
    center_shift = np.zeros_like(weight_in_clusters)

    if sp.issparse(X):
        _labels = lloyd_iter_chunked_sparse
        _inertia = _inertia_sparse
    else:
        _labels = lloyd_iter_chunked_dense
        _inertia = _inertia_dense

    _labels(X, sample_weight, x_squared_norms, centers, centers,
            weight_in_clusters, labels, center_shift, n_threads,
            update_centers=False)

    inertia = _inertia(X, sample_weight, centers, labels)

    return labels, inertia


class KMeans(TransformerMixin, ClusterMixin, BaseEstimator):
    """K-Means clustering.

    Read more in the :ref:`User Guide <k_means>`.

    Parameters
    ----------

    n_clusters : int, default=8
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'k-means++', 'random'}, callable or array-like of shape \
            (n_clusters, n_features), default='k-means++'
        Method for initialization:

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose `n_clusters` observations (rows) at random from data
        for the initial centroids.

        If an array is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

        If a callable is passed, it should take arguments X, n_clusters and a
        random state and return an initialization.

    n_init : int, default=10
        Number of time the k-means algorithm will be run with different
        centroid seeds. The final results will be the best output of
        n_init consecutive runs in terms of inertia.

    max_iter : int, default=300
        Maximum number of iterations of the k-means algorithm for a
        single run.

    tol : float, default=1e-4
        Relative tolerance with regards to Frobenius norm of the difference
        in the cluster centers of two consecutive iterations to declare
        convergence.

    precompute_distances : {'auto', True, False}, default='auto'
        Precompute distances (faster but takes more memory).

        'auto' : do not precompute distances if n_samples * n_clusters > 12
        million. This corresponds to about 100MB overhead per job using
        double precision.

        True : always precompute distances.

        False : never precompute distances.

        .. deprecated:: 0.23
            'precompute_distances' was deprecated in version 0.22 and will be
            removed in 1.0 (renaming of 0.25). It has no effect.

    verbose : int, default=0
        Verbosity mode.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    copy_x : bool, default=True
        When pre-computing distances it is more numerically accurate to center
        the data first. If copy_x is True (default), then the original data is
        not modified. If False, the original data is modified, and put back
        before the function returns, but small numerical differences may be
        introduced by subtracting and then adding the data mean. Note that if
        the original data is not C-contiguous, a copy will be made even if
        copy_x is False. If the original data is sparse, but not in CSR format,
        a copy will be made even if copy_x is False.

    n_jobs : int, default=None
        The number of OpenMP threads to use for the computation. Parallelism is
        sample-wise on the main cython loop which assigns each sample to its
        closest center.

        ``None`` or ``-1`` means using all processors.

        .. deprecated:: 0.23
            ``n_jobs`` was deprecated in version 0.23 and will be removed in
            1.0 (renaming of 0.25).

    algorithm : {"auto", "full", "elkan"}, default="auto"
        K-means algorithm to use. The classical EM-style algorithm is "full".
        The "elkan" variation is more efficient on data with well-defined
        clusters, by using the triangle inequality. However it's more memory
        intensive due to the allocation of an extra array of shape
        (n_samples, n_clusters).

        For now "auto" (kept for backward compatibiliy) chooses "elkan" but it
        might change in the future for a better heuristic.

        .. versionchanged:: 0.18
            Added Elkan algorithm

    Attributes
    ----------
    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers. If the algorithm stops before fully
        converging (see ``tol`` and ``max_iter``), these will not be
        consistent with ``labels_``.

    labels_ : ndarray of shape (n_samples,)
        Labels of each point

    inertia_ : float
        Sum of squared distances of samples to their closest cluster center.

    n_iter_ : int
        Number of iterations run.

    See Also
    --------
    MiniBatchKMeans : Alternative online implementation that does incremental
        updates of the centers positions using mini-batches.
        For large scale learning (say n_samples > 10k) MiniBatchKMeans is
        probably much faster than the default batch implementation.

    Notes
    -----
    The k-means problem is solved using either Lloyd's or Elkan's algorithm.

    The average complexity is given by O(k n T), where n is the number of
    samples and T is the number of iteration.

    The worst case complexity is given by O(n^(k+2/p)) with
    n = n_samples, p = n_features. (D. Arthur and S. Vassilvitskii,
    'How slow is the k-means method?' SoCG2006)

    In practice, the k-means algorithm is very fast (one of the fastest
    clustering algorithms available), but it falls in local minima. That's why
    it can be useful to restart it several times.

    If the algorithm stops before fully converging (because of ``tol`` or
    ``max_iter``), ``labels_`` and ``cluster_centers_`` will not be consistent,
    i.e. the ``cluster_centers_`` will not be the means of the points in each
    cluster. Also, the estimator will reassign ``labels_`` after the last
    iteration to make ``labels_`` consistent with ``predict`` on the training
    set.

    Examples
    --------

    >>> from sklearn.cluster import KMeans
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0]])
    >>> kmeans = KMeans(n_clusters=2, random_state=0).fit(X)
    >>> kmeans.labels_
    array([1, 1, 1, 0, 0, 0], dtype=int32)
    >>> kmeans.predict([[0, 0], [12, 3]])
    array([1, 0], dtype=int32)
    >>> kmeans.cluster_centers_
    array([[10.,  2.],
           [ 1.,  2.]])
    """
    @_deprecate_positional_args
    def __init__(self, n_clusters=8, *, init='k-means++', n_init=10,
                 max_iter=300, tol=1e-4, precompute_distances='deprecated',
                 verbose=0, random_state=None, copy_x=True,
                 n_jobs='deprecated', algorithm='auto'):

        self.n_clusters = n_clusters
        self.init = init
        self.max_iter = max_iter
        self.tol = tol
        self.precompute_distances = precompute_distances
        self.n_init = n_init
        self.verbose = verbose
        self.random_state = random_state
        self.copy_x = copy_x
        self.n_jobs = n_jobs
        self.algorithm = algorithm

    def _check_params(self, X):
        # precompute_distances
        if self.precompute_distances != 'deprecated':
            warnings.warn("'precompute_distances' was deprecated in version "
                          "0.23 and will be removed in 1.0 (renaming of 0.25)"
                          ". It has no effect", FutureWarning)

        # n_jobs
        if self.n_jobs != 'deprecated':
            warnings.warn("'n_jobs' was deprecated in version 0.23 and will be"
                          " removed in 1.0 (renaming of 0.25).", FutureWarning)
            self._n_threads = self.n_jobs
        else:
            self._n_threads = None
        self._n_threads = _openmp_effective_n_threads(self._n_threads)

        # n_init
        if self.n_init <= 0:
            raise ValueError(
                f"n_init should be > 0, got {self.n_init} instead.")
        self._n_init = self.n_init

        # max_iter
        if self.max_iter <= 0:
            raise ValueError(
                f"max_iter should be > 0, got {self.max_iter} instead.")

        # n_clusters
        if X.shape[0] < self.n_clusters:
            raise ValueError(f"n_samples={X.shape[0]} should be >= "
                             f"n_clusters={self.n_clusters}.")

        # tol
        self._tol = _tolerance(X, self.tol)

        # algorithm
        if self.algorithm not in ("auto", "full", "elkan"):
            raise ValueError(f"Algorithm must be 'auto', 'full' or 'elkan', "
                             f"got {self.algorithm} instead.")

        self._algorithm = self.algorithm
        if self._algorithm == "auto":
            self._algorithm = "full" if self.n_clusters == 1 else "elkan"
        if self._algorithm == "elkan" and self.n_clusters == 1:
            warnings.warn("algorithm='elkan' doesn't make sense for a single "
                          "cluster. Using 'full' instead.", RuntimeWarning)
            self._algorithm = "full"

        # init
        if not (hasattr(self.init, '__array__') or callable(self.init)
                or (isinstance(self.init, str)
                    and self.init in ["k-means++", "random"])):
            raise ValueError(
                f"init should be either 'k-means++', 'random', a ndarray or a "
                f"callable, got '{self.init}' instead.")

        if hasattr(self.init, '__array__') and self._n_init != 1:
            warnings.warn(
                f"Explicit initial center position passed: performing only"
                f" one init in {self.__class__.__name__} instead of "
                f"n_init={self._n_init}.", RuntimeWarning, stacklevel=2)
            self._n_init = 1

    def _validate_center_shape(self, X, centers):
        """Check if centers is compatible with X and n_clusters."""
        if centers.shape[0] != self.n_clusters:
            raise ValueError(
                f"The shape of the initial centers {centers.shape} does not "
                f"match the number of clusters {self.n_clusters}.")
        if centers.shape[1] != X.shape[1]:
            raise ValueError(
                f"The shape of the initial centers {centers.shape} does not "
                f"match the number of features of the data {X.shape[1]}.")

    def _check_test_data(self, X):
        X = self._validate_data(X, accept_sparse='csr', reset=False,
                                dtype=[np.float64, np.float32],
                                order='C', accept_large_sparse=False)
        return X

    def _check_mkl_vcomp(self, X, n_samples):
        """Warns when vcomp and mkl are both present"""
        # The BLAS call inside a prange in lloyd_iter_chunked_dense is known to
        # cause a small memory leak when there are less chunks than the number
        # of available threads. It only happens when the OpenMP library is
        # vcomp (microsoft OpenMP) and the BLAS library is MKL. see #18653
        if sp.issparse(X):
            return

        active_threads = int(np.ceil(n_samples / CHUNK_SIZE))
        if active_threads < self._n_threads:
            modules = threadpool_info()
            has_vcomp = "vcomp" in [module["prefix"] for module in modules]
            has_mkl = ("mkl", "intel") in [
                (module["internal_api"], module.get("threading_layer", None))
                for module in modules]
            if has_vcomp and has_mkl:
                if not hasattr(self, "batch_size"):  # KMeans
                    warnings.warn(
                        f"KMeans is known to have a memory leak on Windows "
                        f"with MKL, when there are less chunks than available "
                        f"threads. You can avoid it by setting the environment"
                        f" variable OMP_NUM_THREADS={active_threads}.")
                else:  # MiniBatchKMeans
                    warnings.warn(
                        f"MiniBatchKMeans is known to have a memory leak on "
                        f"Windows with MKL, when there are less chunks than "
                        f"available threads. You can prevent it by setting "
                        f"batch_size >= {self._n_threads * CHUNK_SIZE} or by "
                        f"setting the environment variable "
                        f"OMP_NUM_THREADS={active_threads}")

    def _init_centroids(self, X, x_squared_norms, init, random_state,
                        init_size=None):
        """Compute the initial centroids.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        x_squared_norms : ndarray of shape (n_samples,)
            Squared euclidean norm of each data point. Pass it if you have it
            at hands already to avoid it being recomputed here.

        init : {'k-means++', 'random'}, callable or ndarray of shape \
                (n_clusters, n_features)
            Method for initialization.

        random_state : RandomState instance
            Determines random number generation for centroid initialization.
            See :term:`Glossary <random_state>`.

        init_size : int, default=None
            Number of samples to randomly sample for speeding up the
            initialization (sometimes at the expense of accuracy).

        Returns
        -------
        centers : ndarray of shape (n_clusters, n_features)
        """
        n_samples = X.shape[0]
        n_clusters = self.n_clusters

        if init_size is not None and init_size < n_samples:
            init_indices = random_state.randint(0, n_samples, init_size)
            X = X[init_indices]
            x_squared_norms = x_squared_norms[init_indices]
            n_samples = X.shape[0]

        if isinstance(init, str) and init == 'k-means++':
            centers, _ = _kmeans_plusplus(X, n_clusters,
                                          random_state=random_state,
                                          x_squared_norms=x_squared_norms)
        elif isinstance(init, str) and init == 'random':
            seeds = random_state.permutation(n_samples)[:n_clusters]
            centers = X[seeds]
        elif hasattr(init, '__array__'):
            centers = init
        elif callable(init):
            centers = init(X, n_clusters, random_state=random_state)
            centers = check_array(
                centers, dtype=X.dtype, copy=False, order='C')
            self._validate_center_shape(X, centers)

        if sp.issparse(centers):
            centers = centers.toarray()

        return centers

    def fit(self, X, y=None, sample_weight=None):
        """Compute k-means clustering.

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

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

            .. versionadded:: 0.20

        Returns
        -------
        self
            Fitted estimator.
        """
        X = self._validate_data(X, accept_sparse='csr',
                                dtype=[np.float64, np.float32],
                                order='C', copy=self.copy_x,
                                accept_large_sparse=False)

        self._check_params(X)
        random_state = check_random_state(self.random_state)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        # Validate init array
        init = self.init
        if hasattr(init, '__array__'):
            init = check_array(init, dtype=X.dtype, copy=True, order='C')
            self._validate_center_shape(X, init)

        # subtract of mean of x for more accurate distance computations
        if not sp.issparse(X):
            X_mean = X.mean(axis=0)
            # The copy was already done above
            X -= X_mean

            if hasattr(init, '__array__'):
                init -= X_mean

        # precompute squared norms of data points
        x_squared_norms = row_norms(X, squared=True)

        if self._algorithm == "full":
            kmeans_single = _kmeans_single_lloyd
            self._check_mkl_vcomp(X, X.shape[0])
        else:
            kmeans_single = _kmeans_single_elkan

        best_inertia = None

        for i in range(self._n_init):
            # Initialize centers
            centers_init = self._init_centroids(
                X, x_squared_norms=x_squared_norms, init=init,
                random_state=random_state)
            if self.verbose:
                print("Initialization complete")

            # run a k-means once
            labels, inertia, centers, n_iter_ = kmeans_single(
                X, sample_weight, centers_init, max_iter=self.max_iter,
                verbose=self.verbose, tol=self._tol,
                x_squared_norms=x_squared_norms, n_threads=self._n_threads)

            # determine if these results are the best so far
            if best_inertia is None or inertia < best_inertia:
                best_labels = labels
                best_centers = centers
                best_inertia = inertia
                best_n_iter = n_iter_

        if not sp.issparse(X):
            if not self.copy_x:
                X += X_mean
            best_centers += X_mean

        distinct_clusters = len(set(best_labels))
        if distinct_clusters < self.n_clusters:
            warnings.warn(
                "Number of distinct clusters ({}) found smaller than "
                "n_clusters ({}). Possibly due to duplicate points "
                "in X.".format(distinct_clusters, self.n_clusters),
                ConvergenceWarning, stacklevel=2)

        self.cluster_centers_ = best_centers
        self.labels_ = best_labels
        self.inertia_ = best_inertia
        self.n_iter_ = best_n_iter
        return self

    def fit_predict(self, X, y=None, sample_weight=None):
        """Compute cluster centers and predict cluster index for each sample.

        Convenience method; equivalent to calling fit(X) followed by
        predict(X).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to transform.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        return self.fit(X, sample_weight=sample_weight).labels_

    def fit_transform(self, X, y=None, sample_weight=None):
        """Compute clustering and transform X to cluster-distance space.

        Equivalent to fit(X).transform(X), but more efficiently implemented.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to transform.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_clusters)
            X transformed in the new space.
        """
        # Currently, this just skips a copy of the data if it is not in
        # np.array or CSR format already.
        # XXX This skips _check_test_data, which may change the dtype;
        # we should refactor the input validation.
        return self.fit(X, sample_weight=sample_weight)._transform(X)

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
        return euclidean_distances(X, self.cluster_centers_)

    def predict(self, X, sample_weight=None):
        """Predict the closest cluster each sample in X belongs to.

        In the vector quantization literature, `cluster_centers_` is called
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

        X = self._check_test_data(X)
        x_squared_norms = row_norms(X, squared=True)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        return _labels_inertia(X, sample_weight, x_squared_norms,
                               self.cluster_centers_, self._n_threads)[0]

    def score(self, X, y=None, sample_weight=None):
        """Opposite of the value of X on the K-means objective.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight.

        Returns
        -------
        score : float
            Opposite of the value of X on the K-means objective.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        x_squared_norms = row_norms(X, squared=True)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        return -_labels_inertia(X, sample_weight, x_squared_norms,
                                self.cluster_centers_)[1]

    def _more_tags(self):
        return {
            '_xfail_checks': {
                'check_sample_weights_invariance':
                'zero sample_weight is not equivalent to removing samples',
            },
        }


def _mini_batch_step(X, sample_weight, x_squared_norms, centers, weight_sums,
                     old_center_buffer, compute_squared_diff,
                     distances, random_reassign=False,
                     random_state=None, reassignment_ratio=.01,
                     verbose=False):
    """Incremental update of the centers for the Minibatch K-Means algorithm.

    Parameters
    ----------

    X : ndarray of shape (n_samples, n_features)
        The original data array.

    sample_weight : array-like of shape (n_samples,)
        The weights for each observation in X.

    x_squared_norms : ndarray of shape (n_samples,)
        Squared euclidean norm of each data point.

    centers : ndarray of shape (k, n_features)
        The cluster centers. This array is MODIFIED IN PLACE

    old_center_buffer : int
        Copy of old centers for monitoring convergence.

    compute_squared_diff : bool
        If set to False, the squared diff computation is skipped.

    distances : ndarray of shape (n_samples,), dtype=float, default=None
        If not None, should be a pre-allocated array that will be used to store
        the distances of each sample to its closest center.
        May not be None when random_reassign is True.

    random_reassign : bool, default=False
        If True, centers with very low counts are randomly reassigned
        to observations.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization and to
        pick new clusters amongst observations with uniform probability. Use
        an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    reassignment_ratio : float, default=.01
        Control the fraction of the maximum number of counts for a
        center to be reassigned. A higher value means that low count
        centers are more likely to be reassigned, which means that the
        model will take longer to converge, but should converge in a
        better clustering.

    verbose : bool, default=False
        Controls the verbosity.

    Returns
    -------
    inertia : float
        Sum of squared distances of samples to their closest cluster center.

    squared_diff : ndarray of shape (n_clusters,)
        Squared distances between previous and updated cluster centers.

    """
    # Perform label assignment to nearest centers
    nearest_center, inertia = _labels_inertia(X, sample_weight,
                                              x_squared_norms, centers)

    if random_reassign and reassignment_ratio > 0:
        random_state = check_random_state(random_state)
        # Reassign clusters that have very low weight
        to_reassign = weight_sums < reassignment_ratio * weight_sums.max()
        # pick at most .5 * batch_size samples as new centers
        if to_reassign.sum() > .5 * X.shape[0]:
            indices_dont_reassign = \
                    np.argsort(weight_sums)[int(.5 * X.shape[0]):]
            to_reassign[indices_dont_reassign] = False
        n_reassigns = to_reassign.sum()
        if n_reassigns:
            # Pick new clusters amongst observations with uniform probability
            new_centers = random_state.choice(X.shape[0], replace=False,
                                              size=n_reassigns)
            if verbose:
                print("[MiniBatchKMeans] Reassigning %i cluster centers."
                      % n_reassigns)

            if sp.issparse(X) and not sp.issparse(centers):
                assign_rows_csr(
                        X, new_centers.astype(np.intp, copy=False),
                        np.where(to_reassign)[0].astype(np.intp, copy=False),
                        centers)
            else:
                centers[to_reassign] = X[new_centers]
        # reset counts of reassigned centers, but don't reset them too small
        # to avoid instant reassignment. This is a pretty dirty hack as it
        # also modifies the learning rates.
        weight_sums[to_reassign] = np.min(weight_sums[~to_reassign])

    # implementation for the sparse CSR representation completely written in
    # cython
    if sp.issparse(X):
        return inertia, _mini_batch_update_csr(
            X, sample_weight, x_squared_norms, centers, weight_sums,
            nearest_center, old_center_buffer, compute_squared_diff)

    # dense variant in mostly numpy (not as memory efficient though)
    k = centers.shape[0]
    squared_diff = 0.0
    for center_idx in range(k):
        # find points from minibatch that are assigned to this center
        center_mask = nearest_center == center_idx
        wsum = sample_weight[center_mask].sum()

        if wsum > 0:
            if compute_squared_diff:
                old_center_buffer[:] = centers[center_idx]

            # inplace remove previous count scaling
            centers[center_idx] *= weight_sums[center_idx]

            # inplace sum with new points members of this cluster
            centers[center_idx] += \
                np.sum(X[center_mask] *
                       sample_weight[center_mask, np.newaxis], axis=0)

            # update the count statistics for this center
            weight_sums[center_idx] += wsum

            # inplace rescale to compute mean of all points (old and new)
            # Note: numpy >= 1.10 does not support '/=' for the following
            # expression for a mixture of int and float (see numpy issue #6464)
            centers[center_idx] = centers[center_idx] / weight_sums[center_idx]

            # update the squared diff if necessary
            if compute_squared_diff:
                diff = centers[center_idx].ravel() - old_center_buffer.ravel()
                squared_diff += np.dot(diff, diff)

    return inertia, squared_diff


def _mini_batch_convergence(model, iteration_idx, n_iter, tol,
                            n_samples, centers_squared_diff, batch_inertia,
                            context, verbose=0):
    """Helper function to encapsulate the early stopping logic."""
    # Normalize inertia to be able to compare values when
    # batch_size changes
    batch_inertia /= model.batch_size
    centers_squared_diff /= model.batch_size

    # Compute an Exponentially Weighted Average of the squared
    # diff to monitor the convergence while discarding
    # minibatch-local stochastic variability:
    # https://en.wikipedia.org/wiki/Moving_average
    ewa_diff = context.get('ewa_diff')
    ewa_inertia = context.get('ewa_inertia')
    if ewa_diff is None:
        ewa_diff = centers_squared_diff
        ewa_inertia = batch_inertia
    else:
        alpha = float(model.batch_size) * 2.0 / (n_samples + 1)
        alpha = 1.0 if alpha > 1.0 else alpha
        ewa_diff = ewa_diff * (1 - alpha) + centers_squared_diff * alpha
        ewa_inertia = ewa_inertia * (1 - alpha) + batch_inertia * alpha

    # Log progress to be able to monitor convergence
    if verbose:
        progress_msg = (
            'Minibatch iteration %d/%d:'
            ' mean batch inertia: %f, ewa inertia: %f ' % (
                iteration_idx + 1, n_iter, batch_inertia,
                ewa_inertia))
        print(progress_msg)

    # Early stopping based on absolute tolerance on squared change of
    # centers position (using EWA smoothing)
    if tol > 0.0 and ewa_diff <= tol:
        if verbose:
            print('Converged (small centers change) at iteration %d/%d'
                  % (iteration_idx + 1, n_iter))
        return True

    # Early stopping heuristic due to lack of improvement on smoothed inertia
    ewa_inertia_min = context.get('ewa_inertia_min')
    no_improvement = context.get('no_improvement', 0)
    if ewa_inertia_min is None or ewa_inertia < ewa_inertia_min:
        no_improvement = 0
        ewa_inertia_min = ewa_inertia
    else:
        no_improvement += 1

    if (model.max_no_improvement is not None
            and no_improvement >= model.max_no_improvement):
        if verbose:
            print('Converged (lack of improvement in inertia)'
                  ' at iteration %d/%d'
                  % (iteration_idx + 1, n_iter))
        return True

    # update the convergence context to maintain state across successive calls:
    context['ewa_diff'] = ewa_diff
    context['ewa_inertia'] = ewa_inertia
    context['ewa_inertia_min'] = ewa_inertia_min
    context['no_improvement'] = no_improvement
    return False


class MiniBatchKMeans(KMeans):
    """
    Mini-Batch K-Means clustering.

    Read more in the :ref:`User Guide <mini_batch_kmeans>`.

    Parameters
    ----------

    n_clusters : int, default=8
        The number of clusters to form as well as the number of
        centroids to generate.

    init : {'k-means++', 'random'}, callable or array-like of shape \
            (n_clusters, n_features), default='k-means++'
        Method for initialization:

        'k-means++' : selects initial cluster centers for k-mean
        clustering in a smart way to speed up convergence. See section
        Notes in k_init for more details.

        'random': choose `n_clusters` observations (rows) at random from data
        for the initial centroids.

        If an array is passed, it should be of shape (n_clusters, n_features)
        and gives the initial centers.

        If a callable is passed, it should take arguments X, n_clusters and a
        random state and return an initialization.

    max_iter : int, default=100
        Maximum number of iterations over the complete dataset before
        stopping independently of any early stopping criterion heuristics.

    batch_size : int, default=100
        Size of the mini batches.

    verbose : int, default=0
        Verbosity mode.

    compute_labels : bool, default=True
        Compute label assignment and inertia for the complete dataset
        once the minibatch optimization has converged in fit.

    random_state : int, RandomState instance or None, default=None
        Determines random number generation for centroid initialization and
        random reassignment. Use an int to make the randomness deterministic.
        See :term:`Glossary <random_state>`.

    tol : float, default=0.0
        Control early stopping based on the relative center changes as
        measured by a smoothed, variance-normalized of the mean center
        squared position changes. This early stopping heuristics is
        closer to the one used for the batch variant of the algorithms
        but induces a slight computational and memory overhead over the
        inertia heuristic.

        To disable convergence detection based on normalized center
        change, set tol to 0.0 (default).

    max_no_improvement : int, default=10
        Control early stopping based on the consecutive number of mini
        batches that does not yield an improvement on the smoothed inertia.

        To disable convergence detection based on inertia, set
        max_no_improvement to None.

    init_size : int, default=None
        Number of samples to randomly sample for speeding up the
        initialization (sometimes at the expense of accuracy): the
        only algorithm is initialized by running a batch KMeans on a
        random subset of the data. This needs to be larger than n_clusters.

        If `None`, `init_size= 3 * batch_size`.

    n_init : int, default=3
        Number of random initializations that are tried.
        In contrast to KMeans, the algorithm is only run once, using the
        best of the ``n_init`` initializations as measured by inertia.

    reassignment_ratio : float, default=0.01
        Control the fraction of the maximum number of counts for a
        center to be reassigned. A higher value means that low count
        centers are more easily reassigned, which means that the
        model will take longer to converge, but should converge in a
        better clustering.

    Attributes
    ----------

    cluster_centers_ : ndarray of shape (n_clusters, n_features)
        Coordinates of cluster centers.

    labels_ : int
        Labels of each point (if compute_labels is set to True).

    inertia_ : float
        The value of the inertia criterion associated with the chosen
        partition (if compute_labels is set to True). The inertia is
        defined as the sum of square distances of samples to their nearest
        neighbor.

    n_iter_ : int
        Number of batches processed.

    counts_ : ndarray of shape (n_clusters,)
        Weigth sum of each cluster.

        .. deprecated:: 0.24
           This attribute is deprecated in 0.24 and will be removed in
           1.1 (renaming of 0.26).

    init_size_ : int
        The effective number of samples used for the initialization.

        .. deprecated:: 0.24
           This attribute is deprecated in 0.24 and will be removed in
           1.1 (renaming of 0.26).

    See Also
    --------
    KMeans : The classic implementation of the clustering method based on the
        Lloyd's algorithm. It consumes the whole set of input data at each
        iteration.

    Notes
    -----
    See https://www.eecs.tufts.edu/~dsculley/papers/fastkmeans.pdf

    Examples
    --------
    >>> from sklearn.cluster import MiniBatchKMeans
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [4, 2], [4, 0], [4, 4],
    ...               [4, 5], [0, 1], [2, 2],
    ...               [3, 2], [5, 5], [1, -1]])
    >>> # manually fit on batches
    >>> kmeans = MiniBatchKMeans(n_clusters=2,
    ...                          random_state=0,
    ...                          batch_size=6)
    >>> kmeans = kmeans.partial_fit(X[0:6,:])
    >>> kmeans = kmeans.partial_fit(X[6:12,:])
    >>> kmeans.cluster_centers_
    array([[2. , 1. ],
           [3.5, 4.5]])
    >>> kmeans.predict([[0, 0], [4, 4]])
    array([0, 1], dtype=int32)
    >>> # fit on the whole data
    >>> kmeans = MiniBatchKMeans(n_clusters=2,
    ...                          random_state=0,
    ...                          batch_size=6,
    ...                          max_iter=10).fit(X)
    >>> kmeans.cluster_centers_
    array([[3.95918367, 2.40816327],
           [1.12195122, 1.3902439 ]])
    >>> kmeans.predict([[0, 0], [4, 4]])
    array([1, 0], dtype=int32)
    """
    @_deprecate_positional_args
    def __init__(self, n_clusters=8, *, init='k-means++', max_iter=100,
                 batch_size=100, verbose=0, compute_labels=True,
                 random_state=None, tol=0.0, max_no_improvement=10,
                 init_size=None, n_init=3, reassignment_ratio=0.01):

        super().__init__(
            n_clusters=n_clusters, init=init, max_iter=max_iter,
            verbose=verbose, random_state=random_state, tol=tol, n_init=n_init)

        self.max_no_improvement = max_no_improvement
        self.batch_size = batch_size
        self.compute_labels = compute_labels
        self.init_size = init_size
        self.reassignment_ratio = reassignment_ratio

    @deprecated("The attribute 'counts_' is deprecated in 0.24"  # type: ignore
                " and will be removed in 1.1 (renaming of 0.26).")
    @property
    def counts_(self):
        return self._counts

    @deprecated("The attribute 'init_size_' is deprecated in "  # type: ignore
                "0.24 and will be removed in 1.1 (renaming of 0.26).")
    @property
    def init_size_(self):
        return self._init_size

    @deprecated("The attribute 'random_state_' is deprecated "  # type: ignore
                "in 0.24 and will be removed in 1.1 (renaming of 0.26).")
    @property
    def random_state_(self):
        return getattr(self, "_random_state", None)

    def _check_params(self, X):
        super()._check_params(X)

        # max_no_improvement
        if self.max_no_improvement is not None and self.max_no_improvement < 0:
            raise ValueError(
                f"max_no_improvement should be >= 0, got "
                f"{self.max_no_improvement} instead.")

        # batch_size
        if self.batch_size <= 0:
            raise ValueError(
                f"batch_size should be > 0, got {self.batch_size} instead.")

        # init_size
        if self.init_size is not None and self.init_size <= 0:
            raise ValueError(
                f"init_size should be > 0, got {self.init_size} instead.")
        self._init_size = self.init_size
        if self._init_size is None:
            self._init_size = 3 * self.batch_size
            if self._init_size < self.n_clusters:
                self._init_size = 3 * self.n_clusters
        elif self._init_size < self.n_clusters:
            warnings.warn(
                f"init_size={self._init_size} should be larger than "
                f"n_clusters={self.n_clusters}. Setting it to "
                f"min(3*n_clusters, n_samples)",
                RuntimeWarning, stacklevel=2)
            self._init_size = 3 * self.n_clusters
        self._init_size = min(self._init_size, X.shape[0])

        # reassignment_ratio
        if self.reassignment_ratio < 0:
            raise ValueError(
                f"reassignment_ratio should be >= 0, got "
                f"{self.reassignment_ratio} instead.")

    def fit(self, X, y=None, sample_weight=None):
        """Compute the centroids on X by chunking it into mini-batches.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training instances to cluster. It must be noted that the data
            will be converted to C ordering, which will cause a memory copy
            if the given data is not C-contiguous.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

            .. versionadded:: 0.20

        Returns
        -------
        self
        """
        X = self._validate_data(X, accept_sparse='csr',
                                dtype=[np.float64, np.float32],
                                order='C', accept_large_sparse=False)

        self._check_params(X)
        random_state = check_random_state(self.random_state)
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        # Validate init array
        init = self.init
        if hasattr(init, '__array__'):
            init = check_array(init, dtype=X.dtype, copy=True, order='C')
            self._validate_center_shape(X, init)

        n_samples, n_features = X.shape
        x_squared_norms = row_norms(X, squared=True)

        if self.tol > 0.0:
            tol = _tolerance(X, self.tol)

            # using tol-based early stopping needs the allocation of a
            # dedicated before which can be expensive for high dim data:
            # hence we allocate it outside of the main loop
            old_center_buffer = np.zeros(n_features, dtype=X.dtype)
        else:
            tol = 0.0
            # no need for the center buffer if tol-based early stopping is
            # disabled
            old_center_buffer = np.zeros(0, dtype=X.dtype)

        distances = np.zeros(self.batch_size, dtype=X.dtype)
        n_batches = int(np.ceil(float(n_samples) / self.batch_size))
        n_iter = int(self.max_iter * n_batches)

        self._check_mkl_vcomp(X, self.batch_size)

        validation_indices = random_state.randint(0, n_samples,
                                                  self._init_size)
        X_valid = X[validation_indices]
        sample_weight_valid = sample_weight[validation_indices]
        x_squared_norms_valid = x_squared_norms[validation_indices]

        # perform several inits with random sub-sets
        best_inertia = None
        for init_idx in range(self._n_init):
            if self.verbose:
                print("Init %d/%d with method: %s"
                      % (init_idx + 1, self._n_init, init))
            weight_sums = np.zeros(self.n_clusters, dtype=sample_weight.dtype)

            # TODO: once the `k_means` function works with sparse input we
            # should refactor the following init to use it instead.

            # Initialize the centers using only a fraction of the data as we
            # expect n_samples to be very large when using MiniBatchKMeans
            cluster_centers = self._init_centroids(
                X, x_squared_norms=x_squared_norms,
                init=init,
                random_state=random_state,
                init_size=self._init_size)

            # Compute the label assignment on the init dataset
            _mini_batch_step(
                X_valid, sample_weight_valid,
                x_squared_norms[validation_indices], cluster_centers,
                weight_sums, old_center_buffer, False, distances=None,
                verbose=self.verbose)

            # Keep only the best cluster centers across independent inits on
            # the common validation set
            _, inertia = _labels_inertia(X_valid, sample_weight_valid,
                                         x_squared_norms_valid,
                                         cluster_centers)
            if self.verbose:
                print("Inertia for init %d/%d: %f"
                      % (init_idx + 1, self._n_init, inertia))
            if best_inertia is None or inertia < best_inertia:
                self.cluster_centers_ = cluster_centers
                self._counts = weight_sums
                best_inertia = inertia

        # Empty context to be used inplace by the convergence check routine
        convergence_context = {}

        # Perform the iterative optimization until the final convergence
        # criterion
        for iteration_idx in range(n_iter):
            # Sample a minibatch from the full dataset
            minibatch_indices = random_state.randint(
                0, n_samples, self.batch_size)

            # Perform the actual update step on the minibatch data
            batch_inertia, centers_squared_diff = _mini_batch_step(
                X[minibatch_indices], sample_weight[minibatch_indices],
                x_squared_norms[minibatch_indices],
                self.cluster_centers_, self._counts,
                old_center_buffer, tol > 0.0, distances=distances,
                # Here we randomly choose whether to perform
                # random reassignment: the choice is done as a function
                # of the iteration index, and the minimum number of
                # counts, in order to force this reassignment to happen
                # every once in a while
                random_reassign=((iteration_idx + 1)
                                 % (10 + int(self._counts.min())) == 0),
                random_state=random_state,
                reassignment_ratio=self.reassignment_ratio,
                verbose=self.verbose)

            # Monitor convergence and do early stopping if necessary
            if _mini_batch_convergence(
                    self, iteration_idx, n_iter, tol, n_samples,
                    centers_squared_diff, batch_inertia, convergence_context,
                    verbose=self.verbose):
                break

        self.n_iter_ = iteration_idx + 1

        if self.compute_labels:
            self.labels_, self.inertia_ = \
                    self._labels_inertia_minibatch(X, sample_weight)

        return self

    def _labels_inertia_minibatch(self, X, sample_weight):
        """Compute labels and inertia using mini batches.

        This is slightly slower than doing everything at once but prevents
        memory errors / segfaults.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        sample_weight : array-like of shape (n_samples,)
            The weights for each observation in X.

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Cluster labels for each point.

        inertia : float
            Sum of squared distances of points to nearest cluster.
        """
        if self.verbose:
            print('Computing label assignment and total inertia')
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)
        x_squared_norms = row_norms(X, squared=True)
        slices = gen_batches(X.shape[0], self.batch_size)
        results = [_labels_inertia(X[s], sample_weight[s], x_squared_norms[s],
                                   self.cluster_centers_) for s in slices]
        labels, inertia = zip(*results)
        return np.hstack(labels), np.sum(inertia)

    def partial_fit(self, X, y=None, sample_weight=None):
        """Update k means estimate on a single mini-batch X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Coordinates of the data points to cluster. It must be noted that
            X will be copied if it is not C-contiguous.

        y : Ignored
            Not used, present here for API consistency by convention.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        self
        """
        is_first_call_to_partial_fit = not hasattr(self, 'cluster_centers_')

        X = self._validate_data(X, accept_sparse='csr',
                                dtype=[np.float64, np.float32],
                                order='C', accept_large_sparse=False,
                                reset=is_first_call_to_partial_fit)

        self._random_state = getattr(self, "_random_state",
                                     check_random_state(self.random_state))
        sample_weight = _check_sample_weight(sample_weight, X, dtype=X.dtype)

        x_squared_norms = row_norms(X, squared=True)

        if is_first_call_to_partial_fit:
            # this is the first call to partial_fit on this object
            self._check_params(X)

            # Validate init array
            init = self.init
            if hasattr(init, '__array__'):
                init = check_array(init, dtype=X.dtype, copy=True, order='C')
                self._validate_center_shape(X, init)

            self._check_mkl_vcomp(X, X.shape[0])

            # initialize the cluster centers
            self.cluster_centers_ = self._init_centroids(
                X, x_squared_norms=x_squared_norms,
                init=init,
                random_state=self._random_state,
                init_size=self._init_size)

            self._counts = np.zeros(self.n_clusters,
                                    dtype=sample_weight.dtype)
            random_reassign = False
            distances = None
        else:
            # The lower the minimum count is, the more we do random
            # reassignment, however, we don't want to do random
            # reassignment too often, to allow for building up counts
            random_reassign = self._random_state.randint(
                10 * (1 + self._counts.min())) == 0
            distances = np.zeros(X.shape[0], dtype=X.dtype)

        _mini_batch_step(X, sample_weight, x_squared_norms,
                         self.cluster_centers_, self._counts,
                         np.zeros(0, dtype=X.dtype), 0,
                         random_reassign=random_reassign, distances=distances,
                         random_state=self._random_state,
                         reassignment_ratio=self.reassignment_ratio,
                         verbose=self.verbose)

        if self.compute_labels:
            self.labels_, self.inertia_ = _labels_inertia(
                X, sample_weight, x_squared_norms, self.cluster_centers_)

        return self

    def predict(self, X, sample_weight=None):
        """Predict the closest cluster each sample in X belongs to.

        In the vector quantization literature, `cluster_centers_` is called
        the code book and each value returned by `predict` is the index of
        the closest code in the code book.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            New data to predict.

        sample_weight : array-like of shape (n_samples,), default=None
            The weights for each observation in X. If None, all observations
            are assigned equal weight (default: None).

        Returns
        -------
        labels : ndarray of shape (n_samples,)
            Index of the cluster each sample belongs to.
        """
        check_is_fitted(self)

        X = self._check_test_data(X)
        return self._labels_inertia_minibatch(X, sample_weight)[0]

    def _more_tags(self):
        return {
            '_xfail_checks': {
                'check_sample_weights_invariance':
                'zero sample_weight is not equivalent to removing samples',
            }
        }


def kmeans_plusplus(X, n_clusters, *, x_squared_norms=None,
                    random_state=None, n_local_trials=None):
    """Init n_clusters seeds according to k-means++

    .. versionadded:: 0.24

    Parameters
    ----------
    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        The data to pick seeds from.

    n_clusters : int
        The number of centroids to initialize

    x_squared_norms : array-like of shape (n_samples,), default=None
        Squared Euclidean norm of each data point.

    random_state : int or RandomState instance, default=None
        Determines random number generation for centroid initialization. Pass
        an int for reproducible output across multiple function calls.
        See :term:`Glossary <random_state>`.

    n_local_trials : int, default=None
        The number of seeding trials for each center (except the first),
        of which the one reducing inertia the most is greedily chosen.
        Set to None to make the number of trials depend logarithmically
        on the number of seeds (2+log(k)).

    Returns
    -------
    centers : ndarray of shape (n_clusters, n_features)
        The inital centers for k-means.

    indices : ndarray of shape (n_clusters,)
        The index location of the chosen centers in the data array X. For a
        given index and center, X[index] = center.

    Notes
    -----
    Selects initial cluster centers for k-mean clustering in a smart way
    to speed up convergence. see: Arthur, D. and Vassilvitskii, S.
    "k-means++: the advantages of careful seeding". ACM-SIAM symposium
    on Discrete algorithms. 2007

    Examples
    --------

    >>> from sklearn.cluster import kmeans_plusplus
    >>> import numpy as np
    >>> X = np.array([[1, 2], [1, 4], [1, 0],
    ...               [10, 2], [10, 4], [10, 0]])
    >>> centers, indices = kmeans_plusplus(X, n_clusters=2, random_state=0)
    >>> centers
    array([[10,  4],
           [ 1,  0]])
    >>> indices
    array([4, 2])
    """

    # Check data
    check_array(X, accept_sparse='csr',
                dtype=[np.float64, np.float32])

    if X.shape[0] < n_clusters:
        raise ValueError(f"n_samples={X.shape[0]} should be >= "
                         f"n_clusters={n_clusters}.")

    # Check parameters
    if x_squared_norms is None:
        x_squared_norms = row_norms(X, squared=True)
    else:
        x_squared_norms = check_array(x_squared_norms,
                                      dtype=X.dtype,
                                      ensure_2d=False)

    if x_squared_norms.shape[0] != X.shape[0]:
        raise ValueError(
            f"The length of x_squared_norms {x_squared_norms.shape[0]} should "
            f"be equal to the length of n_samples {X.shape[0]}.")

    if n_local_trials is not None and n_local_trials < 1:
        raise ValueError(
            f"n_local_trials is set to {n_local_trials} but should be an "
            f"integer value greater than zero.")

    random_state = check_random_state(random_state)

    # Call private k-means++
    centers, indices = _kmeans_plusplus(X, n_clusters, x_squared_norms,
                                        random_state, n_local_trials)

    return centers, indices
