import numpy as np
import scipy.spatial.distance
import warnings
from joblib import Parallel, delayed

from collections import defaultdict
from ..utils.validation import check_is_fitted
from ..utils import check_random_state, gen_batches, check_array
from ..base import BaseEstimator, ClusterMixin
from ..neighbors import NearestNeighbors
from ..metrics.pairwise import pairwise_distances_argmin
from .mean_shift_ import estimate_bandwidth
from ..gaussian_process.kernels import RBF

# These are multiples of the bandwidth
CONVERGENCE_THRESHOLD = 1e-7
MERGE_THRESHOLD = 1.
# burn in b/c there is no guarantee that the first iterations will be the largest
BURN_IN = 100

# separate function for each seed's iterative loop
def _mean_shift_single(x, train_X, nbrs, kernel, bandwidth, max_iter):
    # For each x, climb gradient until convergence or max_iter
    radius = nbrs.get_params()['radius']
    stop_thresh = CONVERGENCE_THRESHOLD * bandwidth # when mean has converged
    completed_iterations = 0
    while True:
        # Find mean of points within bandwidth
        d_nbrs, i_nbrs = nbrs.radius_neighbors([x], radius,
                                       return_distance=True)
        d_nbrs = d_nbrs[0]
        i_nbrs = i_nbrs[0]

        points_within = train_X[i_nbrs]
        weights = kernel(d_nbrs.reshape(-1,1), [[0]]).flatten()
    
        if len(points_within) == 0:
            break  # Depending on kernel type and bandwidth this may occur
        old_x = x  # save the old point
        x = np.average(points_within, weights=weights, axis=0)
        distance_traveled = np.squeeze(scipy.spatial.distance.cdist(
            [x],[old_x], metric=nbrs.metric, **(nbrs.metric_params or {})))
        if (    (completed_iterations >= BURN_IN)
            and (  (distance_traveled < stop_thresh)
                or (completed_iterations >= max_iter))):
            break
        completed_iterations += 1
    return x

class MeanShiftv2(ClusterMixin, BaseEstimator):
    """Mean shift clustering

    Mean shift clustering aims to discover "blobs" in a smooth density of
    samples. It is a mode-finding algorithm on the kernel density function:
    sklearn.neighbors.KernelDensity 

    Read more in the :ref:`User Guide <mean_shift>`.

    Parameters
    ----------
    bandwidth : float, optional
        Bandwidth used in the RBF kernel.

        If not given, the bandwidth is estimated using
        sklearn.cluster.estimate_bandwidth; see the documentation for that
        function for hints on scalability (see also the Notes, below).

    nbrs : Object implementing the same interface as sklearn.neighbors.NearestNeighbors
        This is used for fast nearest-neighbor lookup

        If not given, a default parameter sklearn.neighbors.NearestNeighbors is used.

    kernel : Object implementing sklearn.gaussian_process.kernels.Kernel
        This is used to turn distances into weights when computing the
        locally-weighted mean.  

        If not given, an sklearn.gaussian_process.kernels.RBF is used.

    n_jobs : int or None, optional (default=None)
        The number of jobs to use for the computation. This works by computing
        each of the n_init runs in parallel.

        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    max_iter : int, default=800
        Maximum number of iterations, per seed point before the clustering
        operation terminates (for that seed point), if has not converged yet.

        .. versionadded:: 0.22

    Attributes
    ----------
    cluster_centers_ : array, [n_clusters, n_features]
        Coordinates of cluster centers.

    Examples
    --------
    >>> from sklearn.cluster import MeanShiftv2
    >>> import numpy as np
    >>> X = np.array([[1, 1], [2, 1], [1, 0],
    ...               [4, 7], [3, 5], [3, 6]])
    >>> clustering = MeanShiftv2(bandwidth=2).fit(X)
    >>> clustering.predict([[0, 0], [5, 5]])
    array([1, 0])
    >>> clustering
    MeanShiftv2(bandwidth=2)

    Notes
    -----

    Scalability:

    When using a Ball Tree to look up members of each kernel (the default), 
    the complexity will tend towards O(T*n*log(n)) in lower dimensions, 
    with n the number of samples and T the number of points. In higher 
    dimensions the complexity will tend towards O(T*n^2).

    Scalability can be boosted by using a small but representative training
    set, and by using a kernel with finite support.

    Note that the estimate_bandwidth function is much less scalable than the
    mean shift algorithm and will be the bottleneck if it is used.

    References
    ----------

    Dorin Comaniciu and Peter Meer, "Mean Shift: A robust approach toward
    feature space analysis". IEEE Transactions on Pattern Analysis and
    Machine Intelligence. 2002. pp. 603-619.

    """
    def __init__(self, bandwidth=None, nbrs=None, kernel=None, n_jobs=None, max_iter=800):
        self.bandwidth = bandwidth
        self.nbrs = nbrs
        self.kernel = kernel
        self.n_jobs = n_jobs
        self.max_iter = max_iter

    def fit(self, X, y=None):
        """
            Store X.
            Also, set all the defaults if they haven't been set already

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Samples to cluster.

        y : Ignored

        """
        X = check_array(X)
        self.train_X = X
        if self.bandwidth is None:
            self.bandwidth = estimate_bandwidth(X, n_jobs=self.n_jobs)
        elif self.bandwidth <= 0:
            raise ValueError("bandwidth needs to be greater than zero or None,"
                             " got %f" % self.bandwidth)

        if self.nbrs is None:
            #NB: 3 standard deviations is a "most of" an RBF kernel
            #  1 std deviation is a very truncated normal, but would be faster
            #  To get a finite-support kernel, something like a tricube kernel would be good
            self.nbrs = NearestNeighbors(radius=3*self.bandwidth, n_jobs=1).fit(X)

        if self.kernel is None:
            self.kernel = RBF(length_scale=self.bandwidth)

        self.cluster_centers_ = []

    def predict(self, X):

        check_is_fitted(self)

        # execute iterations on all points in parallel
        #TODO this could be sped up a bit by doing stuff in batches, I think
        all_res = Parallel(n_jobs=self.n_jobs)(
            delayed(_mean_shift_single)
            (x, self.train_X, self.nbrs, self.kernel, self.bandwidth, self.max_iter) for x in X)

        # POST PROCESSING: remove near duplicate points
        # If the distance to an already-known center is "small", assign that cluster
        # Else, we have found a new cluster, so label it

        if len(self.cluster_centers_) == 0:
            self.cluster_centers_.append(all_res[0])

        labels = []
        thresh = MERGE_THRESHOLD * self.bandwidth
        for res in all_res:
            pairwise_distances = scipy.spatial.distance.cdist(
                self.cluster_centers_, [res], metric=self.nbrs.metric, **(self.nbrs.metric_params or {}))
            pairwise_distances = pairwise_distances.flatten()
            closest_center = np.argmin(pairwise_distances)
            if pairwise_distances[closest_center] > thresh: # a new cluster
                self.cluster_centers_.append(res)
                labels.append(len(self.cluster_centers_))
            else: # a previously known cluster
                labels.append(closest_center)
        
        return np.array(labels)
