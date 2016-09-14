# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import numpy as np
from scipy.stats import scoreatpercentile

from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import UnsupervisedMixin

from ..utils.validation import check_is_fitted
from ..utils import check_array

__all__ = ["LocalOutlierFactor"]


class LocalOutlierFactor(NeighborsBase, KNeighborsMixin, UnsupervisedMixin):
    """Unsupervised Outlier Detection using Local Outlier Factor (LOF)

    The anomaly score of each sample is called Local Outlier Factor.
    It measures the local deviation of density of a given sample with
    respect to its neighbors.
    It is local in that the anomaly score depends on how isolated the object
    is with respect to the surrounding neighborhood.
    More precisely, locality is given by k-nearest neighbors, whose distance
    is used to estimate the local density.
    By comparing the local density of a sample to the local densities of
    its neighbors, one can identify samples that have a substantially lower
    density than their neighbors. These are considered as outliers.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`k_neighbors` queries.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    p: integer, optional (default = 2)
        Parameter for the Minkowski metric from
        :ref:`sklearn.metrics.pairwise.pairwise_distances`. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric : string or callable, default 'minkowski'
        metric used for the distance computation. Any metric from scikit-learn
        or scipy.spatial.distance can be used.

        If 'precomputed', the training input X is expected to be a distance
        matrix.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Distance matrices are not supported.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
          'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
          'sqeuclidean', 'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics.

    metric_params : dict, optional (default = None)
        Additional keyword arguments for the metric function.

    contamination : float in (0., 0.5), optional (default=0.1)
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set. When fitting it is used to define the
        threshold on the decision function.

    n_jobs : int, optional (default = 1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`k_neighbors` and :meth:`kneighbors_graph` methods.


    Attributes
    ----------
    outlier_factor_ : numpy array, shape (n_samples,)
        The LOF of X. The lower, the more normal.
        Inliers tend to have a LOF score close to 1, while outliers tend
        to have a larger LOF score.

        The local outlier factor (LOF) of a sample captures its
        supposed 'degree of abnormality'.
        It is the average of the ratio of the local reachability density of
        a sample and those of its k-nearest neighbors.

    References
    ----------
    .. [1] Breunig, M. M., Kriegel, H. P., Ng, R. T., & Sander, J. (2000, May).
           LOF: identifying density-based local outliers. In ACM sigmod record.
    """
    def __init__(self, n_neighbors=5, algorithm='auto', leaf_size=30,
                 metric='minkowski', p=2, metric_params=None,
                 contamination=0.1, n_jobs=1):
        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size, metric=metric, p=p,
                          metric_params=metric_params, n_jobs=n_jobs)

        self.contamination = contamination

    def fit_predict(self, X, y=None):
        """Compute the local outlier factor (LOF) on X.
        Return the labels (1 inlier, -1 outlier) of X according to LOF score
        and the contamination parameter.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Factor
            wrt to the training samples.

        Returns
        -------
        is_inlier : array of shape (n_samples,)
            Returns 1 for anomalies/outliers and -1 for inliers.
        """

        return self.fit(X)._predict()

    def fit(self, X, y=None):
        """Fit the model using X as training data.

        Parameters
        ----------
        X : {array-like, sparse matrix, BallTree, KDTree}
            Training data. If array or matrix, shape [n_samples, n_features],
            or [n_samples, n_samples] if metric='precomputed'.

        Returns
        -------
        self : object
            Returns self.
        """
        if not (0. < self.contamination <= .5):
            raise ValueError("contamination must be in (0, 0.5]")

        super(LocalOutlierFactor, self).fit(X)

        self._distances_fit_X_, self._neighbors_indices_fit_X_ = (
            self.kneighbors(None))

        # Compute score over training samples to define threshold_:
        self.outlier_factor_ = self._local_outlier_factor(
            self._distances_fit_X_, self._neighbors_indices_fit_X_)

        self.threshold_ = -scoreatpercentile(
            self.outlier_factor_, 100. * (1. - self.contamination))

        # XXX may be optimized (if X is not None) by only computing it for
        # X_sub_samples = neighbors of neighbors of X ?
        return self

    def _predict(self, X=None):
        """Predict the labels (1 inlier, -1 outlier) of X according to LOF.

        If X is None, returns the same as fit_predict(X_train).
        This method allows to generalize prediction to new observations (not
        in the training set). As LOF originally does not deal with new data,
        this method is kept private.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Factor
            wrt to the training samples. If None, makes prediction on the
            training data without considering them as their own neighbors.

        Returns
        -------
        is_inlier : array of shape (n_samples,)
            Returns 1 for anomalies/outliers and -1 for inliers.
        """
        check_is_fitted(self, ["threshold_", "outlier_factor_",
                               "_distances_fit_X_",
                               "_neighbors_indices_fit_X_"])

        if X is not None:
            X = check_array(X, accept_sparse='csr')
            is_inlier = np.ones(X.shape[0], dtype=int)
            is_inlier[self.decision_function(X) <= self.threshold_] = -1
        else:
            is_inlier = np.ones(self._fit_X.shape[0], dtype=int)
            is_inlier[-self.outlier_factor_ <= self.threshold_] = -1

        return is_inlier

    def decision_function(self, X):
        """Opposite of the Local Outlier Factor of X (as bigger is better).

        The argument X is supposed to contain *new data*: if X contains a
        point from training, it consider the later in its own neighborhood.
        Also, the samples in X are not considered in the neighborhood of any
        point.
        The decision function on training data is available by considering the
        opposite of the outlier_factor_ attribute.

        Parameters
        ----------

        X : array-like, shape (n_samples, n_features)
            The query sample or samples to compute the Local Outlier Factor
            w.r.t. the training samples.

        Returns
        -------
        lof_scores : array of shape (n_samples,)
            The Local Outlier Factor of each input samples. The lower,
            the more abnormal.
        """
        check_is_fitted(self, ["threshold_", "outlier_factor_",
                               "_distances_fit_X_",
                               "_neighbors_indices_fit_X_"])

        X = check_array(X, accept_sparse='csr')

        distances_X, neighbors_indices_X = self.kneighbors(X)
        # as bigger is better:
        return -self._local_outlier_factor(distances_X, neighbors_indices_X)

    def _local_outlier_factor(self, distances_X, neighbors_indices_X):
        """Compute the local outlier factor (LOF)

        It is the average of the ratio of the local reachability density of
        a sample and those of its k-nearest neighbors.

        Parameters
        ----------
        distances_X : Array of shape (n_query, self.n_neighbors).
            Distances to the neighbors (in the training samples self._fit_X) of
            each query point to compute the LRD.

        neighbors_indices : Array of shape (n_query, self.n_neighbors)
            Neighbors indices (of each query point) among training samples
            self._fit_X.

        Returns
        -------
        lof : numpy array of shape (n_samples,)
            The LOF of X. The lower, the more normal.
        """
        n_samples = distances_X.shape[0]

        # Compute the local_reachability_density of samples X:
        X_lrd = self._local_reachability_density(distances_X,
                                                 neighbors_indices_X)
        lrd_ratios_array = np.zeros((n_samples, self.n_neighbors))

        # Avoid re-computing X_lrd if same parameters:
        if not (np.all(distances_X == self._distances_fit_X_) *
                np.all(self._neighbors_indices_fit_X_ == neighbors_indices_X)):
            lrd = self._local_reachability_density(
                self._distances_fit_X_, self._neighbors_indices_fit_X_)
        else:
            lrd = X_lrd

        lrd_ratios_array = lrd[neighbors_indices_X] / X_lrd[:, np.newaxis]

        return np.mean(lrd_ratios_array, axis=1)

    def _local_reachability_density(self, distances_X, neighbors_indices):
        """The local reachability density (LRD)

        The LRD of a sample is the inverse of the average reachability
        distance of its k-nearest neighbors.

        Parameters
        ----------
        distances_X : Array of shape (n_query, self.n_neighbors).
            Distances to the neighbors (in the training samples self._fit_X) of
            each query point to compute the LRD.

        neighbors_indices : Array of shape (n_query, self.n_neighbors)
            Neighbors indices (of each query point) among training samples
            self._fit_X.

        Returns
        -------
        local_reachability_density : array, shape (n_samples,)
            The local reachability density of each sample.
        """
        dist_k = self._distances_fit_X_[neighbors_indices,
                                        self.n_neighbors - 1]
        reach_dist_array = np.maximum(distances_X, dist_k)

        #  1e-10 to avoid `nan' when when nb of duplicates > n_neighbors:
        return 1. / (np.mean(reach_dist_array, axis=1) + 1e-10)
