# Author: Valentino Constantinou <vc@valentino.io>
# License: BSD 3 clause

from math import erf
import numpy as np
from warnings import warn
from scipy.stats import scoreatpercentile

from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import UnsupervisedMixin

from ..utils.validation import check_is_fitted
from ..utils import check_array

__all__ = ["LocalOutlierProbability"]


class LocalOutlierProbability(NeighborsBase, KNeighborsMixin, UnsupervisedMixin):
    """Unsupervised Outlier Detection using Local Outlier Probabilities (LoOP)

    The outlier score of each sample is called the Local Outlier Probability.
    It measures the local deviation of density of a given sample with respect
    to its neighbors as Local Outlier Factor (LOF), but provides normalized
    outlier scores in the range [0,1]. These outlier scores are directly
    interpretable as a probability of an object being an outlier.
    Like LOF, it is local in that the anomaly score depends on how isolated
    the sample is with respect to the surrounding neighborhood. Locality is
    given by k-nearest neighbors, whose distance is used to calculate the
    probabilistic set distance, which is used to estimate the local
    density. By comparing the local density of a sample to the local densities
    of its neighbors, one can identify samples that lie in regions of lower
    density compared to their neighbors and thus identify samples that may be
    outliers according to their Local Outlier Probability.

    Parameters
    ----------
    n_neighbors : int, optional (default=20)
        Number of neighbors to use by default. If n_neighbors is larger than
        the number of samples provided, all samples will be used.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default=30)
        Leaf size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    metric : string or callable, default 'euclidean'
        Metric used for the distance computation.

        If 'precomputed', the training input X is expected to be a distance
        matrix.

        If metric is a callable function, it is called on each
        pair of instances (rows) and the resulting value recorded. The callable
        should take two arrays as input and return one value indicating the
        distance between them. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto',
          'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
          'sqeuclidean', 'yule']

        See the documentation for scipy.spatial.distance for details on these
        metrics:
        http://docs.scipy.org/doc/scipy/reference/spatial.distance.html

    p : integer, optional (default=2)
        Parameter for the Minkowski metric from
        :func:`sklearn.metrics.pairwise.pairwise_distances`. When p = 1, this
        is equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, optional (default=None)
        Additional keyword arguments for the metric function.

    norm_factor : float in (0, 1.0], optional (default=0.95)
        A normalization factor that gives control over the density approximation.
        The norm_factor does not affect the ranking of outliers. norm_factor is
        defined according to the triple sigma rule. The default value of 0.997
        corresponds to a threshold of 3 standard deviations from the mean
        (0.95 for 2 standard deviations, 0.68 for 1 standard deviation).

    n_jobs : int, optional (default=1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.

    Attributes
    ----------
    negative_local_outlier_probability_ : numpy array, shape (n_samples,)
        The negative LoOP of the training samples. The higher, the more normal.
        Outliers tend to have a LoOP score close to -1, while normal values tend
        to have a LoOP score at or close to 0.

        The local outlier probability (LoOP) of a sample captures its
        supposed 'degree of abnormality'.
        It is the probability that a sample is an outlier with respect to
        its n_neighbors.

    n_neighbors_ : integer
        The actual number of neighbors used for :meth:`kneighbors` queries.

    References
    ----------
    .. [1] Breunig, M. M., Kriegel, H. P., Ng, R. T., & Sander, J. (2000, May).
           LOF: identifying density-based local outliers. In ACM sigmod record.
    .. [2] Kriegel H., Kröger P., Schubert E., Zimek A. LoOP: Local Outlier
           Probabilities. 18th ACM conference on Information and knowledge
           management, CIKM (2009).
    .. [3] Goldstein M., Uchida S. A Comparative Evaluation of Unsupervised Anomaly
           Detection Algorithms for Multivariate Data. PLoS ONE 11(4): e0152173 (2016)
    """

    def __init__(self, n_neighbors=20, norm_factor=0.95, algorithm='auto', leaf_size=30,
                 metric='euclidean', p=2, metric_params=None, n_jobs=1):

        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size,
                          metric=metric, p=p,
                          metric_params=metric_params,
                          n_jobs=n_jobs)

        self.norm_factor = norm_factor

    def fit_predict(self, X, y=None):
        """"Fits the model to the training set X and returns the labels
        (1 inlier, -1 outlier) on the training set according to the probabilistic
        set distance and the extent parameter.


        Parameters
        ----------
        X : array-like, shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Probability
            w.r.t. to the training samples.

        Returns
        -------
        is_inlier : array, shape (n_samples,)
            Returns -1 for anomalies/outliers and 1 for inliers.
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
        if not (0. < self.norm_factor <= 1.):
            raise ValueError("extent must be in (0, 1.0]")

        super(LocalOutlierProbability, self).fit(X)

        n_samples = self._fit_X.shape[0]
        if self.n_neighbors > n_samples:
            warn("n_neighbors (%s) is greater than the "
                 "total number of samples (%s). n_neighbors "
                 "will be set to (n_samples - 1) for estimation."
                 % (self.n_neighbors, n_samples))
        self.n_neighbors_ = max(1, min(self.n_neighbors, n_samples - 1))

        self._distances_fit_X_, _neighbors_indices_fit_X_ = (
            self.kneighbors(None, n_neighbors=self.n_neighbors_))

        if not self._distances_fit_X_.any():
            warn('Neighborhood distances all zero. Try using a larger value for n_neighbors.', RuntimeWarning)

        # Compute loop score
        ssd = self._ssd(self._distances_fit_X_)
        self._standard_distances_fit_X_ = self._standard_distances(ssd)
        self._prob_distances_fit_X_ = self._prob_distances(self._standard_distances_fit_X_)
        self._prob_distances_ev_fit_X_ = self._prob_distances_ev(self._prob_distances_fit_X_)
        prob_local_outlier_factors = self._prob_local_outlier_factors(self._prob_distances_fit_X_,
                                                                      self._prob_distances_ev_fit_X_)
        self._prob_local_outlier_factors_ev_fit_X_ = self._prob_local_outlier_factors_ev(prob_local_outlier_factors)
        norm_prob_local_outlier_factors = self._norm_prob_local_outlier_factors(
            self._prob_local_outlier_factors_ev_fit_X_)
        self.negative_local_outlier_probability_ = self._neg_local_outlier_probability(prob_local_outlier_factors,
                                                                                       norm_prob_local_outlier_factors)

        # Compute the local reachability density using the probabilistic set distance to define threshold_:
        self._lrd = self._local_reachability_density(
            self._prob_distances_fit_X_, _neighbors_indices_fit_X_)

        lrd_ratios_array = (self._lrd[_neighbors_indices_fit_X_] /
                            self._lrd[:, np.newaxis])

        self.lrd_ratios_ = -np.mean(lrd_ratios_array, axis=1)

        self.threshold_ = scoreatpercentile(
            self.lrd_ratios_, 100. * (1. - self.norm_factor))

        return self

    def _predict(self, X=None):
        """Predict the labels (1 inlier, -1 outlier) of X according to LoOP.

        If X is None, returns the same as fit_predict(X_train).
        This method allows to generalize prediction to new observations (not
        in the training set). As LoOP originally does not deal with new data,
        this method is kept private.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Probability
            w.r.t. to the training samples. If None, makes prediction on the
            training data without considering them as their own neighbors.

        Returns
        -------
        is_inlier : array, shape (n_samples,)
            Returns -1 for anomalies/outliers and +1 for inliers.
        """
        check_is_fitted(self, ["threshold_", "negative_local_outlier_probability_",
                               "n_neighbors_", "_distances_fit_X_"])

        if X is not None:
            X = check_array(X, accept_sparse='csr')
            is_inlier = np.ones(X.shape[0], dtype=int)
            is_inlier[self._decision_function(X) <= self.threshold_] = -1
        else:
            is_inlier = np.ones(self._fit_X.shape[0], dtype=int)
            is_inlier[self.lrd_ratios_ <= self.threshold_] = -1

            return is_inlier

    @staticmethod
    def _ssd(distances_X):
        """
        Calculates the sum of square distance of X. As this method is used in the
        calculation of the local probability of each sample, this method is kept
        private.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The pairwise distances between each sample in the input data.

        Returns
        -------
        ssd_array : array, shape (n_samples,)
            Returns the sum of square distance of X as an array of
            repeated values.
        """
        ssd_array = np.sum(np.power(distances_X, 2), axis=1)

        return ssd_array

    def _standard_distances(self, ssd_X):
        """
        Calculates the standard distance of each sample in the input data. As
        this method is used in the calculation of the local probability of each sample,
        this method is kept private.

        Parameters
        ----------
        ssd_X : array-like, shape (n_samples,)
            The sum of square distance of X as an array of repeated values.

        Returns
        -------
        standard_distances : array, shape (n_samples,)
            The standard distance of each sample.
        """
        standard_distances = np.sqrt(np.divide(ssd_X, self.n_neighbors))

        return standard_distances

    def _prob_distances(self, standard_distances_X):
        """
        Calculates the probabilistic distance of each sample in the input data.
        As this method is used in the calculation of the local probability of each sample,
        this method is kept private.

        Parameters
        ----------
        standard_distances_X : array-like, shape (n_samples,)
            The standard distance of each sample.

        Returns
        -------
        prob_distances : array, shape(n_samples,)
            The probabilistic distances of each sample.
        """
        prob_distances = (self.norm_factor * standard_distances_X)

        return prob_distances

    @staticmethod
    def _prob_distances_ev(probabilistic_set_distances_X):
        """
        Calculates the expected value of the probabilistic distance. As this method is used in the
        calculation of the local probability of each sample, this method is kept private.

        Parameters
        ----------
        probabilistic_distances_X : array-like, shape (n_samples,)
            The probabilistic set distances of X.

        Returns
        -------
        prob_distance_ev_array : array, shape (n_samples,)
            The expected value of the probabilistic distance of X as an array
            of repeated values.
        """
        prob_set_distance_ev_array = np.array(
            [np.mean(probabilistic_set_distances_X[~np.isinf(probabilistic_set_distances_X)])] *
            probabilistic_set_distances_X.shape[0])

        return prob_set_distance_ev_array

    @staticmethod
    def _prob_local_outlier_factors(probabilistic_set_distances_X, probabilistic_set_distances_ev_X):
        """
        Calculates the probabilistic local outlier factor of each sample. As this method is used in the
        calculation of the local probability of each sample, this method is kept
        private.

        Parameters
        ----------
        probabilistic_distances_X : array-like, shape (n_samples,)
            The probabilistic distances of X.

        probabilistic_distances_ev_X : array-like, shape (n_samples,)
            The expected value of the probabilistic distances of X.

        Returns
        -------
        prob_local_outlier_factors : array, shape (n_samples,)
            The probabilistic local outlier factor of each sample.
        """
        prob_local_outlier_factors = (probabilistic_set_distances_X / probabilistic_set_distances_ev_X) - 1.

        return prob_local_outlier_factors

    @staticmethod
    def _prob_local_outlier_factors_ev(probabilistic_local_outlier_factors_X):
        """
        Calculates the expected value of the probabilistic local outlier factor.
        As this method is used in the calculation of the local probability of each sample,
        this method is kept private.

        Parameters
        ----------
        probabilistic_local_outlier_factors_X : array-like, shape (n_samples,)
            The probabilistic local outlier factor of each sample.

        Returns
        -------
        prob_local_outlier_factors_ev : array, shape (n_samples,)
            The expected value of the probabilistic local outlier factor of X as an array
            of repeated values.
        """
        prob_local_outlier_factors_ev = np.array([np.sum(
            np.power(probabilistic_local_outlier_factors_X[~np.isinf(probabilistic_local_outlier_factors_X)], 2) /
            probabilistic_local_outlier_factors_X.shape[0])] * probabilistic_local_outlier_factors_X.shape[0])

        return prob_local_outlier_factors_ev

    def _norm_prob_local_outlier_factors(self, prob_local_outlier_factors_ev_X):
        """
        Calculates the normalized probabilistic local outlier factor. As this method is used in the
        calculation of the local probability of each sample, this method is kept
        private.

        Parameters
        ----------
        prob_local_outlier_factors_ev_X : array-like, shape (n_samples,)
            The expected value of the probabilistic local outlier factor of X as an
            array of repeated values.

        Returns
        -------
        norm_prob_local_outlier_factors : array, shape (n_samples,)
            The normalized probabilistic local outlier factor of each sample.
        """
        norm_prob_local_outlier_factors = self.norm_factor * np.sqrt(prob_local_outlier_factors_ev_X)

        return norm_prob_local_outlier_factors

    @staticmethod
    def _neg_local_outlier_probability(probabilistic_local_outlier_factors_X,
                                       normalized_probabilistic_local_outlier_factors):
        """
        Calculates the negative local outlier probability of each sample. As this method is used in the
        calculation of the local probability of each sample, this method is kept
        private.

        Parameters
        ----------
        probabilistic_local_outlier_factors_X : array-like, shape (n_samples,)
            The probabilistic local outlier factor of each sample in X.

        normalized_probabilistic_local_outlier_factors : array-like, shape (n_samples,)
            The normalized probabilistic local outlier factor of each sample in X.

        Returns
        -------
        negative_local_outlier_probability : array, shape (n_samples,)
            The negative local outlier probability of each sample in X.

        """
        erf_vec = np.vectorize(erf)
        negative_local_outlier_probability = -1. * np.maximum(0., erf_vec(
            probabilistic_local_outlier_factors_X / (normalized_probabilistic_local_outlier_factors * np.sqrt(2.))))

        return negative_local_outlier_probability

    def _decision_function(self, X):
        """Opposite of the local reachability density of X (as bigger is better,
        i.e. large values correspond to more normal values).

        The argument X is supposed to contain *new data*: if X contains a
        point from training, it consider the later in its own neighborhood.
        Also, the samples in X are not considered in the neighborhood of any
        point.
        The decision function on training data is available by considering the
        opposite of the lrd_ratio_ attribute.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The query sample or samples to compute the Local Outlier Probability
            w.r.t. the training samples.

        Returns
        -------
        lrd_ratio : array, shape (n_samples,)
            The opposite of the local reachability density of each input samples.
            The lower, the more abnormal.
        """
        check_is_fitted(self, ["extent", "n_neighbors", "_distances_fit_X_", "negative_local_outlier_probability_"])

        X = check_array(X, accept_sparse='csr')

        distances_X, neighbors_indices_X = (
            self.kneighbors(X, n_neighbors=self.n_neighbors_))

        ssd = self._ssd(distances_X)
        standard_distances_X = self._standard_distances(ssd)
        prob_distances = self._prob_distances(standard_distances_X)

        X_lrd = self._local_reachability_density(prob_distances,
                                                 neighbors_indices_X)

        lrd_ratios_array = (self._lrd[neighbors_indices_X] /
                            X_lrd[:, np.newaxis])

        # as bigger is better:
        return -np.mean(lrd_ratios_array, axis=1)

    def _local_reachability_density(self, distances_X, neighbors_indices):
        """The local reachability density (LRD)

        The LRD of a sample is the inverse of the probabilistic
        distance of its k-nearest neighbors.

        Parameters
        ----------
        distances_X : array, shape (n_query, self.n_neighbors)
            Probabilistic distances to the neighbors (in the training samples `self._fit_X`)
            of each query point to compute the LRD.

        neighbors_indices : array, shape (n_query, self.n_neighbors)
            Neighbors indices (of each query point) among training samples
            self._fit_X.

        Returns
        -------
        local_reachability_density : array, shape (n_samples,)
            The local reachability density of each sample.
        """
        self._prob_set_distances_fit_X_tiled_ = np.tile(self._prob_distances_fit_X_, (self.n_neighbors_, 1)).T
        dist_k = self._prob_set_distances_fit_X_tiled_[neighbors_indices, self.n_neighbors_ - 1]
        prob_set_distances_tiled = np.tile(distances_X, (self.n_neighbors_, 1)).T

        reach_dist_array = np.maximum(prob_set_distances_tiled, dist_k)

        #  1e-10 to avoid `nan' when nb of duplicates > n_neighbors_:
        return 1. / (np.mean(reach_dist_array, axis=1) + 1e-10)
