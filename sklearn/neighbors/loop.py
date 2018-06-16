# Author: Valentino Constantinou <vc@valentino.io>
# License: BSD 3 clause

from math import erf
import numpy as np

from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import UnsupervisedMixin
from .base import LocalOutlierMixin

from ..utils.validation import check_is_fitted

__all__ = ["LocalOutlierProbability"]


class LocalOutlierProbability(NeighborsBase, KNeighborsMixin, UnsupervisedMixin, LocalOutlierMixin):
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
    .. [1] Breunig, M. M., Kriegel, H.-P., Ng, R. T., & Sander, J. (2000, May).
           LOF: identifying density-based local outliers. In ACM sigmod record.
    .. [2] Kriegel H.-P., Kroeger P., Schubert E., Zimek A. LoOP: Local Outlier
           Probabilities. 18th ACM conference on Information and knowledge
           management, CIKM (2009).
    .. [3] Goldstein M., Uchida S. A Comparative Evaluation of Unsupervised Anomaly
           Detection Algorithms for Multivariate Data. PLoS ONE 11(4): e0152173 (2016).
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
        set distance and the norm_factor parameter.


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
            raise ValueError("norm_factor must be in (0, 1.0]")

        super(LocalOutlierProbability, self).fit(X)

        LocalOutlierMixin.fit(self, X, y, mode='loop')

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

        is_inlier = LocalOutlierMixin._assign_label(self, X, mode='loop')
        if X is None:
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
        check_is_fitted(self,
                        ["norm_factor", "n_neighbors", "_distances_fit_X_", "negative_local_outlier_probability_"])

        return LocalOutlierMixin._decision_function(self, X, mode='loop')
