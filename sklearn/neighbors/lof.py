# Authors: Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import UnsupervisedMixin
from .base import LocalOutlierMixin

import warnings

from .base import NeighborsBase
from .base import KNeighborsMixin
from .base import UnsupervisedMixin
from ..base import OutlierMixin

from ..utils.validation import check_is_fitted

__all__ = ["LocalOutlierFactor"]


class LocalOutlierFactor(NeighborsBase, KNeighborsMixin, UnsupervisedMixin,
                         OutlierMixin, LocalOutlierMixin):
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
    density than their neighbors. These are considered outliers.

    Parameters
    ----------
    n_neighbors : int, optional (default=20)
        Number of neighbors to use by default for :meth:`kneighbors` queries.
        If n_neighbors is larger than the number of samples provided,
        all samples will be used.

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

        Valid values for metric are:

        - from scikit-learn: ['cityblock', 'cosine', 'euclidean', 'l1', 'l2',
          'manhattan']

        - from scipy.spatial.distance: ['braycurtis', 'canberra', 'chebyshev',
          'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski',
          'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao',
          'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean',
          'yule']

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

    contamination : float in (0., 0.5), optional (default=0.1)
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set. When fitting this is used to define the
        threshold on the decision function. If "auto", the decision function
        threshold is determined as in the original paper.

    n_jobs : int, optional (default=1)
        The number of parallel jobs to run for neighbors search.
        If ``-1``, then the number of jobs is set to the number of CPU cores.
        Affects only :meth:`kneighbors` and :meth:`kneighbors_graph` methods.


    Attributes
    ----------
    negative_outlier_factor_ : numpy array, shape (n_samples,)
        The opposite LOF of the training samples. The higher, the more normal.
        Inliers tend to have a LOF score close to 1 (negative_outlier_factor_
        close to -1), while outliers tend to have a larger LOF score.

        The local outlier factor (LOF) of a sample captures its
        supposed 'degree of abnormality'.
        It is the average of the ratio of the local reachability density of
        a sample and those of its k-nearest neighbors.

    n_neighbors_ : integer
        The actual number of neighbors used for :meth:`kneighbors` queries.

    offset_ : float
        Offset used to obtain binary labels from the raw scores.
        Observations having a negative_outlier_factor smaller than offset_ are
        detected as abnormal.
        The offset is set to -1.5 (inliers score around -1), except when a
        contamination parameter different than "auto" is provided. In that
        case, the offset is defined in such a way we obtain the expected
        number of outliers in training.

    References
    ----------
    .. [1] Breunig, M. M., Kriegel, H. P., Ng, R. T., & Sander, J. (2000, May).
           LOF: identifying density-based local outliers. In ACM sigmod record.
    """
    def __init__(self, n_neighbors=20, algorithm='auto', leaf_size=30,
                 metric='minkowski', p=2, metric_params=None,
                 contamination="legacy", n_jobs=1):
        super(LocalOutlierFactor, self).__init__(
              n_neighbors=n_neighbors,
              algorithm=algorithm,
              leaf_size=leaf_size, metric=metric, p=p,
              metric_params=metric_params, n_jobs=n_jobs)

        if contamination == "legacy":
            warnings.warn('default contamination parameter 0.1 will change '
                          'in version 0.22 to "auto". This will change the '
                          'predict method behavior.',
                          DeprecationWarning)
        self.contamination = contamination

    def fit_predict(self, X, y=None):
        """"Fits the model to the training set X and returns the labels
        (1 inlier, -1 outlier) on the training set according to the LOF score
        and the contamination parameter.


        Parameters
        ----------
        X : array-like, shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Factor
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
        """
        if self.contamination not in ["auto", "legacy"]:  # rm legacy in 0.22
            if not(0. < self.contamination <= .5):
                raise ValueError("contamination must be in (0, 0.5], "
                                 "got: %f" % self.contamination)

        super(LocalOutlierFactor, self).fit(X)

        LocalOutlierMixin.fit(self, X, y)

        return self

    def _predict(self, X=None):
        """Predict the labels (1 inlier, -1 outlier) of X according to LOF.

        If X is None, returns the same as fit_predict(X_train).
        This method allows to generalize prediction to new observations (not
        in the training set). As LOF originally does not deal with new data,
        this method is kept private. In particular, fit(X)._predict(X) is not
        the same as fit_predict(X).

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Factor
            w.r.t. to the training samples. If None, makes prediction on the
            training data without considering them as their own neighbors.

        Returns
        -------
        is_inlier : array, shape (n_samples,)
            Returns -1 for anomalies/outliers and +1 for inliers.
        """
        check_is_fitted(self, ["offset_", "negative_outlier_factor_",
                               "n_neighbors_", "_distances_fit_X_"])

        is_inlier = LocalOutlierMixin._score_samples(self, X)

        return is_inlier

    def _decision_function(self, X):
        """Shifted opposite of the Local Outlier Factor of X

        Bigger is better, i.e. large values correspond to inliers.

        The shift offset allows a zero threshold for being an outlier.
        The argument X is supposed to contain *new data*: if X contains a
        point from training, it consider the later in its own neighborhood.
        Also, the samples in X are not considered in the neighborhood of any
        point.
        This method is kept private as the predict method is.
        The decision function on training data is available by considering the
        the negative_outlier_factor_ attribute.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The query sample or samples to compute the Local Outlier Factor
            w.r.t. the training samples.

        Returns
        -------
        shifted_opposite_lof_scores : array, shape (n_samples,)
            The shifted opposite of the Local Outlier Factor of each input
            samples. The lower, the more abnormal. Negative scores represent
            outliers, positive scores represent inliers.
        """
        return self._score_samples(X) - self.offset_

    def _score_samples(self, X):
        """Opposite of the Local Outlier Factor of X (as bigger is
        better, i.e. large values correspond to inliers).

        The argument X is supposed to contain *new data*: if X contains a
        point from training, it consider the later in its own neighborhood.
        Also, the samples in X are not considered in the neighborhood of any
        point.
        This method is kept private as the predict method is.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The query sample or samples to compute the Local Outlier Factor
            w.r.t. the training samples.

        Returns
        -------
        opposite_lof_scores : array, shape (n_samples,)
            The opposite of the Local Outlier Factor of each input samples.
            The lower, the more abnormal.
        """
        check_is_fitted(self, ["offset_", "negative_outlier_factor_",
                               "_distances_fit_X_"])

        return LocalOutlierMixin._decision_function(self, X)


