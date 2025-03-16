# Authors: Erich Schubert <erich.schubert@tu-dortmund.de>
#          Nicolas Goix <nicolas.goix@telecom-paristech.fr>
#          Alexandre Gramfort <alexandre.gramfort@telecom-paristech.fr>
# License: BSD 3 clause

import warnings
from numbers import Real

import numpy as np

from ..base import OutlierMixin, _fit_context
from ..utils import check_array
from ..utils._param_validation import Interval, StrOptions
from ..utils.metaestimators import available_if
from ..utils.validation import check_is_fitted
from ._base import KNeighborsMixin, NeighborsBase

__all__ = ["NearestNeighborOutlierDetection"]


class NearestNeighborOutlierDetection(KNeighborsMixin, OutlierMixin, NeighborsBase):
    """Unsupervised Outlier Detection using the kNN distance.

    One of the earliest and simplest, yet effective, outlier detection
    methods is the distance (or the sum of distances) to the k-nearest neighbor.

    Parameters
    ----------
    n_neighbors : int, default=10
        Number of neighbors to use by default for :meth:`kneighbors` queries.
        If n_neighbors is larger than the number of samples provided,
        all samples will be used.

    weighting : {'max', 'sum'}, default='max'
        How the score is computed:

        - 'max' or 'knn' uses the distance to the kth nearest neighbor,
          the standard kNN outlier detection version of Ramaswamy et al.
        - 'sum', 'weight' or 'knnweight' uses the sum of distances to the
          k nearest neighbors, a variant due to Angiulli and Pizzuti

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, default='auto'
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDTree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, default=30
        Leaf is size passed to :class:`BallTree` or :class:`KDTree`. This can
        affect the speed of the construction and query, as well as the memory
        required to store the tree. The optimal value depends on the
        nature of the problem.

    metric : str or callable, default='minkowski'
        Metric to use for distance computation. Default is "minkowski", which
        results in the standard Euclidean distance when p = 2. See the
        documentation of `scipy.spatial.distance
        <https://docs.scipy.org/doc/scipy/reference/spatial.distance.html>`_ and
        the metrics listed in
        :class:`~sklearn.metrics.pairwise.distance_metrics` for valid metric
        values.

        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square during fit. X may be a :term:`sparse graph`, in which
        case only "nonzero" elements may be considered neighbors.

        If metric is a callable function, it takes two arrays representing 1D
        vectors as inputs and must return one value indicating the distance
        between those vectors. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

    p : float, default=2
        Parameter for the Minkowski metric from
        :func:`sklearn.metrics.pairwise_distances`. When p = 1, this
        is equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    metric_params : dict, default=None
        Additional keyword arguments for the metric function.

    contamination : float, default='0.01'
        The amount of contamination of the data set, i.e. the proportion
        of outliers in the data set. When fitting this is used to define the
        threshold on the scores of the samples.

        Because the kNN distances (in contrast to LOF and LoOP) do not
        provide a heuristic for this parameter, it defaults to 0.01
        corresponding to labeling 1% of the data as outliers.

        The contamination should usually be in the range (0, 0.5].

    novelty : bool, default=False
        By default, NearestNeighborOutlierDetection is only meant to be used
        for outlier detection (novelty=False). Set novelty to True if you want
        to use NearestNeighborOutlierDetection for novelty detection. In this
        case be aware that you should only use predict, decision_function and
        score_samples on new unseen data and not on the training set;
        and note that the results obtained this way may differ from the
        standard kNN results.

    n_jobs : int, default=None
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    negative_knn_distances_ : ndarray of shape (n_samples,)
        The negative kNN distances, for consistency with other parts
        of the scikit-learn API (e.g., OCSVM). A value close to 0 indicates
        a normal point, a negative value indicates outliers. The numeric range
        depends on the data set, preprocessing, and metric, and cannot be
        compared across data sets.

    n_neighbors_ : int
        The actual number of neighbors used for :meth:`kneighbors` queries.

    offset_ : float
        Offset used to obtain binary labels from the raw scores.
        Observations having negative_knn_distances_ smaller than
        `offset_` are detected as abnormal.
        The offset is defined such that the expected amount of outliers
        in training specified by the contamination parameter is flagged.

    effective_metric_ : str
        The effective metric used for the distance computation.

    effective_metric_params_ : dict
        The effective additional keyword arguments for the metric function.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    n_samples_fit_ : int
        It is the number of samples in the fitted data.

    See Also
    --------
    sklearn.neighbors.LocalOutlierFactory: Unsupervised Outlier Detection
        using the Local Outlier Factor.
    sklearn.neighbors.LocalOutlierProbabilities: Unsupervised Outlier Detection
        using Local Outlier Probabilities (LoOP).
    sklearn.svm.OneClassSVM: Unsupervised Outlier Detection using
        Support Vector Machine.

    References
    ----------
    .. [1] Ramaswamy, S., Rastogi, R. and Shim, K. (2000).
           Efficient algorithms for mining outliers from large data sets.
           In Proc. SIGMOD 2000, pp. 427–438.
       [2] Angiulli, F., & Pizzuti, C. (2002).
           Fast outlier detection in high dimensional spaces.
           In Proc. PKDD 2002, pp. 15-27.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.neighbors import NearestNeighborOutlierDetection
    >>> X = [[-1.1], [0.2], [101.1], [0.3]]
    >>> clf = NearestNeighborOutlierDetection(n_neighbors=2)
    >>> clf.fit_predict(X)
    array([ 1,  1, -1,  1])
    >>> -1 * clf.negative_knn_distances_
    array([ ?, ?, ?, ?])
    """

    _parameter_constraints: dict = {
        **NeighborsBase._parameter_constraints,
        "weighting": [StrOptions({"max", "knn", "sum", "weight", "knnweight"})],
        "contamination": [Interval(Real, 0, 0.5, closed="right")],
        "novelty": ["boolean"],
    }
    _parameter_constraints.pop("radius")

    def __init__(
        self,
        n_neighbors=20,
        weighting="max",
        *,
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        metric_params=None,
        contamination=0.01,
        novelty=False,
        n_jobs=None,
    ):
        super().__init__(
            n_neighbors=n_neighbors,
            algorithm=algorithm,
            leaf_size=leaf_size,
            metric=metric,
            p=p,
            metric_params=metric_params,
            n_jobs=n_jobs,
        )
        self.weighting = weighting
        self.contamination = contamination
        self.novelty = novelty

    def _check_novelty_fit_predict(self):
        if self.novelty:
            msg = (
                "fit_predict is not available when novelty=True. Use "
                "novelty=False if you want to predict on the training set."
            )
            raise AttributeError(msg)
        return True

    @available_if(_check_novelty_fit_predict)
    def fit_predict(self, X, y=None):
        """Fit the model to the training set X and return the labels.

        **Not available for novelty detection (when novelty is set to True).**
        Label is 0 for an inlier and -1 for an outlier according to the kNN
        score and the contamination parameter.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), default=None
            The query sample or samples to compute the kNN outlier scores
            w.r.t. the training samples.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        is_inlier : ndarray of shape (n_samples,)
            Returns -1 for anomalies/outliers and 1 for inliers.
        """

        # As fit_predict would be different from fit.predict, fit_predict is
        # only available for outlier detection (novelty=False)

        return self.fit(X)._predict()

    @_fit_context(
        # NearestNeighborOutlierDetection.metric is not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y=None):
        """Fit the kNN outlier scores detector from the training dataset.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features) or \
                (n_samples, n_samples) if metric='precomputed'
            Training data.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        self : NearestNeighborOutlierDetection
            The fitted nearest neighbor outlier detector.
        """
        self._fit(X)

        n_samples = self.n_samples_fit_
        if self.n_neighbors > n_samples:
            warnings.warn(
                "n_neighbors (%s) is greater than the "
                "total number of samples (%s). n_neighbors "
                "will be set to (n_samples - 1) for estimation."
                % (self.n_neighbors, n_samples)
            )
        self.n_neighbors_ = max(1, min(self.n_neighbors, n_samples - 1))

        distances_fit_X, _ = self.kneighbors(n_neighbors=self.n_neighbors_)

        if self._fit_X.dtype == np.float32:
            distances_fit_X = distances_fit_X.astype(self._fit_X.dtype, copy=False)

        # knn distances
        if self.weighting in ["max", "knn"]:
            self.negative_knn_distances_ = -1 * distances_fit_X[:, -1]
        # knn weight (c.f., Angiulli and Pizzuti)
        elif self.weighting in ["sum", "weight", "knnweight"]:
            self.negative_knn_distances_ = -1 * distances_fit_X.sum(axis=1)
        else:
            raise AttributeError("The weighting parameter was set to an invalid value.")

        # apply contamination to define a threshold for binary labeling
        self.offset_ = np.percentile(
            self.negative_knn_distances_, 100.0 * self.contamination
        )

        return self

    def _check_novelty_predict(self):
        if not self.novelty:
            msg = (
                "predict is not available when novelty=False, use "
                "fit_predict if you want to predict on training data. Use "
                "novelty=True if you want to use kNN for novelty detection "
                "and predict on new unseen data."
            )
            raise AttributeError(msg)
        return True

    @available_if(_check_novelty_predict)
    def predict(self, X=None):
        """Predict the labels (1 inlier, -1 outlier) of X according to kNN.

        **Only available for novelty detection (when novelty is set to True).**
        This method allows to generalize prediction to *new observations* (not
        in the training set). Note that the result of ``clf.fit(X)`` then
        ``clf.predict(X)`` with ``novelty=True`` may differ from the result
        obtained by ``clf.fit_predict(X)`` with ``novelty=False``.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The query sample or samples to compute the kNN outlier scores
            w.r.t. the training samples.

        Returns
        -------
        is_inlier : ndarray of shape (n_samples,)
            Returns -1 for anomalies/outliers and +1 for inliers.
        """
        return self._predict(X)

    def _predict(self, X=None):
        """Predict the labels (1 inlier, -1 outlier) of X according to kNN.

        If X is None, returns the same as fit_predict(X_train).

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), default=None
            The query sample or samples to compute the kNN outlier scores
            w.r.t. the training samples. If None, makes prediction on the
            training data without considering them as their own neighbors.

        Returns
        -------
        is_inlier : ndarray of shape (n_samples,)
            Returns -1 for anomalies/outliers and +1 for inliers.
        """
        check_is_fitted(self)

        if X is not None:
            X = check_array(X, accept_sparse="csr")
            is_inlier = np.ones(X.shape[0], dtype=int)
            is_inlier[self.decision_function(X) < 0] = -1
        else:
            is_inlier = np.ones(self.n_samples_fit_, dtype=int)
            is_inlier[self.negative_knn_distances_ < self.offset_] = -1

        return is_inlier

    def _check_novelty_decision_function(self):
        if not self.novelty:
            msg = (
                "decision_function is not available when novelty=False. "
                "Use novelty=True if you want to use kNN for novelty "
                "detection and compute decision_function for new unseen "
                "data. Note that the opposite kNN of the training samples "
                "is always available by considering the "
                "negative_knn_distances_ attribute."
            )
            raise AttributeError(msg)
        return True

    @available_if(_check_novelty_decision_function)
    def decision_function(self, X):
        """Shifted opposite of the kNN outlier scores of X.

        Bigger is better, i.e., large values correspond to inliers.

        **Only available for novelty detection (when novelty is set to True).**
        The shift offset allows a zero threshold for being an outlier.
        The argument X is supposed to contain *new data*: if X contains a
        point from training, it considers the later in its own neighborhood.
        Also, the samples in X are not considered in the neighborhood of any
        point.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The query sample or samples to compute the kNN outlier scores
            w.r.t. the training samples.

        Returns
        -------
        shifted_opposite_lof_scores : ndarray of shape (n_samples,)
            The shifted opposite of the kNN outlier scores of each input
            samples. The lower, the more abnormal. Negative scores represent
            outliers, positive scores represent inliers.
        """
        return self.score_samples(X) - self.offset_

    def _check_novelty_score_samples(self):
        if not self.novelty:
            msg = (
                "score_samples is not available when novelty=False. The "
                "scores of the training samples are always available "
                "through the negative_knn_distances_ attribute. Use "
                "novelty=True if you want to use kNN for novelty detection "
                "and compute score_samples for new unseen data."
            )
            raise AttributeError(msg)
        return True

    @available_if(_check_novelty_score_samples)
    def score_samples(self, X):
        """Opposite of the kNN outlier scores of X.

        It is the opposite as bigger is better, i.e., large values correspond
        to inliers.

        **Only available for novelty detection (when novelty is set to True).**
        The argument X is supposed to contain *new data*: if X contains a
        point from training, it considers the later in its own neighborhood.
        Also, the samples in X are not considered in the neighborhood of any
        point. Because of this, the scores obtained via ``score_samples`` may
        differ from the standard kNN scores.
        The standard kNN scores for the training data is available via the
        ``negative_knn_distances_`` attribute.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The query sample or samples to compute the kNN outlier scores
            w.r.t. the training samples.

        Returns
        -------
        opposite_lof_scores : ndarray of shape (n_samples,)
            The opposite of the kNN outlier scores of each input samples.
            The lower, the more abnormal.
        """
        check_is_fitted(self)
        X = check_array(X, accept_sparse="csr")

        distances_X, _ = self.kneighbors(X, n_neighbors=self.n_neighbors_)

        if X.dtype == np.float32:
            distances_X = distances_X.astype(X.dtype, copy=False)

        # knn distances
        if self.weighting in ["max", "knn"]:
            X_negative_knn_distances = -1 * distances_X[:, -1]
        # knn weight (c.f., Angiulli and Pizzuti)
        elif self.weighting in ["sum", "weight", "knnweight"]:
            X_negative_knn_distances = -1 * distances_X.sum(axis=1)
        else:
            raise AttributeError("The weighting parameter was set to an invalid value.")
        return X_negative_knn_distances

    def _more_tags(self):
        return {
            "preserves_dtype": [np.float64, np.float32],
        }
