# Authors: Pierrick Pochelu <pierrick.pochelu@gmail.com>

import numpy as np
import warnings

from ._base import NeighborsBase
from ._base import KNeighborsMixin
from ..base import OutlierMixin
from numbers import Real

from ..utils._param_validation import Interval, StrOptions
from ..utils.metaestimators import available_if
from ..utils.validation import check_is_fitted
from ..utils import check_array

__all__ = ["KNeighborsOutlierDetection"]


class KNeighborsOutlierDetection(KNeighborsMixin, OutlierMixin, NeighborsBase):
    """
        Unsupervised Outlier Detection using K Nearest Neighbors.

        Return the anomaly score of each sample using the KNN algorithm.
        The anomaly score of one sample is the sum of the distance to its k-th nearest neighbor.
        By comparing the scores, one can identify samples that have a substantially larger distance than the others.
        These are considered outliers.


        Parameters
        ----------
        n_neighbors : int, default=20
            Number of neighbors to use by default for :meth:`kneighbors` queries.
            If n_neighbors is larger than the number of samples provided,
            all samples will be used.

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

        p : int, default=2
            Parameter for the Minkowski metric from
            :func:`sklearn.metrics.pairwise.pairwise_distances`. When p = 1, this
            is equivalent to using manhattan_distance (l1), and euclidean_distance
            (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

        metric_params : dict, default=None
            Additional keyword arguments for the metric function.

        contamination : 'auto' or float, default='auto'
            The amount of contamination of the data set, i.e. the proportion
            of outliers in the data set. When fitting this is used to define the
            threshold on the scores of the samples.

            - if 'auto', the threshold is determined to fit the other anomaly detector implementations.
            - if a float, the contamination should be in the range (0, 0.5].

            .. versionchanged:: 0.22
               The default value of ``contamination`` changed from 0.1
               to ``'auto'``.


        n_jobs : int, default=None
            The number of parallel jobs to run for neighbors search.
            ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
            ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
            for more details.

        Attributes
        ----------
        cached_knn_distances : ndarray of shape (n_samples,)
            The higher, the more abnormal.
            Inliers tend to have a LOF score close to 1
            (``negative_outlier_factor_`` close to -1), while outliers tend to have
            a larger LOF score.

        n_neighbors : int
            The number of neighbors desired by the user. To avoid mistakes in the case where n_neighbors>length(X),
            we temporary considers that n_neighbors=length(X).

        offset_ : float
            Offset used to obtain binary labels from the raw scores.
            Observations having a negative_outlier_factor smaller than `offset_`
            are detected as abnormal.
            The offset is set to -1.5 (inliers score around -1), except when a
            contamination parameter different than "auto" is provided. In that
            case, the offset is defined in such a way we obtain the expected
            number of outliers in training.

            .. versionadded:: 0.20

        effective_metric_ : str
            The effective metric used for the distance computation.

        effective_metric_params_ : dict
            The effective additional keyword arguments for the metric function.

        n_features_in_ : int
            Number of features seen during :term:`fit`.

            .. versionadded:: 0.24

        feature_names_in_ : ndarray of shape (`n_features_in_`,)
            Names of features seen during :term:`fit`. Defined only when `X`
            has feature names that are all strings.

            .. versionadded:: 1.0

        n_samples_fit_ : int
            It is the number of samples in the fitted data.

        See Also
        --------

        sklearn.neighbors.LOF:  Unsupervised Outlier Detection using
        Local Outlier Factor

        sklearn.svm.OneClassSVM: Unsupervised Outlier Detection using
            Support Vector Machine.

        sklean.ensemble.IsolationForest: Unsupervised Outlier Detection using
        Isolation Forest



        References
        ----------
        .. [1] S. Ramaswamy, R. Rastogi, and K. Shim. Efficient algorithms for mining outliers from large data sets. In
    SIGMOD, pages 427–438, 2000.

        Examples
        --------
        >>> import numpy as np
        >>> from sklearn.neighbors import
        >>> X = [[0.1], [0.2], [0.3], [999.9], [0.4]]
        >>> clf = (n_neighbors=2)
        >>> clf.fit_predict(X)
        array([ 1, 1, 1, -1,  1])
        >>> clf.cached_knn_distances
        array([0.1, 0.1, 0.1, 999.5, 0.1])
    """

    _parameter_constraints: dict = {
        **NeighborsBase._parameter_constraints,
        "contamination": [
            StrOptions({"auto"}),
            Interval(Real, 0, 0.5, closed="right"),
        ],
    }
    _parameter_constraints.pop("radius")

    def __init__(
        self,
        n_neighbors=20,
        *,
        algorithm="auto",
        leaf_size=30,
        metric="minkowski",
        p=2,
        metric_params=None,
        contamination="auto",
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
        self.contamination = contamination
        self.cached_knn_distances = None

    def fit_predict(self, X, y=None):
        """Fit the model to the training set X and return the labels.
        fit_predict(X) is faster than computing separately fit(X) and predict(X). The distances computed during the fit(X) call are re-used in the predict(X) one.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), default=None
            The query sample or samples to compute the Local Outlier Factor
            w.r.t. to the training samples.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        is_inlier : ndarray of shape (n_samples,)
            Returns -1 for outliers (or "anomalies") and 1 for inliers.
        """
        return self.fit(X)._predict()

    def fit(self, X, y=None):
        """Fit the KNN detector from the training dataset.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features) or \
                (n_samples, n_samples) if metric='precomputed'
            Training data.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        self : KNeighborsOutlierDetection
            The fitted KNN outliers detector
        """
        self._validate_params()

        self._fit(X)

        # Calibrate offset_
        self.cached_knn_distances = self.knn_distance(X)
        if self.contamination == "auto":
            self.offset_ = np.quantile(self.cached_knn_distances, 1 - 0.1)
        else:
            self.offset_ = np.quantile(
                self.cached_knn_distances, 1 - self.contamination
            )

        return self

    def safe_k(self, desired_k, n_samples):
        """Compute the effective K value, taking into account the desired k value and the number of samples
         Parameters
        ----------
        desired_k : integer
        The hyperparameter K given by the user. According scikit-learn convention, K=2 means 1 neighbors and itself.

        n_samples : integer
        The number of samples feeding the KNN model.

        Returns
        -------
        safe_k : integer
        The safe K value which is correct
        """
        # edge cases
        if n_samples == 1:
            return 1

        # normal cases
        if desired_k > n_samples:
            warnings.warn(
                "n_neighbors (%s) is greater than the "
                "total number of samples (%s). n_neighbors "
                "will be set to n_samples for estimation."
                % (desired_k, n_samples)
            )
            safe_k = n_samples
        else:
            safe_k = desired_k
        return safe_k

    def predict(self, X=None):
        """
        Predict the labels of X. 1 for inliers and -1 for outliers.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)

        Returns
        -------
        is_inlier : ndarray of shape (n_samples,)
            Returns -1 for outliers (or "anomalies") and +1 for inliers.
        """
        return self._predict(X)

    def _predict(self, X=None):
        """Predict the labels (1 inlier, -1 outlier) of X according to LOF.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features), default=None
        If X is None we return previous computed samples. Caching allows to accelerate fit_predict(X).

        Returns
        -------
        is_inlier : ndarray of shape (n_samples,)
            Returns -1 for outliers (or "anomalies") and +1 for inliers.
        """
        check_is_fitted(self)

        if X is None:
            is_inlier = np.ones(self.cached_knn_distances.shape[0], dtype=int)
            is_inlier[self.cached_knn_distances > self.offset_] = -1
        else:
            X = check_array(X, accept_sparse="csr")
            is_inlier = np.ones(X.shape[0], dtype=int)
            is_inlier[self.knn_distance(X) > self.offset_] = -1

        return is_inlier

    def decision_function(self, X):
        """
        The anomaly score of X is centered.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        scores : ndarray of shape (n_samples,)
            The anomaly score of the input samples.
            The lower, the more abnormal. Negative scores represent outliers,
            positive scores represent inliers.
        """
        # 0 is the threshold value for being an outlier
        return self.offset_ - self.knn_distance(X)

    def knn_distance(self, X):
        """
        Opposite of the anomaly score defined in the original paper.

        The anomaly score of an input sample is computed as
        the mean anomaly score of the trees in the forest.

        The measure of normality of an observation given a tree is the depth
        of the leaf containing this observation, which is equivalent to
        the number of splittings required to isolate this point. In case of
        several observations n_left in the leaf, the average path length of
        a n_left samples isolation tree is added.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        Returns
        -------
        scores : ndarray of shape (n_samples,)
            The anomaly score of the input samples.
            The bigger, the more abnormal.
        """
        # code structure from ForestClassifier/predict_proba
        check_is_fitted(self)

        # Check data
        X = self._validate_data(X, accept_sparse="csr", reset=False)

        safe_k = self.safe_k(self.n_neighbors, X.shape[0])
        k_neighs_distances, k_neighs_identifiers = self.kneighbors(X, safe_k)
        sum_of_k_nearest_neighbours = np.sum(k_neighs_distances, axis=1)
        return sum_of_k_nearest_neighbours
