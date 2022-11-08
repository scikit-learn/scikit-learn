"""Nearest Neighbor Classification"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck
#          Multi-output support by Arnaud Joly <a.joly@ulg.ac.be>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam
from numbers import Integral

import numpy as np
from ..utils.fixes import _mode
from ..utils.extmath import weighted_mode
from ..utils.validation import _is_arraylike, _num_samples

import warnings
from ._base import _get_weights
from ._base import NeighborsBase, KNeighborsMixin, RadiusNeighborsMixin
from ..base import ClassifierMixin
from ..utils._param_validation import StrOptions
from .._engine import get_engine_classes


class KNeighborsClassifierCythonEngine:
    def __init__(self, estimator):
        self.estimator = estimator

    def accepts(self, X, y=None):
        # The default engine accepts everything
        return True

    def fit(self, X, y=None):
        return self.estimator._fit(X, y)

    def predict(self, X):
        if self.estimator.weights == "uniform":
            # In that case, we do not need the distances to perform
            # the weighting so we do not compute them.
            neigh_ind = self.estimator.kneighbors(X, return_distance=False)
            neigh_dist = None
        else:
            neigh_dist, neigh_ind = self.estimator.kneighbors(X)

        classes_ = self.estimator.classes_
        _y = self.estimator._y
        if not self.estimator.outputs_2d_:
            _y = self.estimator._y.reshape((-1, 1))
            classes_ = [self.estimator.classes_]

        n_outputs = len(classes_)
        n_queries = _num_samples(X)
        weights = _get_weights(neigh_dist, self.estimator.weights)

        y_pred = np.empty((n_queries, n_outputs), dtype=classes_[0].dtype)
        for k, classes_k in enumerate(classes_):
            if weights is None:
                mode, _ = _mode(_y[neigh_ind, k], axis=1)
            else:
                mode, _ = weighted_mode(_y[neigh_ind, k], weights, axis=1)

            mode = np.asarray(mode.ravel(), dtype=np.intp)
            y_pred[:, k] = classes_k.take(mode)

        if not self.estimator.outputs_2d_:
            y_pred = y_pred.ravel()

        return y_pred

    def predict_proba(self, X):
        if self.estimator.weights == "uniform":
            # In that case, we do not need the distances to perform
            # the weighting so we do not compute them.
            neigh_ind = self.estimator.kneighbors(X, return_distance=False)
            neigh_dist = None
        else:
            neigh_dist, neigh_ind = self.estimator.kneighbors(X)

        classes_ = self.estimator.classes_
        _y = self.estimator._y
        if not self.estimator.outputs_2d_:
            _y = self.estimator._y.reshape((-1, 1))
            classes_ = [self.estimator.classes_]

        n_queries = _num_samples(X)

        weights = _get_weights(neigh_dist, self.estimator.weights)
        if weights is None:
            weights = np.ones_like(neigh_ind)

        all_rows = np.arange(n_queries)
        probabilities = []
        for k, classes_k in enumerate(classes_):
            pred_labels = _y[:, k][neigh_ind]
            proba_k = np.zeros((n_queries, classes_k.size))

            # a simple ':' index doesn't work right
            for i, idx in enumerate(pred_labels.T):  # loop is O(n_neighbors)
                proba_k[all_rows, idx] += weights[:, i]

            # normalize 'votes' into real [0,1] probabilities
            normalizer = proba_k.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            proba_k /= normalizer

            probabilities.append(proba_k)

        if not self.estimator.outputs_2d_:
            probabilities = probabilities[0]

        return probabilities

    def kneighbors(self, X=None, n_neighbors=None, return_distance=True):
        return self.estimator._kneighbors(
            X=X, n_neighbors=n_neighbors, return_distance=return_distance
        )


class KNeighborsClassifier(KNeighborsMixin, ClassifierMixin, NeighborsBase):
    """Classifier implementing the k-nearest neighbors vote.

    Read more in the :ref:`User Guide <classification>`.

    Parameters
    ----------
    n_neighbors : int, default=5
        Number of neighbors to use by default for :meth:`kneighbors` queries.

    weights : {'uniform', 'distance'}, callable or None, default='uniform'
        Weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

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
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    p : int, default=2
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

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

    metric_params : dict, default=None
        Additional keyword arguments for the metric function.

    n_jobs : int, default=None
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.
        Doesn't affect :meth:`fit` method.

    Attributes
    ----------
    classes_ : array of shape (n_classes,)
        Class labels known to the classifier

    effective_metric_ : str or callble
        The distance metric used. It will be same as the `metric` parameter
        or a synonym of it, e.g. 'euclidean' if the `metric` parameter set to
        'minkowski' and `p` parameter set to 2.

    effective_metric_params_ : dict
        Additional keyword arguments for the metric function. For most metrics
        will be same with `metric_params` parameter, but may also contain the
        `p` parameter value if the `effective_metric_` attribute is set to
        'minkowski'.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_samples_fit_ : int
        Number of samples in the fitted data.

    outputs_2d_ : bool
        False when `y`'s shape is (n_samples, ) or (n_samples, 1) during fit
        otherwise True.

    See Also
    --------
    RadiusNeighborsClassifier: Classifier based on neighbors within a fixed radius.
    KNeighborsRegressor: Regression based on k-nearest neighbors.
    RadiusNeighborsRegressor: Regression based on neighbors within a fixed radius.
    NearestNeighbors: Unsupervised learner for implementing neighbor searches.

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    .. warning::

       Regarding the Nearest Neighbors algorithms, if it is found that two
       neighbors, neighbor `k+1` and `k`, have identical distances
       but different labels, the results will depend on the ordering of the
       training data.

    https://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> neigh = KNeighborsClassifier(n_neighbors=3)
    >>> neigh.fit(X, y)
    KNeighborsClassifier(...)
    >>> print(neigh.predict([[1.1]]))
    [0]
    >>> print(neigh.predict_proba([[0.9]]))
    [[0.666... 0.333...]]
    """

    _parameter_constraints: dict = {**NeighborsBase._parameter_constraints}
    _parameter_constraints.pop("radius")
    _parameter_constraints.update(
        {"weights": [StrOptions({"uniform", "distance"}), callable, None]}
    )

    def __init__(
        self,
        n_neighbors=5,
        *,
        weights="uniform",
        algorithm="auto",
        leaf_size=30,
        p=2,
        metric="minkowski",
        metric_params=None,
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
        self.weights = weights

    def fit(self, X, y):
        """Fit the k-nearest neighbors classifier from the training dataset.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features) or \
                (n_samples, n_samples) if metric='precomputed'
            Training data.

        y : {array-like, sparse matrix} of shape (n_samples,) or \
                (n_samples, n_outputs)
            Target values.

        Returns
        -------
        self : KNeighborsClassifier
            The fitted k-nearest neighbors classifier.
        """
        self._validate_params()

        engine = self._get_engine(X, y)

        if hasattr(engine, "pre_fit"):
            X, y, sample_weight = engine.pre_fit(
                X,
                y=y,
            )

        engine.fit(
            X,
            y=y,
        )

        if hasattr(engine, "post_fit"):
            engine.post_fit(
                X,
                y=y,
            )

        return self

    def predict(self, X):
        """Predict the class labels for the provided data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_queries, n_features), \
                or (n_queries, n_indexed) if metric == 'precomputed'
            Test samples.

        Returns
        -------
        y : ndarray of shape (n_queries,) or (n_queries, n_outputs)
            Class labels for each data sample.
        """
        engine = self._get_engine(X)

        if hasattr(engine, "pre_predict"):
            X = engine.pre_predict(
                X,
            )

        y_pred = engine.predict(X)

        if hasattr(engine, "post_predict"):
            engine.post_predict(
                X,
            )

        return y_pred

    def predict_proba(self, X):
        """Return probability estimates for the test data X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_queries, n_features), \
                or (n_queries, n_indexed) if metric == 'precomputed'
            Test samples.

        Returns
        -------
        p : ndarray of shape (n_queries, n_classes), or a list of n_outputs \
                of such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are ordered
            by lexicographic order.
        """
        engine = self._get_engine(X)

        if hasattr(engine, "pre_predict_proba"):
            X = engine.pre_predict_proba(
                X,
            )

        probabilities = engine.predict_proba(X)

        if hasattr(engine, "post_predict_proba"):
            engine.post_predict_proba(
                X,
            )

        return probabilities

    def kneighbors(self, X=None, n_neighbors=None, return_distance=True):
        """Find the K-neighbors of a point.

        Returns indices of and distances to the neighbors of each point.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_queries, n_features), \
            or (n_queries, n_indexed) if metric == 'precomputed', default=None
            The query point or points.
            If not provided, neighbors of each indexed point are returned.
            In this case, the query point is not considered its own neighbor.

        n_neighbors : int, default=None
            Number of neighbors required for each sample. The default is the
            value passed to the constructor.

        return_distance : bool, default=True
            Whether or not to return the distances.

        Returns
        -------
        neigh_dist : ndarray of shape (n_queries, n_neighbors)
            Array representing the lengths to points, only present if
            return_distance=True.

        neigh_ind : ndarray of shape (n_queries, n_neighbors)
            Indices of the nearest points in the population matrix.

        Examples
        --------
        In the following example, we construct a NearestNeighbors
        class from an array representing our data set and ask who's
        the closest point to [1,1,1]

        >>> samples = [[0., 0., 0.], [0., .5, 0.], [1., 1., .5]]
        >>> from sklearn.neighbors import NearestNeighbors
        >>> neigh = NearestNeighbors(n_neighbors=1)
        >>> neigh.fit(samples)
        NearestNeighbors(n_neighbors=1)
        >>> print(neigh.kneighbors([[1., 1., 1.]]))
        (array([[0.5]]), array([[2]]))

        As you can see, it returns [[0.5]], and [[2]], which means that the
        element is at distance 0.5 and is the third element of samples
        (indexes start at 0). You can also query for multiple points:

        >>> X = [[0., 1., 0.], [1., 0., 1.]]
        >>> neigh.kneighbors(X, return_distance=False)
        array([[1],
               [2]]...)
        """
        engine = self._get_engine(X)

        if hasattr(engine, "pre_kneighbors"):
            X = engine.pre_kneighbors(
                X=X, n_neighbors=n_neighbors, return_distance=return_distance
            )

        probabilities = engine.kneighbors(
            X=X, n_neighbors=n_neighbors, return_distance=return_distance
        )

        if hasattr(engine, "post_kneighbors"):
            engine.post_kneighbors(
                X=X, n_neighbors=n_neighbors, return_distance=return_distance
            )

        return probabilities

    def _more_tags(self):
        return {"multilabel": True}

    def _get_engine(self, X, y=None, sample_weight=None):
        for engine_class in get_engine_classes(
            "kneigborsclassifier", default=KNeighborsClassifierCythonEngine
        ):
            engine = engine_class(self)
            if engine.accepts(X, y=y):
                return engine


class RadiusNeighborsClassifier(RadiusNeighborsMixin, ClassifierMixin, NeighborsBase):
    """Classifier implementing a vote among neighbors within a given radius.

    Read more in the :ref:`User Guide <classification>`.

    Parameters
    ----------
    radius : float, default=1.0
        Range of parameter space to use by default for :meth:`radius_neighbors`
        queries.

    weights : {'uniform', 'distance'}, callable or None, default='uniform'
        Weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

        Uniform weights are used by default.

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
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    p : int, default=2
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

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

    outlier_label : {manual label, 'most_frequent'}, default=None
        Label for outlier samples (samples with no neighbors in given radius).

        - manual label: str or int label (should be the same type as y)
          or list of manual labels if multi-output is used.
        - 'most_frequent' : assign the most frequent label of y to outliers.
        - None : when any outlier is detected, ValueError will be raised.

    metric_params : dict, default=None
        Additional keyword arguments for the metric function.

    n_jobs : int, default=None
        The number of parallel jobs to run for neighbors search.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    Attributes
    ----------
    classes_ : ndarray of shape (n_classes,)
        Class labels known to the classifier.

    effective_metric_ : str or callable
        The distance metric used. It will be same as the `metric` parameter
        or a synonym of it, e.g. 'euclidean' if the `metric` parameter set to
        'minkowski' and `p` parameter set to 2.

    effective_metric_params_ : dict
        Additional keyword arguments for the metric function. For most metrics
        will be same with `metric_params` parameter, but may also contain the
        `p` parameter value if the `effective_metric_` attribute is set to
        'minkowski'.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

        .. versionadded:: 1.0

    n_samples_fit_ : int
        Number of samples in the fitted data.

    outlier_label_ : int or array-like of shape (n_class,)
        Label which is given for outlier samples (samples with no neighbors
        on given radius).

    outputs_2d_ : bool
        False when `y`'s shape is (n_samples, ) or (n_samples, 1) during fit
        otherwise True.

    See Also
    --------
    KNeighborsClassifier : Classifier implementing the k-nearest neighbors
        vote.
    RadiusNeighborsRegressor : Regression based on neighbors within a
        fixed radius.
    KNeighborsRegressor : Regression based on k-nearest neighbors.
    NearestNeighbors : Unsupervised learner for implementing neighbor
        searches.

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    https://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import RadiusNeighborsClassifier
    >>> neigh = RadiusNeighborsClassifier(radius=1.0)
    >>> neigh.fit(X, y)
    RadiusNeighborsClassifier(...)
    >>> print(neigh.predict([[1.5]]))
    [0]
    >>> print(neigh.predict_proba([[1.0]]))
    [[0.66666667 0.33333333]]
    """

    _parameter_constraints: dict = {
        **NeighborsBase._parameter_constraints,
        "weights": [StrOptions({"uniform", "distance"}), callable, None],
        "outlier_label": [Integral, str, "array-like", None],
    }
    _parameter_constraints.pop("n_neighbors")

    def __init__(
        self,
        radius=1.0,
        *,
        weights="uniform",
        algorithm="auto",
        leaf_size=30,
        p=2,
        metric="minkowski",
        outlier_label=None,
        metric_params=None,
        n_jobs=None,
    ):
        super().__init__(
            radius=radius,
            algorithm=algorithm,
            leaf_size=leaf_size,
            metric=metric,
            p=p,
            metric_params=metric_params,
            n_jobs=n_jobs,
        )
        self.weights = weights
        self.outlier_label = outlier_label

    def fit(self, X, y):
        """Fit the radius neighbors classifier from the training dataset.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features) or \
                (n_samples, n_samples) if metric='precomputed'
            Training data.

        y : {array-like, sparse matrix} of shape (n_samples,) or \
                (n_samples, n_outputs)
            Target values.

        Returns
        -------
        self : RadiusNeighborsClassifier
            The fitted radius neighbors classifier.
        """
        self._validate_params()
        self._fit(X, y)

        classes_ = self.classes_
        _y = self._y
        if not self.outputs_2d_:
            _y = self._y.reshape((-1, 1))
            classes_ = [self.classes_]

        if self.outlier_label is None:
            outlier_label_ = None

        elif self.outlier_label == "most_frequent":
            outlier_label_ = []
            # iterate over multi-output, get the most frequent label for each
            # output.
            for k, classes_k in enumerate(classes_):
                label_count = np.bincount(_y[:, k])
                outlier_label_.append(classes_k[label_count.argmax()])

        else:
            if _is_arraylike(self.outlier_label) and not isinstance(
                self.outlier_label, str
            ):
                if len(self.outlier_label) != len(classes_):
                    raise ValueError(
                        "The length of outlier_label: {} is "
                        "inconsistent with the output "
                        "length: {}".format(self.outlier_label, len(classes_))
                    )
                outlier_label_ = self.outlier_label
            else:
                outlier_label_ = [self.outlier_label] * len(classes_)

            for classes, label in zip(classes_, outlier_label_):
                if _is_arraylike(label) and not isinstance(label, str):
                    # ensure the outlier label for each output is a scalar.
                    raise TypeError(
                        "The outlier_label of classes {} is "
                        "supposed to be a scalar, got "
                        "{}.".format(classes, label)
                    )
                if np.append(classes, label).dtype != classes.dtype:
                    # ensure the dtype of outlier label is consistent with y.
                    raise TypeError(
                        "The dtype of outlier_label {} is "
                        "inconsistent with classes {} in "
                        "y.".format(label, classes)
                    )

        self.outlier_label_ = outlier_label_

        return self

    def predict(self, X):
        """Predict the class labels for the provided data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_queries, n_features), \
                or (n_queries, n_indexed) if metric == 'precomputed'
            Test samples.

        Returns
        -------
        y : ndarray of shape (n_queries,) or (n_queries, n_outputs)
            Class labels for each data sample.
        """

        probs = self.predict_proba(X)
        classes_ = self.classes_

        if not self.outputs_2d_:
            probs = [probs]
            classes_ = [self.classes_]

        n_outputs = len(classes_)
        n_queries = probs[0].shape[0]
        y_pred = np.empty((n_queries, n_outputs), dtype=classes_[0].dtype)

        for k, prob in enumerate(probs):
            # iterate over multi-output, assign labels based on probabilities
            # of each output.
            max_prob_index = prob.argmax(axis=1)
            y_pred[:, k] = classes_[k].take(max_prob_index)

            outlier_zero_probs = (prob == 0).all(axis=1)
            if outlier_zero_probs.any():
                zero_prob_index = np.flatnonzero(outlier_zero_probs)
                y_pred[zero_prob_index, k] = self.outlier_label_[k]

        if not self.outputs_2d_:
            y_pred = y_pred.ravel()

        return y_pred

    def predict_proba(self, X):
        """Return probability estimates for the test data X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_queries, n_features), \
                or (n_queries, n_indexed) if metric == 'precomputed'
            Test samples.

        Returns
        -------
        p : ndarray of shape (n_queries, n_classes), or a list of \
                n_outputs of such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are ordered
            by lexicographic order.
        """

        n_queries = _num_samples(X)

        neigh_dist, neigh_ind = self.radius_neighbors(X)
        outlier_mask = np.zeros(n_queries, dtype=bool)
        outlier_mask[:] = [len(nind) == 0 for nind in neigh_ind]
        outliers = np.flatnonzero(outlier_mask)
        inliers = np.flatnonzero(~outlier_mask)

        classes_ = self.classes_
        _y = self._y
        if not self.outputs_2d_:
            _y = self._y.reshape((-1, 1))
            classes_ = [self.classes_]

        if self.outlier_label_ is None and outliers.size > 0:
            raise ValueError(
                "No neighbors found for test samples %r, "
                "you can try using larger radius, "
                "giving a label for outliers, "
                "or considering removing them from your dataset." % outliers
            )

        weights = _get_weights(neigh_dist, self.weights)
        if weights is not None:
            weights = weights[inliers]

        probabilities = []
        # iterate over multi-output, measure probabilities of the k-th output.
        for k, classes_k in enumerate(classes_):
            pred_labels = np.zeros(len(neigh_ind), dtype=object)
            pred_labels[:] = [_y[ind, k] for ind in neigh_ind]

            proba_k = np.zeros((n_queries, classes_k.size))
            proba_inl = np.zeros((len(inliers), classes_k.size))

            # samples have different size of neighbors within the same radius
            if weights is None:
                for i, idx in enumerate(pred_labels[inliers]):
                    proba_inl[i, :] = np.bincount(idx, minlength=classes_k.size)
            else:
                for i, idx in enumerate(pred_labels[inliers]):
                    proba_inl[i, :] = np.bincount(
                        idx, weights[i], minlength=classes_k.size
                    )
            proba_k[inliers, :] = proba_inl

            if outliers.size > 0:
                _outlier_label = self.outlier_label_[k]
                label_index = np.flatnonzero(classes_k == _outlier_label)
                if label_index.size == 1:
                    proba_k[outliers, label_index[0]] = 1.0
                else:
                    warnings.warn(
                        "Outlier label {} is not in training "
                        "classes. All class probabilities of "
                        "outliers will be assigned with 0."
                        "".format(self.outlier_label_[k])
                    )

            # normalize 'votes' into real [0,1] probabilities
            normalizer = proba_k.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            proba_k /= normalizer

            probabilities.append(proba_k)

        if not self.outputs_2d_:
            probabilities = probabilities[0]

        return probabilities

    def _more_tags(self):
        return {"multilabel": True}
