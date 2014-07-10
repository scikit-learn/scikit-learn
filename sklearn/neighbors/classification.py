"""Nearest Neighbor Classification"""

# Authors: Jake Vanderplas <vanderplas@astro.washington.edu>
#          Fabian Pedregosa <fabian.pedregosa@inria.fr>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
#          Sparseness support by Lars Buitinck <L.J.Buitinck@uva.nl>
#          Multi-output support by Arnaud Joly <a.joly@ulg.ac.be>
#
# License: BSD 3 clause (C) INRIA, University of Amsterdam

import numpy as np
import scipy.sparse as sp
from scipy import stats
from ..utils.extmath import weighted_mode
from ..utils.sparsefuncs_fast import csr_row_mode

from .base import \
    _check_weights, _get_weights, \
    NeighborsBase, KNeighborsMixin,\
    RadiusNeighborsMixin, SupervisedIntegerMixin
from ..base import ClassifierMixin
from ..utils import check_array


class KNeighborsClassifier(NeighborsBase, KNeighborsMixin,
                           SupervisedIntegerMixin, ClassifierMixin):
    """Classifier implementing the k-nearest neighbors vote.

    Parameters
    ----------
    n_neighbors : int, optional (default = 5)
        Number of neighbors to use by default for :meth:`k_neighbors` queries.

    weights : str or callable
        weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

        Uniform weights are used by default.

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
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    metric : string or DistanceMetric object (default='minkowski')
        the distance metric to use for the tree.  The default metric is
        minkowski, and with p=2 is equivalent to the standard Euclidean
        metric. See the documentation of the DistanceMetric class for a
        list of available metrics.

    p : integer, optional (default = 2)
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    **kwargs :
        additional keyword arguments are passed to the distance function as
        additional arguments.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import KNeighborsClassifier
    >>> neigh = KNeighborsClassifier(n_neighbors=3)
    >>> neigh.fit(X, y) # doctest: +ELLIPSIS
    KNeighborsClassifier(...)
    >>> print(neigh.predict([[1.1]]))
    [0]
    >>> print(neigh.predict_proba([[0.9]]))
    [[ 0.66666667  0.33333333]]

    See also
    --------
    RadiusNeighborsClassifier
    KNeighborsRegressor
    RadiusNeighborsRegressor
    NearestNeighbors

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    .. warning::

       Regarding the Nearest Neighbors algorithms, if it is found that two
       neighbors, neighbor `k+1` and `k`, have identical distances but
       but different labels, the results will depend on the ordering of the
       training data.

    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, n_neighbors=5,
                 weights='uniform', algorithm='auto', leaf_size=30,
                 p=2, metric='minkowski', **kwargs):
        self._init_params(n_neighbors=n_neighbors,
                          algorithm=algorithm,
                          leaf_size=leaf_size, metric=metric, p=p, **kwargs)
        self.weights = _check_weights(weights)

    def predict(self, X):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            A 2-D array representing the test points.

        Returns
        -------
        y : array of shape [n_samples] or [n_samples, n_outputs]
            Class labels for each data sample.
        """
        X = check_array(X, accept_sparse='csr')

        # XXX Current
        # X = atleast2d_or_csr(X)

        # neigh_dist, neigh_ind = self.kneighbors(X)

        # classes_ = self.classes_
        # _y = self._y
        # if not self.outputs_2d_:
        #     _y = self._y.reshape((-1, 1))
        #     classes_ = [self.classes_]

        # n_outputs = len(classes_)
        # n_samples = X.shape[0]
        # weights = _get_weights(neigh_dist, self.weights)

        # y_pred = np.empty((n_samples, n_outputs), dtype=classes_[0].dtype)
        # for k, classes_k in enumerate(classes_):
        #     if weights is None:
        #         mode, _ = stats.mode(_y[neigh_ind, k], axis=1)
        #     else:
        #         mode, _ = weighted_mode(_y[neigh_ind, k], weights, axis=1)

        #     mode = np.asarray(mode.ravel(), dtype=np.intp)
        #     y_pred[:, k] = classes_k.take(mode)

        # if not self.outputs_2d_:
        #     y_pred = y_pred.ravel()

        # return y_pred

        # XXX sparse
        X = atleast2d_or_csr(X)

        neigh_dist, neigh_ind = self.kneighbors(X)

        classes_ = self.classes_
        _y = self._y
        if not self.outputs_2d_:
            # XXX sparsiify this
            _y = self._y.reshape((-1, 1))
            classes_ = [self.classes_]

        n_outputs = len(classes_)
        n_samples = X.shape[0]
        weights = _get_weights(neigh_dist, self.weights)

        if not self.sparse_target_input_:
            _y = sp.csc_matrix(_y) # XXX Remove when sparse matrix comes from fit in one continous pass

        data = np.array([], dtype=classes_[0].dtype)
        indices = np.array([])
        indptr = np.array([0])
        for k, classes_k in enumerate(classes_):
            neigh_lbls_k = _y[neigh_ind,k]
            if weights is None:
                mode = csr_row_mode(sp.csr_matrix(neigh_lbls_k.T))
                mode = sp.csc_matrix(mode ,dtype=np.intp)
            else:
                # XXX make a sparse weighted_mode function, do not densify neigh_lbls_k
                mode, _ = weighted_mode(neigh_lbls_k.toarray().T, weights, axis=1)
                mode = sp.csc_matrix(mode, dtype=np.intp)

            # XXX Changes appends to concats
            data = np.append(data, mode.data)
            indices = np.append(indices, mode.indices)
            indptr = np.append(indptr, len(indices))

        y_pred = sp.csc_matrix((data, indices, indptr), (n_samples, n_outputs), dtype=np.intp)

        if not self.sparse_target_input_:
            y_pred = y_pred.toarray() # XXX remove this
            y_pred = classes_k[y_pred] # XXX remove this for sparse cases, sparse matricies cant have str data

        if not self.outputs_2d_:
            y_pred = y_pred.ravel()

        return y_pred

    def predict_proba(self, X):
        """Return probability estimates for the test data X.

        Parameters
        ----------
        X : array, shape = (n_samples, n_features)
            A 2-D array representing the test points.

        Returns
        -------
        p : array of shape = [n_samples, n_classes], or a list of n_outputs
            of such arrays if n_outputs > 1.
            The class probabilities of the input samples. Classes are ordered
            by lexicographic order.
        """
        X = check_array(X, accept_sparse='csr')

        neigh_dist, neigh_ind = self.kneighbors(X)

        classes_ = self.classes_
        _y = self._y

        # XXX temporary patch
        if self.sparse_target_input_:
            _y = _y.toarray()

        if not self.outputs_2d_:
            _y = self._y.reshape((-1, 1))
            classes_ = [self.classes_]

        n_samples = X.shape[0]

        weights = _get_weights(neigh_dist, self.weights)
        if weights is None:
            weights = np.ones_like(neigh_ind)
        else:
            # Some weights may be infinite (zero distance), which can cause
            # downstream NaN values when used for normalization.
            weights[np.isinf(weights)] = np.finfo('f').max

        all_rows = np.arange(X.shape[0])
        probabilities = []
        for k, classes_k in enumerate(classes_):
            pred_labels = _y[:, k][neigh_ind]
            proba_k = np.zeros((n_samples, classes_k.size))

            # a simple ':' index doesn't work right
            for i, idx in enumerate(pred_labels.T):  # loop is O(n_neighbors)
                proba_k[all_rows, idx] += weights[:, i]

            # normalize 'votes' into real [0,1] probabilities
            normalizer = proba_k.sum(axis=1)[:, np.newaxis]
            normalizer[normalizer == 0.0] = 1.0
            proba_k /= normalizer

            probabilities.append(proba_k)

        if not self.outputs_2d_:
            probabilities = probabilities[0]

        return probabilities


class RadiusNeighborsClassifier(NeighborsBase, RadiusNeighborsMixin,
                                SupervisedIntegerMixin, ClassifierMixin):
    """Classifier implementing a vote among neighbors within a given radius

    Parameters
    ----------
    radius : float, optional (default = 1.0)
        Range of parameter space to use by default for :meth`radius_neighbors`
        queries.

    weights : str or callable
        weight function used in prediction.  Possible values:

        - 'uniform' : uniform weights.  All points in each neighborhood
          are weighted equally.
        - 'distance' : weight points by the inverse of their distance.
          in this case, closer neighbors of a query point will have a
          greater influence than neighbors which are further away.
        - [callable] : a user-defined function which accepts an
          array of distances, and returns an array of the same shape
          containing the weights.

        Uniform weights are used by default.

    algorithm : {'auto', 'ball_tree', 'kd_tree', 'brute'}, optional
        Algorithm used to compute the nearest neighbors:

        - 'ball_tree' will use :class:`BallTree`
        - 'kd_tree' will use :class:`KDtree`
        - 'brute' will use a brute-force search.
        - 'auto' will attempt to decide the most appropriate algorithm
          based on the values passed to :meth:`fit` method.

        Note: fitting on sparse input will override the setting of
        this parameter, using brute force.

    leaf_size : int, optional (default = 30)
        Leaf size passed to BallTree or KDTree.  This can affect the
        speed of the construction and query, as well as the memory
        required to store the tree.  The optimal value depends on the
        nature of the problem.

    metric : string or DistanceMetric object (default='minkowski')
        the distance metric to use for the tree.  The default metric is
        minkowski, and with p=2 is equivalent to the standard Euclidean
        metric. See the documentation of the DistanceMetric class for a
        list of available metrics.

    p : integer, optional (default = 2)
        Power parameter for the Minkowski metric. When p = 1, this is
        equivalent to using manhattan_distance (l1), and euclidean_distance
        (l2) for p = 2. For arbitrary p, minkowski_distance (l_p) is used.

    outlier_label : int, optional (default = None)
        Label, which is given for outlier samples (samples with no
        neighbors on given radius).
        If set to None, ValueError is raised, when outlier is detected.

    **kwargs :
        additional keyword arguments are passed to the distance function as
        additional arguments.

    Examples
    --------
    >>> X = [[0], [1], [2], [3]]
    >>> y = [0, 0, 1, 1]
    >>> from sklearn.neighbors import RadiusNeighborsClassifier
    >>> neigh = RadiusNeighborsClassifier(radius=1.0)
    >>> neigh.fit(X, y) # doctest: +ELLIPSIS
    RadiusNeighborsClassifier(...)
    >>> print(neigh.predict([[1.5]]))
    [0]

    See also
    --------
    KNeighborsClassifier
    RadiusNeighborsRegressor
    KNeighborsRegressor
    NearestNeighbors

    Notes
    -----
    See :ref:`Nearest Neighbors <neighbors>` in the online documentation
    for a discussion of the choice of ``algorithm`` and ``leaf_size``.

    http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm
    """

    def __init__(self, radius=1.0, weights='uniform',
                 algorithm='auto', leaf_size=30, p=2, metric='minkowski',
                 outlier_label=None, **kwargs):
        self._init_params(radius=radius,
                          algorithm=algorithm,
                          leaf_size=leaf_size,
                          metric=metric, p=p, **kwargs)
        self.weights = _check_weights(weights)
        self.outlier_label = outlier_label

    def predict(self, X):
        """Predict the class labels for the provided data

        Parameters
        ----------
        X : array of shape [n_samples, n_features]
            A 2-D array representing the test points.

        Returns
        -------
        y : array of shape [n_samples] or [n_samples, n_outputs]
            Class labels for each data sample.

        """
        X = check_array(X, accept_sparse='csr')
        n_samples = X.shape[0]

        neigh_dist, neigh_ind = self.radius_neighbors(X)
        inliers = [i for i, nind in enumerate(neigh_ind) if len(nind) != 0]
        outliers = [i for i, nind in enumerate(neigh_ind) if len(nind) == 0]

        classes_ = self.classes_
        _y = self._y
        if not self.outputs_2d_:
            _y = self._y.reshape((-1, 1))
            classes_ = [self.classes_]
        n_outputs = len(classes_)

        if self.outlier_label is not None:
            neigh_dist[outliers] = 1e-6
        elif outliers:
            raise ValueError('No neighbors found for test samples %r, '
                             'you can try using larger radius, '
                             'give a label for outliers, '
                             'or consider removing them from your dataset.'
                             % outliers)

        weights = _get_weights(neigh_dist, self.weights)

        y_pred = np.empty((n_samples, n_outputs), dtype=classes_[0].dtype)
        for k, classes_k in enumerate(classes_):
            pred_labels = np.array([_y[ind, k] for ind in neigh_ind],
                                   dtype=object)
            if weights is None:
                mode = np.array([stats.mode(pl)[0]
                                 for pl in pred_labels[inliers]], dtype=np.int)
            else:
                mode = np.array([weighted_mode(pl, w)[0]
                                 for (pl, w)
                                 in zip(pred_labels[inliers], weights)],
                                dtype=np.int)

            mode = mode.ravel()

            y_pred[inliers, k] = classes_k.take(mode)

        if outliers:
            y_pred[outliers, :] = self.outlier_label

        if not self.outputs_2d_:
            y_pred = y_pred.ravel()

        return y_pred
