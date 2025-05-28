"""
Classical multi-dimensional scaling (classical MDS).
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from numbers import Integral

import numpy as np
from scipy import linalg

from ..base import BaseEstimator, _fit_context
from ..metrics import pairwise_distances
from ..utils import check_symmetric
from ..utils._param_validation import Interval
from ..utils.extmath import svd_flip
from ..utils.validation import validate_data


class ClassicalMDS(BaseEstimator):
    """Classical multidimensional scaling.

    Read more in the :ref:`User Guide <multidimensional_scaling>`.

    Parameters
    ----------
    n_components : int, default=2
        Number of embedding dimensions.

    dissimilarity : str or callable, default='euclidean'
        Metric to use for dissimilarity computation. Default is "euclidean".
        See the documentation of `scipy.spatial.distance
        <https://docs.scipy.org/doc/scipy/reference/spatial.distance.html>`_ and
        the metrics listed in
        :class:`~sklearn.metrics.pairwise.distance_metrics` for valid metric
        values.

        If metric is "precomputed", X is assumed to be a distance matrix and
        must be square during fit.

        If metric is a callable function, it takes two arrays representing 1D
        vectors as inputs and must return one value indicating the distance
        between those vectors. This works for Scipy's metrics, but is less
        efficient than passing the metric name as a string.

    metric_params : dict, default=None
        Additional keyword arguments for the metric function.

    Attributes
    ----------
    embedding_ : ndarray of shape (n_samples, n_components)
        Stores the position of the dataset in the embedding space.

    dissimilarity_matrix_ : ndarray of shape (n_samples, n_samples)
        Pairwise dissimilarities between the points.

    eigenvalues_ : ndarray of shape (n_components,)
        Eigenvalues of the double-centered dissimilarity matrix, corresponding
        to each of the selected components. They are equal to the squared 2-norms
        of the `n_components` variables in the embedding space.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when `X`
        has feature names that are all strings.

    See Also
    --------
    sklearn.decomposition.PCA : Principal component analysis.
    MDS : Metric and non-metric MDS.

    References
    ----------
    .. [1] "Modern Multidimensional Scaling - Theory and Applications" Borg, I.;
       Groenen P. Springer Series in Statistics (1997)

    Examples
    --------
    >>> from sklearn.datasets import load_digits
    >>> from sklearn.manifold import ClassicalMDS
    >>> X, _ = load_digits(return_X_y=True)
    >>> X.shape
    (1797, 64)
    >>> cmds = ClassicalMDS(n_components=2)
    >>> X_emb = cmds.fit_transform(X[:100])
    >>> X_emb.shape
    (100, 2)
    """

    _parameter_constraints: dict = {
        "n_components": [Interval(Integral, 1, None, closed="left")],
        "dissimilarity": [str, callable],
        "metric_params": [dict, None],
    }

    def __init__(
        self,
        n_components=2,
        *,
        dissimilarity="euclidean",
        metric_params=None,
    ):
        self.n_components = n_components
        self.dissimilarity = dissimilarity
        self.metric_params = metric_params

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.input_tags.pairwise = self.dissimilarity == "precomputed"
        return tags

    def fit(self, X, y=None):
        """
        Compute the embedding positions.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or \
                (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        self.fit_transform(X)
        return self

    @_fit_context(prefer_skip_nested_validation=True)
    def fit_transform(self, X, y=None):
        """
        Compute and return the embedding positions.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or \
                (n_samples, n_samples)
            Input data. If ``dissimilarity=='precomputed'``, the input should
            be the dissimilarity matrix.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            The embedding coordinates.
        """

        X = validate_data(self, X)

        if self.dissimilarity == "precomputed":
            self.dissimilarity_matrix_ = X
            self.dissimilarity_matrix_ = check_symmetric(
                self.dissimilarity_matrix_, raise_exception=True
            )
        else:
            self.dissimilarity_matrix_ = pairwise_distances(
                X, metric=self.dissimilarity, **self.metric_params
            )

        # Double centering
        B = self.dissimilarity_matrix_**2
        B = B.astype(np.float64)
        B -= np.mean(B, axis=0)
        B -= np.mean(B, axis=1, keepdims=True)
        B *= -0.5

        # Eigendecomposition
        w, U = linalg.eigh(B)

        # Reversing the order of the eigenvalues/eigenvectors to put
        # the eigenvalues in decreasing order
        w = w[::-1][: self.n_components]
        U = U[:, ::-1][:, : self.n_components]

        # Set the signs of eigenvectors to enforce deterministic output
        U, _ = svd_flip(U, None)

        self.embedding_ = np.sqrt(w) * U
        self.eigenvalues_ = w

        return self.embedding_
