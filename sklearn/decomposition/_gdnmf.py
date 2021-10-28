"""
Graph regularized discriminative non-negative matrix factorization.
"""

# Author: Jason Zhang <jason.zhang@ucb.com>
# License: BSD 3 clause

from typing import Optional, Tuple
from scipy.sparse import csgraph
from sklearn.neighbors import NearestNeighbors
import numbers
import numpy as np
import time
import warnings

from ._nmf import trace_dot
from ..base import BaseEstimator, TransformerMixin, clone
from ..exceptions import ConvergenceWarning
from ..utils import check_random_state
from ..utils.extmath import safe_sparse_dot
from ..utils.validation import check_is_fitted, check_non_negative
from .._config import config_context

EPSILON = np.finfo(np.float32).eps


def _gdnmf_objective(X, C, B, S, W, H, A, graph_coeff, label_coeff) -> float:
    """
    Compute the objective function (loss) for GDNMF.
    """
    L = B - C
    reconstruction_error = np.trace(safe_sparse_dot(X - W @ H, (X - W @ H).T))
    graph_regularization = graph_coeff * np.trace(W.T @ L @ W)
    label_regularization = label_coeff * trace_dot((S - A @ W.T).T, S - A @ W.T)
    return reconstruction_error + graph_regularization + label_regularization


def _init_w_h_a(X, n_components: int, n_classes: int, random_state=None):
    check_non_negative(X, "NMF initialization")
    n_samples, n_features = X.shape

    avg = np.sqrt(X.mean() / n_components)
    rng = check_random_state(random_state)
    H = avg * rng.randn(n_components, n_features).astype(X.dtype, copy=False)
    W = avg * rng.randn(n_samples, n_components).astype(X.dtype, copy=False)
    A = avg * rng.randn(n_classes, n_components).astype(X.dtype, copy=False)
    np.abs(H, out=H)
    np.abs(W, out=W)
    np.abs(A, out=A)
    return W, H, A


def _multiplicative_update_w(X, C, B, S, W, H, A, graph_coeff, label_coeff):
    nominator = label_coeff * A.T @ S + safe_sparse_dot(H, X.T) + graph_coeff * W.T @ C
    nominator[nominator <= 0] = EPSILON
    denominator = H @ H.T @ W.T + label_coeff * A.T @ A @ W.T + graph_coeff * W.T @ B
    denominator[denominator <= 0] = 1.0
    return nominator / denominator


def _multiplicative_update_h(X, W, H):
    nominator = safe_sparse_dot(X.T, W)
    nominator[nominator <= 0] = EPSILON
    denominator = H.T @ W.T @ W
    denominator[denominator <= 0] = 1.0
    return nominator / denominator


def _multiplicative_update_a(S, W, A):
    nominator = S @ W
    nominator[nominator <= 0] = EPSILON
    denominator = A @ W.T @ W
    denominator[denominator <= 0] = 1.0
    return nominator / denominator


def _fit_multiplicative_update(
    X,
    C,
    B,
    S,
    W,
    H,
    A,
    graph_coeff: float = 0.0,
    label_coeff: float = 0.0,
    max_iter: int = 200,
    tol: float = 1e-4,
    verbose: int = 0,
) -> Tuple[np.array, np.array, int]:
    """
    Compute Graph Regularized Discriminative Non-negative Matrix
    Factorization with Multiplicative Update.

    The objective function is _beta_divergence(X, WH) and is minimized with an
    alternating minimization of W and H. Each minimization is done with a
    Multiplicative Update.

    Parameters
    ----------
    X : array-like of shape (n_samples, n_features)
        Constant input matrix.

    C : array-like of shape (n_samples, n_samples)
        Samples neighboring relationship graph.

    B : array-like of shape (n_samples, n_samples)
        Row sum of the weight matrix C, used for graph laplacian calculation.

    S : array-like of shape (n_classes, n_samples)
        Label information matrix.

    W : array-like of shape (n_samples, n_components)
        Initial guess for the solution.

    H : array-like of shape (n_components, n_features)
        Initial guess for the solution.

    A : array-like of shape (n_classes, n_components)
        Initial guess for the solution.

    graph_coeff : float, default=0.0
        Coefficient for the graph Laplacian regularization term.

    label_coeff : float, default=0.0
        Coefficient for the label information regularization term.

    max_iter : int, default=200
        Number of iterations.

    tol : float, default=1e-4
        Tolerance of the stopping condition.

    verbose : int, default=0
        The verbosity level.

    Returns
    -------
    W : ndarray of shape (n_samples, n_components)
        Left factor matrix in :math:`X = WH`.

    H : ndarray of shape (n_components, n_features)
        Right factor matrix in :math:`X = WH`.

    n_iter : int
        The number of iterations done by the algorithm.

    References
    ----------
    Long, X., Lu, H., Peng, Y., & Li, W. (2014). Graph regularized
    discriminative non-negative matrix factorization for face recognition.
    Multimedia tools and applications, 72(3), 2679-2699.
    """
    start_time = time.time()
    previous_error = error = error_at_init = _gdnmf_objective(
        X=X,
        C=C,
        B=B,
        S=S,
        W=W,
        H=H,
        A=A,
        graph_coeff=graph_coeff,
        label_coeff=label_coeff,
    )

    for n_iter in range(1, max_iter + 1):
        # update W
        # H_sum, HHt and XHt are saved and reused if not update_H
        delta_W = _multiplicative_update_w(
            X=X,
            C=C,
            B=B,
            S=S,
            W=W,
            H=H,
            A=A,
            graph_coeff=graph_coeff,
            label_coeff=label_coeff,
        ).T
        W *= delta_W

        # update H
        delta_H = _multiplicative_update_h(X=X, W=W, H=H).T
        H *= delta_H

        # update A
        delta_A = _multiplicative_update_a(S=S, W=W, A=A)
        A *= delta_A

        # test convergence criterion every 10 iterations
        if tol > 0 and n_iter % 10 == 0:
            error = _gdnmf_objective(
                X=X,
                C=C,
                B=B,
                S=S,
                W=W,
                H=H,
                A=A,
                graph_coeff=graph_coeff,
                label_coeff=label_coeff,
            )

            if verbose:
                iter_time = time.time()
                print(
                    f"Epoch {n_iter:02d} reached after "
                    f"{iter_time - start_time:.3f} seconds, error: {error}"
                )

            if (previous_error - error) / error_at_init < tol:
                break
            previous_error = error

    # do not print if we have already printed in the convergence test
    if verbose and (tol == 0 or n_iter % 10 != 0):
        end_time = time.time()
        print(f"Epoch {n_iter:02d} reached after {end_time - start_time:.3f} seconds.")

    return W, H, n_iter


class GDNMF(TransformerMixin, BaseEstimator):
    """
    Graph regularized discriminative non-negative matrix factorization.

    Parameters
    ----------
    n_components : int, default=None
        Number of components, if n_components is not set all features
        are kept.

    n_classes : int, default=None
        Number of unique classes, if n_classes is not set we will
        derive from y.

    tol : float, default=1e-4
        Tolerance of the stopping condition.

    max_iter : int, default=200
        Maximum number of iterations before timing out.

    graph_coeff : float, default=0.0
        Coefficient for the graph Laplacian regularization term.

    label_coeff : float, default=0.0
        Coefficient for the label information regularization term.

    random_state : int, RandomState instance or None, default=None
        Pass an int for reproducible results across multiple function calls.
        See :term:`Glossary <random_state>`.

    knn : NearestNeighbors, default=None
        Pre-configured nearest-neighbors calculator.

    verbose : int, default=0
        The verbosity level.

    Attributes
    ----------
    components_ : ndarray of shape (n_components, n_features)
        Factorization matrix, sometimes called 'dictionary'.

    n_components_ : int
        The number of components. It is same as the `n_components` parameter
        if it was given. Otherwise, it will be same as the number of
        features.

    labels_ : array-like of shape (n_samples)
        Vector of the class labels for each sample.

    n_classes_ : int
        The number of unique classes in the labels. It is same as the
        `n_classes` parameter if it was given. Otherwise, it will be
        same as the number of classes.

    reconstruction_err_ : float
        Frobenius norm of the matrix difference, between
        the training data ``X`` and the reconstructed data ``WH`` from
        the fitted model.

    n_iter_ : int
        Actual number of iterations.

    n_features_in_ : int
        Number of features seen during :term:`fit`.

    See Also
    --------
    NMF : Non-negative matrix factorization without labels.
    DictionaryLearning : Find a dictionary that sparsely encodes data.
    MiniBatchSparsePCA : Mini-batch Sparse Principal Components Analysis.
    PCA : Principal component analysis.
    SparseCoder : Find a sparse representation of data from a fixed,
        precomputed dictionary.
    SparsePCA : Sparse Principal Components Analysis.
    TruncatedSVD : Dimensionality reduction using truncated SVD.

    References
    ----------
    Long, X., Lu, H., Peng, Y., & Li, W. (2014). Graph regularized
    discriminative non-negative matrix factorization for face recognition.
    Multimedia tools and applications, 72(3), 2679-2699.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[1, 1], [2, 1], [3, 1.2], [4, 1], [5, 0.8], [6, 1]])
    >>> y = np.array([0, 0, 1, 0, 1, 1])
    >>> from sklearn.decomposition import GDNMF
    >>> model = GDNMF(n_components=2, n_classes=2,
    ... graph_coeff=0.1, label_coeff=0.1, random_state=0)
    >>> W = model.fit_transform(X, y)
    >>> H = model.components_
    """

    def __init__(
        self,
        n_components: Optional[int] = None,
        n_classes: Optional[int] = None,
        tol: Optional[float] = 1e-4,
        max_iter: Optional[int] = 200,
        graph_coeff: Optional[float] = 0.0,
        label_coeff: Optional[float] = 0.0,
        random_state=None,
        knn: Optional[NearestNeighbors] = None,
        verbose=0,
    ):
        self.n_components = n_components
        self.n_classes = n_classes
        self.tol = tol
        self.max_iter = max_iter
        self.graph_coeff = graph_coeff
        self.label_coeff = label_coeff
        self.random_state = random_state
        self.knn = knn
        self.verbose = verbose

    def fit_transform(self, X, y):
        """
        Learn a GDNMF model for the data X, y and returns the transformed data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : array-like of shape (n_samples)
            Vector of the class labels for each sample.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        X, y = self._validate_data(
            X=X, y=y, accept_sparse=("csr", "csc"), dtype=[np.float64, np.float32]
        )

        with config_context(assume_finite=True):
            W, H, n_iter = self._fit_transform(X, y)

        self.n_components_ = H.shape[0]
        self.components_ = H
        self.labels_ = y
        self.n_iter_ = n_iter

        return W

    def _fit_transform(self, X, y):
        """
        Learn a GDNMF model for the data X and returns the transformed data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Data matrix to be decomposed

        y : array-like of shape (n_samples)
            Vector of the class labels for each sample.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.

        H : ndarray of shape (n_components, n_features)
            Factorization matrix, sometimes called 'dictionary'.

        n_iter_ : int
            Actual number of iterations.
        """
        check_non_negative(X, "NMF (input X)")

        # check parameters
        self._check_params(X, y)

        # calculate graph laplacian
        self.knn_.fit(X)
        C = self.knn_.kneighbors_graph(X, mode="connectivity").toarray()
        L = csgraph.laplacian(C)
        B = L + C

        # calculate label information matrix
        S = np.array(
            [
                np.array([1 if label == cls else 0 for label in y])
                for cls in range(self.n_classes_)
            ]
        )

        # initialize or check W and H
        W, H, A = _init_w_h_a(X, self.n_components_, self.n_classes_, self.random_state)

        W, H, n_iter = _fit_multiplicative_update(
            X=X,
            C=C,
            B=B,
            S=S,
            W=W,
            H=H,
            A=A,
            graph_coeff=self.graph_coeff,
            label_coeff=self.label_coeff,
            max_iter=self.max_iter,
            tol=self.tol,
            verbose=self.verbose,
        )
        self.reconstruction_err_ = _gdnmf_objective(
            X=X, C=C, B=B, S=S, W=W, H=H, A=A, graph_coeff=0, label_coeff=0
        )

        if n_iter == self.max_iter and self.tol > 0:
            warnings.warn(
                f"Maximum number of iterations {self.max_iter} reached. Increase "
                "it to improve convergence.",
                ConvergenceWarning,
            )

        return W, H, n_iter

    def fit(self, X, y):
        """
        Learn a NMF model for the data X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : array-like of shape (n_samples)
            Vector of the class labels for each sample.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        self.fit_transform(X, y)
        return self

    def transform(self, X):
        """
        Transform the data X according to the fitted GDNMF model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        W : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        check_is_fitted(self)
        X = self._validate_data(
            X,
            accept_sparse=("csr", "csc"),
            dtype=[np.float64, np.float32],
            reset=False,
        )
        y = self._validate_data(y=self.labels_[: X.shape[0]])
        if X.shape[1] != self.components_.shape[1]:
            raise ValueError(
                "Number of features to transform must match that fitted"
                f"(n_to_transform={X.shape[1]}, "
                f"n_fitted={self.components_.shape[1]})"
            )

        with config_context(assume_finite=True):
            W, *_ = self._fit_transform(X, y)

        return W

    def inverse_transform(self, W):
        """
        Transform data back to its original space.

        Parameters
        ----------
        W : {ndarray, sparse matrix} of shape (n_samples, n_components)
            Transformed data matrix.

        Returns
        -------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            Returns a data matrix of the original shape.
        """
        check_is_fitted(self)
        return np.dot(W, self.components_)

    def _check_params(self, X, y):
        # n_components
        self.n_components_ = self.n_components
        if self.n_components_ is None:
            self.n_components_ = X.shape[1]
        if (
            not isinstance(self.n_components_, numbers.Integral)
            or self.n_components_ <= 0
        ):
            raise ValueError(
                "Number of components must be a positive integer; got "
                f"(n_components={self.n_components_!r})"
            )

        # n_classes
        self.n_classes_ = self.n_classes
        if self.n_classes_ is None:
            self.n_classes_ = len(set(y))
        if not isinstance(self.n_classes_, numbers.Integral) or self.n_classes_ <= 0:
            raise ValueError(
                "Number of classes must be a positive integer; got "
                f"(n_classes={self.n_classes_!r})"
            )

        # max_iter
        if not isinstance(self.max_iter, numbers.Integral) or self.max_iter < 0:
            raise ValueError(
                "Maximum number of iterations must be a positive "
                f"integer; got (max_iter={self.max_iter!r})"
            )

        # tol
        if not isinstance(self.tol, numbers.Number) or self.tol < 0:
            raise ValueError(
                "Tolerance for stopping criteria must be a positive number; got "
                f"(tol={self.tol!r})"
            )

        # graph_coeff
        if not isinstance(self.graph_coeff, numbers.Number) or self.graph_coeff < 0:
            raise ValueError(
                "Coefficient for the graph Laplacian regularization term must be "
                f"a positive number; got (graph_coeff={self.graph_coeff!r})"
            )

        # label_coeff
        if not isinstance(self.label_coeff, numbers.Number) or self.label_coeff < 0:
            raise ValueError(
                "Coefficient for the label information regularization term must be "
                f"a positive number; got (label_coeff={self.label_coeff!r})"
            )

        # knn
        self.knn_ = (
            clone(self.knn)
            if self.knn
            else NearestNeighbors(n_neighbors=min(5, X.shape[0]))
        )
        if not isinstance(self.knn_, NearestNeighbors):
            raise ValueError(
                "Pre-configured nearest-neighbors calculator must be of type"
                f" NearestNeighbors; got (knn={self.knn_!r})"
            )

        return self

    @property
    def _n_features_out(self):
        """Number of transformed output features."""
        return self.components_.shape[0]
