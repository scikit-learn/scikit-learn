# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#
# License: BSD Style.
"""Implementation of Stochastic Gradient Descent (SGD) with sparse data."""

import numpy as np
from scipy import sparse

from ...externals.joblib import Parallel, delayed
from ..base import BaseSGD
from ..sgd_fast_sparse import plain_sgd


class SGD(BaseSGD):
    """Linear model fitted by minimizing a regularized empirical loss with SGD

    SGD stands for Stochastic Gradient Descent: the gradient of the loss is
    estimated each sample at a time and the model is updated along the way with
    a decreasing strength schedule (aka learning rate).

    The regularizer is a penalty added to the loss function that shrinks model
    parameters towards the zero vector using either the squared euclidean norm
    L2 or the absolute norm L1 or a combination of both (Elastic Net). If the
    parameter update crosses the 0.0 value because of the regularizer, the
    update is truncated to 0.0 to allow for learning sparse models and achieve
    online feature selection.

    This implementation works on scipy.sparse X and dense coef_.

    Parameters
    ----------
    loss : str, 'hinge' or 'log' or 'modifiedhuber'
        The loss function to be used. Defaults to 'hinge'. The hinge loss is a
        margin loss used by standard linear SVM models. The 'log' loss is the
        loss of logistic regression models and can be used for probability
        estimation in binary classifiers. 'modifiedhuber' is another smooth loss
        that brings tolerance to outliers.

    penalty : str, 'l2' or 'l1' or 'elasticnet'
        The penalty (aka regularization term) to be used. Defaults to 'l2' which
        is the standard regularizer for linear SVM models. 'l1' and 'elasticnet'
        migh bring sparsity to the model (feature selection) not achievable with
        'l2'.

    alpha : float
        Constant that multiplies the regularization term. Defaults to 0.0001

    rho : float
        The Elastic Net mixing parameter, with 0 < rho <= 1.
        Defaults to 0.85.

    init_coef_ : array, shape = [n_features] if n_classes == 2 else [n_classes,
    n_features]
        The initial coeffients to warm-start the optimization.

    init_intercept_ : array, shape = [1] if n_classes == 2 else [n_classes]
        The initial intercept to warm-start the optimization.

    fit_intercept: bool
        Whether the intercept should be estimated or not. If False, the
        data is assumed to be already centered. Defaults to True.

    n_iter: int
        The number of passes over the training data (aka epochs).
        Defaults to 5.

    shuffle: bool
        Whether or not the training data should be shuffled after each epoch.
        Defaults to False.

    verbose: integer, optional
        The verbosity level

    n_jobs: integer, optional
        The number of CPUs to use to do the OVA (One Versus All, for
        multi-class problems) computation. -1 means 'all CPUs'. Defaults
        to 1.

    Attributes
    ----------
    `coef_` : array, shape = [n_features] if n_classes == 2 else [n_classes,
    n_features]
        Weights asigned to the features.

    `intercept_` : array, shape = [1] if n_classes == 2 else [n_classes]
        Constants in decision function.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> Y = np.array([1, 1, 2, 2])
    >>> from scikits.learn.sgd.sparse import SGD
    >>> clf = SGD()
    >>> clf.fit(X, Y)
    SGD(loss='hinge', init_coef_=None, shuffle=False, verbose=0,
      fit_intercept=True, n_iter=5, penalty='l2', n_jobs=1, rho=1.0,
      init_intercept_=None, alpha=0.0001)
    >>> print clf.predict([[-0.8, -1]])
    [ 1.]

    See also
    --------
    LinearSVC

    """

    def _set_coef(self, coef_):
        self.coef_ = coef_
        if coef_ is None:
            self.sparse_coef_ = None
        else:
            # sparse representation of the fitted coef for the predict method
            self.sparse_coef_ = sparse.csr_matrix(coef_)

    def fit(self, X, Y, **params):
        """Fit current model with Stochastic Gradient Descent

        X is expected to be a sparse matrix. For maximum efficiency, use a
        sparse matrix in CSR format (scipy.sparse.csr_matrix)
        """
        self._set_params(**params)
        X = sparse.csr_matrix(X)
        Y = np.asanyarray(Y, dtype=np.float64)
        # largest class id is positive class
        classes = np.unique(Y)
        self.classes = classes
        if classes.shape[0] > 2:
            self._fit_multiclass(X, Y)
        elif classes.shape[0] == 2:
            self._fit_binary(X, Y)
        else:
            raise ValueError("The number of class labels must be " \
                             "greater than one. ")
        # return self for chaining fit and predict calls
        return self

    def _fit_binary(self, X, Y):
        """Fit a binary classifier.
        """
        # encode original class labels as 1 (classes[1]) or -1 (classes[0]).
        Y_new = np.ones(Y.shape, dtype=np.float64) * -1.0
        Y_new[Y == self.classes[1]] = 1.0
        Y = Y_new

        n_samples, n_features = X.shape[0], X.shape[1]
        self.coef_ = np.zeros(n_features, dtype=np.float64, order="C")
        self.intercept_ = np.zeros(1, dtype=np.float64, order="C")
        if self.init_coef_ is not None:
            init_coef_ = np.asanyarray(self.init_coef_)
            if init_coef_.shape != (n_features,):
                raise ValueError("Provided init_coef_ does not match dataset.")
            self.coef_ = init_coef_
        if self.init_intercept_ is not None:
            init_intercept_ = np.asanyarray(self.init_intercept_)
            if init_intercept_.shape != (1,):
                raise ValueError("Provided init_intercept_ " \
                                 "does not match dataset.")
            else:
                self.intercept_ = init_intercept_

        X_data = np.array(X.data, dtype=np.float64, order="C")
        X_indices = X.indices
        X_indptr = X.indptr
        coef_, intercept_ = plain_sgd(self.coef_,
                                      self.intercept_,
                                      self.loss_function,
                                      self.penalty_type,
                                      self.alpha, self.rho,
                                      X_data,
                                      X_indices, X_indptr, Y,
                                      self.n_iter,
                                      int(self.fit_intercept),
                                      int(self.verbose),
                                      int(self.shuffle))

        # update self.coef_ and self.sparse_coef_ consistently
        self._set_coef(self.coef_)
        self.intercept_ = np.asarray(intercept_)

    def _fit_multiclass(self, X, Y):
        """Fit a multi-class classifier with a combination
        of binary classifiers, each predicts one class versus
        all others (OVA: One Versus All).
        """
        n_classes = self.classes.shape[0]
        n_samples, n_features = X.shape[0], X.shape[1]
        self.coef_ = np.zeros((n_classes, n_features),
                         dtype=np.float64, order="C")
        self.intercept_ = np.zeros(n_classes, dtype=np.float64, order="C")

        if self.init_coef_ is not None:
            init_coef_ = np.asanyarray(self.init_coef_)
            if init_coef_.shape != (n_classes, n_features):
                raise ValueError("Provided coef_ does not match dataset. ")
            else:
                self.coef_ = init_coef_

        if self.init_intercept_ is not None:
            init_intercept_ = np.asanyarray(self.init_intercept_)
            if init_intercept_.shape != (n_classes, ):
                raise ValueError("Provided init_intercept_ " \
                                 "does not match dataset.")
            else:
                self.intercept_ = init_intercept_

        X_data = np.array(X.data, dtype=np.float64, order="C")
        X_indices = X.indices
        X_indptr = X.indptr

        res = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
                delayed(_train_ova_classifier)(i, c, X_data, X_indices,
                                               X_indptr, Y, self.coef_[i],
                                               self.intercept_[i],
                                               self.loss_function,
                                               self.penalty_type, self.alpha,
                                               self.rho, self.n_iter,
                                               self.fit_intercept,
                                               self.verbose, self.shuffle)
            for i, c in enumerate(self.classes))

        for i, coef, intercept in res:
            self.coef_[i] = coef
            self.intercept_[i] = intercept

        self._set_coef(self.coef_)
        self.intercept_ = self.intercept_

    def predict_margin(self, X):
        """Predict signed 'distance' to the hyperplane (aka confidence score).

        Parameters
        ----------
        X : scipy.sparse matrix of shape [n_samples, n_features]

        Returns
        -------
        array, shape = [n_samples] if n_classes == 2 else [n_samples, n_classes]
          The signed 'distances' to the hyperplane(s).
        """
        # np.dot only works correctly if both arguments are sparse matrices
        if not sparse.issparse(X):
            X = sparse.csr_matrix(X)
        scores = np.asarray(np.dot(X, self.sparse_coef_.T).todense()
                            + self.intercept_)
        if self.classes.shape[0] == 2:
            return np.ravel(scores)
        else:
            return scores


def _train_ova_classifier(i, c, X_data, X_indices, X_indptr, Y, coef_,
                          intercept_, loss_function, penalty_type, alpha,
                          rho, n_iter, fit_intercept, verbose, shuffle):
    """Inner loop for One-vs.-All scheme"""
    Y_i = np.ones(Y.shape, dtype=np.float64) * -1.0
    Y_i[Y == c] = 1.0
    coef, intercept = plain_sgd(coef_, intercept_,
                                loss_function, penalty_type,
                                alpha, rho, X_data, X_indices,
                                X_indptr, Y_i, n_iter,
                                fit_intercept, verbose,
                                shuffle)
    return (i, coef, intercept)
