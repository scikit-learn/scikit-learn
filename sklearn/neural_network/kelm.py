import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils.multiclass import unique_labels
from sklearn.gaussian_process import kernels
import abc

# Authors: Carlos Perales <sir.perales@gmail.com>
# License: BSD 3 clause

kernel_dict = {}
for k_str in dir(kernels):
    k = getattr(kernels, k_str)
    if isinstance(k, abc.ABCMeta):
        kernel_dict[k_str.lower()] = k
kernel_dict['linear'] = kernel_dict['dotproduct']


class KernelELM(BaseEstimator, ClassifierMixin):
    """ Kernel implementation of Extreme Learning Machine.

    Parameters
    ----------
    C : float, optional (default=1.0)
        Penalty parameter C of the error term.

    kernel : string, optional (default='rbf')
         Specifies the kernel type to be used in the algorithm.
         It must be one of 'linear', 'poly', 'rbf', 'sigmoid', 'precomputed' or
         a callable.
         If none is given, 'rbf' will be used. If a callable is given it is
         used to pre-compute the kernel matrix from data matrices; that matrix
         should be an array of shape ``(n_samples, n_samples)``.

    gamma : float, optional (default='auto')
        Kernel coefficient for 'rbf', 'poly' and 'sigmoid'.
        If gamma is 'auto' then 1/n_features will be used instead.
    """
    def __init__(self, C=1.0, kernel='rbf', gamma='auto'):
        self.C = C
        self.kernel = kernel
        self.gamma = gamma

    def fit(self, X, y):
        """ Fitting method, which lies in obtaining beta and kernel data X_.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = [n_samples, n_features]
            The training input samples.
        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (class labels in classification, real numbers in
            regression).

        Returns
        -------
        self : object
            Returns self.
        """
        if self.C <= 0:
            raise ValueError('C must be positive')
        X, y = check_X_y(X, y)
        self.classes_ = unique_labels(y)

        # y must be encoded by zero arrays with 1 in a position assigned
        # to a label
        n = X.shape[0]
        if self.classes_.shape[0] == 1:
            self.n_classes_ = int(self.classes_[0] + 1)
            self.classes_ = np.arange(self.n_classes_)
        else:
            self.n_classes_ = len(self.classes_)

        T = np.zeros((n, self.n_classes_), dtype=np.float64)
        # It is essential in order to adapt it to string labels
        self.class_corr_ = {}
        for i in range(n):
            row = [y[i] == self.classes_]
            T[i] = np.array(row, dtype=np.float64)
            pos = np.argmax(row)
            self.class_corr_[str(pos)] = self.classes_[pos]
        # T is already encoded

        # Kernels
        if self.gamma == 'auto' and self.kernel == 'linear':
            self.gamma_ = 0.0
        elif self.gamma == 'auto' and self.kernel != 'linear':
            self.gamma_ = 1.0 / n
        else:
            self.gamma_ = self.gamma

        # Essential fitting
        if self.kernel == 'precomputed':
            self.gamma_ = 0.0
            omega_train = X
        else:
            try:  # Compatibility between unicode and str, for Python 2 and 3
                self.kernel_fun_ = kernel_dict[str(self.kernel)](self.gamma_)
            except KeyError:
                self.kernel_fun_ = self.kernel
            omega_train = self.kernel_fun_(X, X)

        alpha = np.eye(n) / self.C + omega_train
        self.beta_ = np.linalg.solve(alpha, T)
        self.h_ = self.beta_.shape[0]

        # Saving training data
        self.X_ = X
        self.y_ = y
        # Return the estimator
        return self

    def predict(self, X):
        """ Predicting method, which lies in applying kernel
        trick to Extreme Learning Machine.

        Parameters
        ----------
        X : array-like of shape = [n_samples, n_features]
            The input samples.

        Returns
        -------
        y : array of shape = [n_samples]
            Returns predicted targets.
        """
        check_is_fitted(self, ['X_', 'y_'])

        # Input validation
        try:
            X = check_array(X)
        except TypeError:
            raise ValueError('Predict with sparse input '
                             'when trained with dense')

        # Fitting
        if self.kernel == 'precomputed':
            omega_test = X.T
        else:
            omega_test = self.kernel_fun_(self.X_, X)

        indicator = np.dot(omega_test.T, self.beta_)

        # Decoding
        y = []
        for y_i in indicator:
            pos = np.argmax(y_i)
            y.append(self.class_corr_[str(pos)])
        return np.array(y)
