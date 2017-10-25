import numpy as np
import scipy as sp

from sklearn.covariance import empirical_covariance
from ..base import BaseEstimator, TransformerMixin
from ..utils.validation import check_array, check_is_fitted


class SFA(BaseEstimator, TransformerMixin):
    """ Slow Feature Analysis (SFA)

    Linear projection extracting the slowly-varying components from the input
    data. This implementation uses the formulation as generalized eigenvalue
    problem introduced in Berkes and Wiskott, 2005.

    More information about Slow Feature Analysis can be found in
    Wiskott, L. and Sejnowski, T.J., Slow Feature Analysis: Unsupervised
    Learning of Invariances, Neural Computation, 14(4):715-770 (2002).

    Parameters
    ----------

    copy : bool (default True)
        If False, data passed to fit are overwritten and running
        fit(X).transform(X) will not yield the expected results,
        use fit_transform(X) instead.

    Attributes
    ----------
    n_input_dim_ : int
        The number of input dimensions of the data passed to :meth:`fit`
    """
    def __init__(self, copy=True):
        self.copy = copy

    def fit(self, X, y=None):
        """Fit the SFA parameters with the `X`.

        Parameters
        ----------
        X : array-like, shape (n_training_samples, n_input_dim)
            The training input samples.
        y : None
            There is no need of a target in a transformer, yet the pipeline API
            requires this parameter.

        Returns
        -------
        self : SFA
            Returns self
        """
        X = check_array(X, dtype=[np.float64, np.float32], copy=self.copy,
                        ensure_min_samples=2)
        dX_dt = X[1:, :] - X[:-1, :]
        self.n_input_dim_ = X.shape[1]

        self.mean_ = np.mean(X, axis=0)
        X -= self.mean_

        cov_mtx = empirical_covariance(X, assume_centered=True)
        dcov_mtx = empirical_covariance(dX_dt, assume_centered=True)

        d, self.components_ = sp.linalg.eigh(dcov_mtx, cov_mtx)
        self.bias_ = np.dot(self.mean_, self.components_)

        # Return the transformer
        return self

    def transform(self, X):
        """ Project the input on the slow feature components.

        X is projected on the first show feature components previous extracted
        from the training data.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_input_dim)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)
        """
        # Check that fit has been called
        check_is_fitted(self, ['n_input_dim_', 'components_', 'bias_'])

        # Input validation
        X = check_array(X, dtype=[np.float64, np.float32],
                        ensure_min_samples=2)

        # Check that the input dimenstionality is the same as the one passed
        # during fit
        if X.shape[1] != self.n_input_dim_:
            raise ValueError(
                'Dimensionality of input is different from that seen '
                'in `fit`')

        X_transformed = np.dot(X, self.components_) - self.bias_
        return X_transformed
