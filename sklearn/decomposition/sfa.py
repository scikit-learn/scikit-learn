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

    n_components : int or None
        Number of components to keep. SFA will keep the `n_components` slowest
        components.
        If `n_components` is `None`, all components are kept, and the output
        dim is the same as the input dim.

    Attributes
    ----------
    n_input_dim_ : int
        The number of input dimensions of the data passed to :meth:`fit`

    components_ : array, shape (n_components, n_features)
        Slowest directions in inputs space. The components are sorted
        by slowest to fastest.

    mean_ : array, shape (n_components)
        Per-feature empirical mean, estimated from the training set.
        Equal to `X.mean(axis=0)`.

    bias_ : array, shape (n_components)
        Displacement necessary to make sure that the output signals have
        zero mean. This is a convenience quantity that is equal to
        `np.dot(self.mean_, self.components_)`.

    eta_ : array, shape (n_components)
        The eta value of the components. This is a measure of the slowness of
        a signal, defined as `eta(x) = 1/(2*pi) * sqrt(dx/dt)`. It can be
        interpreted as the number of oscillations per time step (or in
        other words in `n_samples` time points, a sine wave `x` makes
        `eta(x) * n_samples` oscillations).

    """
    def __init__(self, copy=True, n_components=None):
        self.copy = copy
        self.n_components = n_components

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
        # Time derivative
        dX_dt = X[1:, :] - X[:-1, :]
        self.n_input_dim_ = X.shape[1]
        self.mean_ = np.mean(X, axis=0)

        X -= self.mean_
        cov_mtx = empirical_covariance(X, assume_centered=True)

        # The time derivative is not centered around the mean:
        # we want the second moment matrix (centered around 0) and
        # not the second central moment matrix (centered around the mean), i.e.
        # the covariance matrix
        dt_cov_mtx = empirical_covariance(dX_dt, assume_centered=True)

        if self.n_components is None:
            eigvals = None
        else:
            eigvals = (0, self.n_components - 1)
        d, self.components_ = sp.linalg.eigh(dt_cov_mtx, cov_mtx,
                                             eigvals=eigvals)
        self.bias_ = np.dot(self.mean_, self.components_)
        self.eta_ = 0.5 / np.pi * np.sqrt(d)

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
        check_is_fitted(self, ['n_input_dim_', 'components_', 'bias_', 'eta_'])

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
