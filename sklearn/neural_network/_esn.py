"""Echo State Network
"""

# Authors: Giovanni Azua Garcia
# License: BSD 3 clause

import numpy as np
import pandas as pd
from numbers import Integral, Real

from ..base import RegressorMixin
from ..utils._param_validation import Interval
from ..utils.validation import check_is_fitted


class EchoStateNetwork(RegressorMixin):
    """Echo State Network (ESN).

    # TODO: document
    Echo State Network implementation ...

    # TODO: document
    The time complexity of this implementation is ``O(d ** 2)`` assuming
    d ~ n_features ~ n_components.

    Read more in the :ref:`User Guide <rbm>`.

    # TODO: document
    Parameters
    ----------
    n_components : int, default=256
        Number of binary hidden units.

    # TODO: document
    Attributes
    ----------
    reservoir_size_ : array-like of shape (n_components,)
        Biases of the hidden units.

    See Also
    --------
    sklearn.neural_network.MLPRegressor : Multi-layer Perceptron regressor.
    sklearn.neural_network.MLPClassifier : Multi-layer Perceptron classifier.
    sklearn.decomposition.PCA : An unsupervised linear dimensionality
        reduction model.

    # TODO: document
    References
    ----------

    [1] Hinton, G. E., Osindero, S. and Teh, Y. A fast learning algorithm for
        deep belief nets. Neural Computation 18, pp 1527-1554.
        https://www.cs.toronto.edu/~hinton/absps/fastnc.pdf

    # TODO: document
    Examples
    --------
    """

    _parameter_constraints: dict = {
        "reservoir-size": [Interval(Integral, 50, None, closed="left")],
        "leaking-rate": [Interval(Real, 0, None, closed="neither")],
        "spectral-radius": [Interval(Real, 0, None, closed="neither")],
        "input-scaling": [Interval(Real, 0, None, closed="neither")],
        "regularization": [Interval(Real, 0, None, closed="neither")],
        "num-warmup-steps": [Interval(Integral, 0, None, closed="left")],
        "verbose": ["verbose"],
        "random_state": ["random_state"]
    }

    def __init__(self, W_in, W, W_fb, solver='ridge', regularization=1e-8, leaking_rate=0.25, spectral_radius=1.25,
                 input_scaling=1.0, sparsity_fraction=0.0, num_warn_up_steps=50):
        """Constructor

        Parameters
        ----------
        W_in : ndarray first block of input weights. in_size and res_size are derived from the dimensions of this matrix.

        W: ndarray second block of input weights, must be pre-normalized by the largest eigenvalue.

        W_fb: ndarray output feedback weights.

        solver: solver implementation.

        regularization: float if applicable the regularization parameter

        leaking_rate: float the leaking rate.

        spectral_radius: float the spectral radius.

        input_scaling: float the input scaling.

        sparsity_fraction: float reservoir sparsity fraction.

        param num_warm_up_steps: int number of warm up steps.
        """
        super().__init__()

        self._in_size = W_in.shape[1] - 1
        self._res_size = W.shape[0]

        # run assertions
        assert W.shape[0] == W.shape[1]
        assert W.shape[0] == W_in.shape[0]
        assert W_in.shape[1] == self._in_size + 1
        assert W_fb.shape[0] == self._res_size
        assert not np.isna(W).any().any()
        assert not np.isna(W_in).any().any()
        assert not np.isna(W_fb).any().any()
        assert not pd.isnull(leaking_rate)
        assert not pd.isnull(spectral_radius) and spectral_radius > 0.0
        assert not not pd.isnull(input_scaling) and input_scaling > 0.0

        self._W_in = W_in

        # apply spectral radius once
        W *= spectral_radius
        self._W = W

        self._W_fb = W_fb

        self._solver = solver
        self._regularization = regularization
        self._leaking_rate = leaking_rate
        self._spectral_radius = spectral_radius
        self._input_scaling = input_scaling
        self._sparsity_fraction = sparsity_fraction
        self._num_warm_up_steps = num_warn_up_steps

    def fit(self, X, y=None):
        """Fit the model to the data X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs), default=None
            Target values (None for unsupervised transformations).

        Returns
        -------
        self : EchoStateNetwork
            The fitted model.
        """
        pass

    def predict(self, X):
        """Predict using the multi-layer perceptron model.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : ndarray of shape (n_samples, n_outputs)
            The predicted values.
        """
        check_is_fitted(self)
        return self._predict(X)

    def predict(self, X):
        """
        TODO: implement

        Parameters
        ----------
        X

        Returns
        -------

        """
        pass
