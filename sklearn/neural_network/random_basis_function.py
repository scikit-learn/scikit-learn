import numpy as np

from ._base import ACTIVATIONS

from ..base import BaseEstimator, TransformerMixin
from ..utils import check_random_state
from ..utils.extmath import safe_sparse_dot
from ..utils.validation import check_array


class RandomBasisFunction(BaseEstimator, TransformerMixin):
    """Random basis activation.

    The algorithm uses the number of features of the input data to randomly 
    generate coefficient parameters based on the uniform probability 
    distribution.

    Using these coefficient parameters, the algorithm can transform 
    the data into a different dimensional space.

    Parameters
    ----------
    n_outputs : int, default 10
        The number of output features to generate.

    weight_scale : float, default 'auto'
        If 'auto', `coef_` and `intercept_` get values ranging between 
        plus/minus 'sqrt(6. / (n_features + n_outputs))' based on 
        the uniform distribution; otherwise, between +weight_scale and 
        -weight_scale.

    activation : {'logistic', 'tanh', 'relu'}, default 'tanh'
        Activation function for the output features.

         - 'logistic', the logistic sigmoid function,
            returns f(x) = 1 / (1 + exp(x)).

         - 'tanh', the hyperbolic tan function,
            returns f(x) = tanh(x).

         - 'relu', the rectified linear unit function,
            returns f(x) = max(0, x).

    intercept : boolean, default True
        Whether to randomly generate an intercept. 

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    Attributes
    ----------
    `coef_` : array-like, shape (n_features, n_outputs)
        The coefficient parameters used to generate the output features.

    `intercept_` : array-like, shape (n_outputs,)
        An intercept parameter added to the output features.

    References
    ----------
    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.

    Schmidt, Wouter F., Martin A. Kraaijveld, and Robert PW Duin.
        "Feedforward neural networks with random weights." Pattern Recognition,
        1992. Vol. II. Conference B: Pattern Recognition Methodology and Systems,
        Proceedings., 11th IAPR International Conference on. IEEE, 1992.

    See also
    --------
    `sklearn.random_projection` and `sklearn.kernel_approximation` contains algorithms
    that are similar to `RandomBasisFunction` in that they transform the input features to another
    dimensional space. However, `RandomBasisFunction` is more general in that
    the user defines the number of features to generate and the function to 
    apply on these output features.

    """
    def __init__(self, n_outputs=10, weight_scale='auto',
                 activation='tanh', intercept=True, random_state=None):
        self.n_outputs = n_outputs
        self.weight_scale = weight_scale
        self.activation = activation
        self.intercept = intercept
        self.random_state = random_state

    def fit(self, X, y=None):
        """
        Generate random parameters based on the number of features the input
        data has.

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape (n_samples, n_features).

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        self

        """
        # Sanity checks
        if self.n_outputs <= 0:
            raise ValueError("n_outputs must be > 0, got %s." %
                             self.n_outputs)

        if self.activation not in ACTIVATIONS:
            raise ValueError("The activation %s is not supported. Supported "
                             "activations are %s." % (self.activation,
                                                      ACTIVATIONS))

        X = check_array(X, accept_sparse=['csr', 'csc'])

        n_samples, n_features = X.shape

        rng = check_random_state(self.random_state)

        if self.weight_scale == 'auto':
            weight_init_bound = np.sqrt(6. / (n_features + self.n_outputs))
        else:
            weight_init_bound = self.weight_scale

        self.coef_ = rng.uniform(-weight_init_bound, weight_init_bound,
                                 (n_features, self.n_outputs))

        if self.intercept:
            self.intercept_ = rng.uniform(-weight_init_bound, weight_init_bound,
                                          self.n_outputs)

        self.activation_function = ACTIVATIONS[self.activation]

        return self

    def transform(self, X, y=None):
        """
        Transform input data to another space using the randomly generated 
        parameters

        Parameters
        ----------
        X : numpy array or scipy.sparse of shape (n_samples, n_features).

        y : is not used: placeholder to allow for usage in a Pipeline.

        Returns
        -------
        X_new : numpy array or scipy sparse of shape [n_samples, n_components]
            Projected array.

        """
        X = check_array(X, accept_sparse=['csr', 'csc'])

        X_new = safe_sparse_dot(X, self.coef_)

        if self.intercept:
            X_new += self.intercept_

        return self.activation_function(X_new)
