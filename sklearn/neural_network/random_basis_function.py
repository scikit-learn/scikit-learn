import numpy as np

from base import ACTIVATIONS

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
    intercept : boolean, default True
        Whether to randomly generate an intercept. 

    weight_scale : float, default 'auto'
        Initializes and scales the input-to-hidden weights.
        The weight values will range between plus and minus
        'sqrt(weight_scale * 6. / (n_features + n_hidden))' based on the
        uniform distribution.

    n_activated_features : int, default 10
        The number of units in the hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'tanh'
        Activation function for the hidden layer.

         - 'logistic', the logistic sigmoid function,
            returns f(x) = 1 / (1 + exp(x)).

         - 'tanh', the hyperbolic tan function,
            returns f(x) = tanh(x).

         - 'relu', the rectified linear unit function,
            returns f(x) = max(0, x).

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    Attributes
    ----------
    `coef_` : array-like, shape (n_features, n_hidden)
        The input-to-hidden weights.

    `intercept_` : array-like, shape (n_hidden,)
        The bias added to the hidden layer neurons.

    References
    ----------
    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """

    def __init__(self, n_activated_features=10, weight_scale='auto',
                 activation='tanh', intercept=True, random_state=None):
        self.n_activated_features = n_activated_features
        self.weight_scale = weight_scale
        self.intercept = intercept
        self.activation = activation
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
        if self.n_activated_features <= 0:
            raise ValueError("n_activated_features must be > 0, got %s." %
                             self.n_activated_features)

        if self.activation not in ACTIVATIONS:
            raise ValueError("The activation %s is not supported. Supported "
                             "activations are %s." % (self.activation,
                                                      ACTIVATIONS))

        X = check_array(X, accept_sparse=['csr', 'csc'])

        n_samples, n_features = X.shape

        rng = check_random_state(self.random_state)

        if self.weight_scale == 'auto':
            weight_init_bound = np.sqrt(6. / (n_features +
                                              self.n_activated_features))
        else:
            weight_init_bound = self.weight_scale

        self.coef_ = rng.uniform(-weight_init_bound, weight_init_bound,
                                 (n_features, self.n_activated_features))

        if self.intercept:
            self.intercept_ = rng.uniform(-weight_init_bound, weight_init_bound,
                                          self.n_activated_features)

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
