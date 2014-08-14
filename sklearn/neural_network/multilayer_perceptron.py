"""Multi-layer Perceptron
"""

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# Licence: BSD 3 clause

import numpy as np

from abc import ABCMeta, abstractmethod
from scipy.optimize import fmin_l_bfgs_b
import warnings

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..externals import six
from ..preprocessing import LabelBinarizer
from ..utils import gen_batches
from ..utils import shuffle
from ..utils import check_array, check_X_y, column_or_1d
from ..utils import ConvergenceWarning
from ..utils.extmath import safe_sparse_dot
from .base import logistic, softmax, init_weights
from .base import clear_layer_lists
from .base import ACTIVATIONS, DERIVATIVES, LOSS_FUNCTIONS


def _pack(layers_coef_, layers_intercept_):
    """Pack the parameters into a single vector."""
    return np.hstack([l.ravel() for l in layers_coef_ + layers_intercept_])


class BaseMultilayerPerceptron(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for MLP classification and regression.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """

    @abstractmethod
    def __init__(self, n_hidden, activation, algorithm,
                 alpha, batch_size, learning_rate, learning_rate_init, power_t,
                 max_iter, loss, shuffle, random_state, tol, verbose,
                 warm_start):
        self.activation = activation
        self.algorithm = algorithm
        self.alpha = alpha
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.learning_rate_init = learning_rate_init
        self.power_t = power_t
        self.max_iter = max_iter
        self.loss = loss
        self.n_hidden = n_hidden
        self.shuffle = shuffle
        self.random_state = random_state
        self.tol = tol
        self.verbose = verbose
        self.warm_start = warm_start

        self.layers_coef_ = None
        self.layers_intercept_ = None
        self.cost_ = None
        self.n_iter_ = None
        self.learning_rate_ = None
        self.classes_ = None

    def _unpack(self, packed_parameters):
        """Extract the coefficients and intercepts from packed_parameters."""
        for i in range(self.n_layers_ - 1):
            start, end, shape = self._coef_indptr[i]
            self.layers_coef_[i] = np.reshape(packed_parameters[start:end],
                                              shape)

            start, end = self._intercept_indptr[i]
            self.layers_intercept_[i] = packed_parameters[start:end]

    def _forward_pass(self, with_output_activation=True):
        """Perform a forward pass on the network by computing the values
        of the neurons in the hidden layers and the output layer.

        Parameters
        ----------
        with_output_activation : if True, the output passes through
                                 the output activation function, which
                                 is either the softmax function or the
                                 logistic function
        """
        last_layer = self.n_layers_ - 1

        # Iterate over the layers
        for i in range(last_layer):
            self._a_layers[i + 1] = safe_sparse_dot(self._a_layers[i],
                                                    self.layers_coef_[i])
            self._a_layers[i + 1] += self.layers_intercept_[i]

            # For the hidden layers
            if i + 1 != last_layer:
                activation = ACTIVATIONS[self.activation]
                self._a_layers[i + 1] = activation(self._a_layers[i + 1])

            # For the last layer
            elif with_output_activation:
                out_activation = ACTIVATIONS[self.out_activation_]
                self._a_layers[i + 1] = out_activation(self._a_layers[i + 1])

    def _compute_cost_grad(self, layer, n_samples):
        """Compute the cost gradient for the layer."""
        self._coef_grads[layer] = safe_sparse_dot(self._a_layers[layer].T,
                                                  self._deltas[layer])
        self._coef_grads[layer] += (self.alpha * self.layers_coef_[layer])
        self._coef_grads[layer] /= n_samples

        self._intercept_grads[layer] = np.mean(self._deltas[layer], 0)

    def _backprop_sgd(self, X, y):
        """Update the weights using the computed gradients.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Subset of the target values.

        """
        self._a_layers[0] = X

        cost = self._backprop(X, y)

        # update weights
        for i in range(self.n_layers_ - 1):
            self.layers_coef_[i] -= (self.learning_rate_ * self._coef_grads[i])
            self.layers_intercept_[i] -= (self.learning_rate_ *
                                          self._intercept_grads[i])

        if self.learning_rate == 'invscaling':
            self.learning_rate_ = self.learning_rate_init / \
                pow(self.n_iter_ + 1, self.power_t)

        return cost

    def _backprop_lbfgs(self, X, y):
        """Apply the quasi-Newton optimization methods that uses a l_BFGS
        to train the weights.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Subset of the target values.

        """
        packed_coef_inter = _pack(self.layers_coef_, self.layers_intercept_)

        if self.verbose is True or self.verbose >= 1:
            iprint = 1
        else:
            iprint = -1

        optimal_parameters, f, d = fmin_l_bfgs_b(func=self._cost_grad_lbfgs,
                                                 x0=packed_coef_inter,
                                                 maxfun=self.max_iter,
                                                 iprint=iprint,
                                                 pgtol=self.tol,
                                                 args=(X, y))
        self.cost_ = f
        self._unpack(optimal_parameters)

    def _cost_grad_lbfgs(self, packed_coef_inter, X, y):
        """Compute the MLP cost  function and its
        corresponding derivatives with respect to the
        different parameters given in the initialization.

        Parameters
        ----------
        packed_parameters : array-like
                    A vector comprising the  flattened coefficients and
                    intercepts

        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Subset of the target values.

        Returns
        -------
        cost : float
        grad : array-like, shape (number of nodes of all layers,)

        """
        self._unpack(packed_coef_inter)
        cost = self._backprop(X, y)

        self.n_iter_ += 1
        grad = _pack(self._coef_grads, self._intercept_grads)
        return cost, grad

    def _backprop(self, X, y):
        """Compute the MLP cost function and its corresponding derivatives
        with respect to each parameter: weights and bias vectors.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Subset of the target values.

        Returns
        -------
        cost : float
        """
        n_samples = X.shape[0]
        # Step (1/3): Forward propagate
        self._forward_pass()

        # Step (2/3): Get cost
        cost = LOSS_FUNCTIONS[self.loss](y, self._a_layers[-1])
        # Add L2 regularization term to cost
        values = np.sum(np.array([np.sum(s ** 2) for s in self.layers_coef_]))
        cost += (0.5 * self.alpha) * values / n_samples

        # Step (3/3): Backward propagate
        last = len(self.n_hidden)

        diff = y - self._a_layers[-1]
        self._deltas[last] = -diff

        # Compute gradient for the last layer
        self._compute_cost_grad(last, n_samples)

        # Iterate over the hidden layers
        for i in range(self.n_layers_ - 2, 0, -1):
            self._deltas[i - 1] = safe_sparse_dot(self._deltas[i],
                                                  self.layers_coef_[i].T)
            derivative = DERIVATIVES[self.activation]
            self._deltas[i - 1] *= derivative(self._a_layers[i])

            self._compute_cost_grad(i - 1, n_samples)

        return cost

    def _fit(self, X, y, warm_start=False):
        # Ensure self.n_hidden is a list
        if not hasattr(self.n_hidden, "__iter__"):
            self.n_hidden = [self.n_hidden]

        # Validate input parameters.
        if np.any(np.array(self.n_hidden) <= 0):
            raise ValueError("n_hidden must be > 0, got %s." % self.n_hidden)
        if not isinstance(self.shuffle, bool):
            raise ValueError("shuffle must be either True or False, got %s." %
                             self.shuffle)
        if self.max_iter <= 0:
            raise ValueError("max_iter must be > 0, got %s." % self.max_iter)
        if self.alpha < 0.0:
            raise ValueError("alpha must be >= 0, got %s." % self.alpha)
        if (self.learning_rate in ["constant", "invscaling"] and
           self.learning_rate_init <= 0.0):
            raise ValueError("learning_rate_init must be > 0, got %s." %
                             self.learning_rate)

        # raise ValueError if not registered
        if self.activation not in ACTIVATIONS:
            raise ValueError("The activation %s is not supported. Supported "
                             "activations are %s." % (self.activation,
                                                      ACTIVATIONS))
        if self.learning_rate not in ["constant", "invscaling"]:
            raise ValueError("learning rate %s is not supported. " %
                             self.learning_rate)
        if self.algorithm not in ["sgd", "l-bfgs"]:
            raise ValueError("The algorithm %s is not supported. " %
                             self.algorithm)

        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'],
                         multi_output=True)

        # This outputs a warning when a 1d array is expected
        if y.ndim == 2 and y.shape[1] == 1:
            y = column_or_1d(y, warn=True)

        n_samples, n_features = X.shape

        # Classification
        if isinstance(self, ClassifierMixin):
            self.label_binarizer_.fit(y)

            if self.classes_ is None or not warm_start:
                self.classes_ = self.label_binarizer_.classes_
            else:
                classes = self.label_binarizer_.classes_
                if not np.all(np.in1d(classes, self.classes_)):
                    raise ValueError("`y` has classes not in `self.classes_`."
                                     " `self.classes_` has %s. 'y' has %s." %
                                     (self.classes_, classes))

            y = self.label_binarizer_.transform(y)

        # Ensure y is 2D
        if y.ndim == 1:
            y = y.reshape((-1, 1))

        self.n_outputs_ = y.shape[1]

        layer_units = ([n_features] + self.n_hidden + [self.n_outputs_])

        # First time training the model
        if self.layers_coef_ is None or (not self.warm_start and
                                         not warm_start):
            # Initialize parameters
            self.n_iter_ = 0
            self.learning_rate_ = self.learning_rate_init
            self.n_outputs_ = y.shape[1]

            # Compute the number of layers
            self.n_layers_ = len(layer_units)

            # Output for regression
            if not isinstance(self, ClassifierMixin):
                self.out_activation_ = 'identity'
            # Output for multi class
            elif self.label_binarizer_.y_type_ == 'multiclass':
                self.out_activation_ = 'softmax'
            # Output for binary class and multi-label
            else:
                self.out_activation_ = 'logistic'

            # Initialize coefficient and intercept layers
            self.layers_coef_ = []
            self.layers_intercept_ = []

            for i in range(self.n_layers_ - 1):
                fan_in = layer_units[i]
                fan_out = layer_units[i + 1]

                # Get scaled weights based on Glorot et al. randomization
                if self.activation == 'tanh':
                    weight_scale = np.sqrt(6. / (fan_in + fan_out))

                elif self.activation == 'logistic':
                    weight_scale = 4. * np.sqrt(6. / (fan_in + fan_out))
                # For other activation functions
                else:
                    weight_scale = np.sqrt(1. / fan_in)

                coef_intercept_weights = init_weights(weight_scale, fan_in,
                                                      fan_out,
                                                      self.random_state)
                # Initialize weights
                self.layers_coef_.append(coef_intercept_weights[0])
                self.layers_intercept_.append(coef_intercept_weights[1])

        if self.shuffle:
            X, y = shuffle(X, y, random_state=self.random_state)

        # l-bfgs does not support mini-batches
        if self.algorithm == 'l-bfgs':
            batch_size = n_samples
        else:
            batch_size = np.clip(self.batch_size, 1, n_samples)

        # Initialize lists
        self._a_layers = [X]
        self._a_layers.extend(np.empty((batch_size, layer_units[i + 1]))
                              for i in range(self.n_layers_ - 1))
        self._deltas = [np.empty((batch_size, layer_units[i + 1]))
                        for i in range(self.n_layers_ - 1)]
        self._coef_grads = [np.empty((layer_units[i], layer_units[i + 1]))
                            for i in range(self.n_layers_ - 1)]
        self._intercept_grads = [np.empty(layer_units[i + 1])
                                 for i in range(self.n_layers_ - 1)]

        # Run the Stochastic Gradient Descent algorithm
        if self.algorithm == 'sgd':
            prev_cost = np.inf
            cost_increase_count = 0

            for i in range(self.max_iter):
                for batch_slice in gen_batches(n_samples, batch_size):
                    self.cost_ = self._backprop_sgd(X[batch_slice],
                                                    y[batch_slice])

                self.n_iter_ += 1

                if self.verbose:
                    print("Iteration %d, cost = %.8f" % (i, self.cost_))

                if self.cost_ > prev_cost:
                    cost_increase_count += 1
                    if cost_increase_count == 0.2 * self.max_iter:
                        warnings.warn('Cost is increasing for more than 20%%'
                                      ' of the iterations. Consider reducing'
                                      ' learning_rate_init and preprocessing'
                                      ' your data with StandardScaler or '
                                      ' MinMaxScaler.'
                                      % self.cost_, ConvergenceWarning)

                elif prev_cost - self.cost_ < self.tol or warm_start:
                    break

                prev_cost = self.cost_

        # Run the LBFGS algorithm
        elif self.algorithm == 'l-bfgs':
            # Store meta information for the parameters
            self._coef_indptr = []
            self._intercept_indptr = []
            start = 0

            # Save sizes and indices of coefficients for faster unpacking
            for i in range(self.n_layers_ - 1):
                fan_in, fan_out = layer_units[i], layer_units[i + 1]

                end = start + (fan_in * fan_out)
                self._coef_indptr.append((start, end, (fan_in, fan_out)))
                start = end

            # Save sizes and indices of intercepts for faster unpacking
            for i in range(self.n_layers_ - 1):
                end = start + layer_units[i + 1]
                self._intercept_indptr.append((start, end))
                start = end

            # Run LBFGS
            self._backprop_lbfgs(X, y)

        # Clear the lists
        clear_layer_lists(self._a_layers, self._deltas, self._coef_grads,
                          self._intercept_grads)

        return self

    def fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Target values.

        Returns
        -------
        self : returns an instance of self.
        """
        return self._fit(X, y, warm_start=False)

    def partial_fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Subset of training data.

        y : array-like, shape (n_samples,)
            Subset of target values.

        Returns
        -------
        self : returns an instance of self.
        """
        if self.algorithm != 'sgd':
            raise ValueError("only SGD algorithm supports partial fit")

        return self._fit(X, y, warm_start=True)

    def _decision_scores(self, X):
        """Predict using the trained model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y_pred : array-like, shape (n_samples,) or (n_samples, n_outputs)
                 The predicted values.
        """
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])

        layer_units = [X.shape[1]] + self.n_hidden + [self.n_outputs_]

        # Initialize layers
        self._a_layers = []
        self._a_layers.append(X)

        for i in range(self.n_layers_ - 1):
            self._a_layers.append(np.empty((X.shape[0],
                                            layer_units[i + 1])))
        # forward propagate
        self._forward_pass(with_output_activation=False)
        y_pred = self._a_layers[-1]

        # Clear the list
        clear_layer_lists(self._a_layers)

        return y_pred


class MultilayerPerceptronClassifier(BaseMultilayerPerceptron,
                                     ClassifierMixin):
    """Multi-layer Perceptron classifier.

    Under a loss function, the algorithm trains either by l-bfgs or gradient
    descent. The training is iterative, in that at each time step the
    partial derivatives of the loss function with respect to the model
    parameters are computed to update the parameters.

    It has a regularizer as a penalty term added to the loss function that
    shrinks model parameters towards zero.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    Parameters
    ----------
    n_hidden : python list, length = n_layers - 2, default [100]
        A list containing the number of neurons in the ith hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'tanh'
        Activation function for the hidden layer. It only applies to
        kernel='random'.

         - 'logistic' returns f(x) = 1 / (1 + exp(x)).

         - 'tanh' returns f(x) = tanh(x).

         - 'relu' returns f(x) = max(0, x)

    algorithm : {'l-bfgs', 'sgd'}, default 'l-bfgs'
        The algorithm for weight optimization.  Defaults to 'l-bfgs'

        - 'l-bfgs' is an optimization algorithm in the family of
          quasi-Newton methods.

        - 'sgd' refers to stochastic gradient descent.

    alpha : float, optional, default 0.00001
        L2 penalty (regularization term) parameter.

    batch_size : int, optional, default 200
        Size of minibatches in SGD optimizer.
        If you select the algorithm as 'l-bfgs',
        then the classifier will not use minibatches.

    learning_rate : {'constant', 'invscaling'}, default 'constant'
        Base learning rate for weight updates.

        -'constant', as it stands,  keeps the learning rate
         'learning_rate_init' constant throughout training.
         learning_rate_ = learning_rate_init

        -'invscaling' gradually decreases the learning rate 'learning_rate_' at
          each time step 't' using an inverse scaling exponent of 'power_t'.
          learning_rate_ = learning_rate_init / pow(t, power_t)

    max_iter : int, optional, default 200
        Maximum number of iterations. The algorithm
        iterates until convergence (determined by 'tol') or
        this number of iterations.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    shuffle : bool, optional, default False
        Whether to shuffle samples in each iteration before extracting
        minibatches.

    tol : float, optional, default 1e-5
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached and the algorithm exits.

    learning_rate_init : double, optional, default 0.5
        The initial learning rate used. It controls the step-size
        in updating the weights.

    power_t : double, optional, default 0.5
        The exponent for inverse scaling learning rate.
        It is used in updating learning_rate_init when the learning_rate
        is set to 'invscaling'.

    verbose : bool, optional, default False
        Whether to print progress messages to stdout.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.

    Attributes
    ----------
    `classes_` : array or list of array of shape (n_classes,)
        Class labels for each output.

    `cost_` : float
        The current cost value computed by the loss function.

    `label_binarizer_` : LabelBinarizer
        A LabelBinarizer object trained on the training set.

    `layers_coef_` : python list, length n_layers - 1
        The ith element in the list represents the weight matrix corresponding
        to layer i.

    `layers_intercept_` : python list, length n_layers - 1
        The ith element in the list represents the bias vector corresponding to
        layer i + 1.

    `learning_rate_` : float
        The current learning rate.

    n_iter_ : int,
        The current number of iterations the algorithm has ran.

    n_layers_ : int
        Number of layers.

    `n_outputs_` : int
        Number of outputs.

    `out_activation_` : string
        Name of the output activation function.

    References
    ----------
    Hinton, Geoffrey E.
        "Connectionist learning procedures." Artificial intelligence 40.1
        (1989): 185-234.

    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """
    def __init__(self, n_hidden=[100], activation="tanh",
                 algorithm='l-bfgs', alpha=0.00001,
                 batch_size=200, learning_rate="constant",
                 learning_rate_init=0.5, power_t=0.5, max_iter=200,
                 shuffle=False, random_state=None, tol=1e-5,
                 verbose=False, warm_start=False):

        sup = super(MultilayerPerceptronClassifier, self)
        sup.__init__(n_hidden=n_hidden, activation=activation,
                     algorithm=algorithm, alpha=alpha, batch_size=batch_size,
                     learning_rate=learning_rate,
                     learning_rate_init=learning_rate_init, power_t=power_t,
                     max_iter=max_iter, loss='log_loss', shuffle=shuffle,
                     random_state=random_state, tol=tol,
                     verbose=verbose, warm_start=warm_start)

        self.label_binarizer_ = LabelBinarizer()

    def decision_function(self, X):
        """Decision function of the elm model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array-like, shape (n_samples,) or (n_samples, n_classes)
            The predict values.
        """
        y_scores = self._decision_scores(X)

        if self.n_outputs_ == 1:
            return y_scores.ravel()
        else:
            return y_scores

    def predict(self, X):
        """Predict using the extreme learning machines model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array-like, shape (n_samples,) or (n_samples, n_classes)
            The predicted classes, or the predict values.
        """
        y_scores = self.decision_function(X)
        y_scores = ACTIVATIONS[self.out_activation_](y_scores)

        return self.label_binarizer_.inverse_transform(y_scores)

    def partial_fit(self, X, y, classes=None):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        y : array-like, shape (n_samples,)
            Subset of the target values.

        classes : array, shape (n_classes)
            Classes across all calls to partial_fit.
            Can be obtained by via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.
            Note that y doesn't need to contain all labels in `classes`.

        Returns
        -------
        self : returns an instance of self.
        """
        self.classes_ = classes

        super(MultilayerPerceptronClassifier, self).partial_fit(X, y)

        return self

    def predict_log_proba(self, X):
        """Return the log of probability estimates.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y_prob : array-like, shape (n_samples, n_classes)
                 The predicted log-probability of the sample for each class
                 in the model, where classes are ordered as they are in
                 `self.classes_`. Equivalent to log(predict_proba(X))
        """
        y_prob = self.predict_proba(X)
        return np.log(y_prob, out=y_prob)

    def predict_proba(self, X):
        """Probability estimates.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y_prob : array-like, shape (n_samples, n_classes)
                 The predicted probability of the sample for each class in the
                 model, where classes are ordered as they are in
                 `self.classes_`.
        """
        y_scores = self.decision_function(X)

        if y_scores.ndim == 1:
            y_scores = logistic(y_scores)
            return np.vstack([1 - y_scores, y_scores]).T
        else:
            return softmax(y_scores)


class MultilayerPerceptronRegressor(BaseMultilayerPerceptron, RegressorMixin):
    """Multi-layer Perceptron regressor.

    Under a loss function, the algorithm trains either by l-bfgs or gradient
    descent. The training is iterative, in that at each time step the
    partial derivatives of the loss function with respect to the model
    parameters are computed to update the parameters.

    It has a regularizer as a penalty term added to the loss function that
    shrinks model parameters towards zero.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    Parameters
    ----------
    n_hidden : python list, length = n_layers - 2, default [100]
        A list containing the number of neurons in the ith hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'tanh'
        Activation function for the hidden layer. It only applies to
        kernel='random'.

         - 'logistic' returns f(x) = 1 / (1 + exp(x)).

         - 'tanh' returns f(x) = tanh(x).

         - 'relu' returns f(x) = max(0, x)

    algorithm : {'l-bfgs', 'sgd'}, default 'l-bfgs'
        The algorithm for weight optimization.  Defaults to 'l-bfgs'

        - 'l-bfgs' is an optimization algorithm in the family of
          quasi-Newton methods.

        - 'sgd' refers to stochastic gradient descent.

    alpha : float, optional, default 0.00001
        L2 penalty (regularization term) parameter.

    batch_size : int, optional, default 200
        Size of minibatches in SGD optimizer.
        If you select the algorithm as 'l-bfgs',
        then the classifier will not use minibatches.

    learning_rate : {'constant', 'invscaling'}, default 'constant'
        Base learning rate for weight updates.

        -'constant', as it stands,  keeps the learning rate
         'learning_rate_init' constant throughout training.
         learning_rate_ = learning_rate_init

        -'invscaling' gradually decreases the learning rate 'learning_rate_' at
          each time step 't' using an inverse scaling exponent of 'power_t'.
          learning_rate_ = learning_rate_init / pow(t, power_t)

    max_iter : int, optional, default 200
        Maximum number of iterations. The algorithm
        iterates until convergence (determined by 'tol') or
        this number of iterations.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    shuffle : bool, optional, default False
        Whether to shuffle samples in each iteration before extracting
        minibatches.

    tol : float, optional, default 1e-5
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached and the algorithm exits.

    learning_rate_init : double, optional, default 0.5
        The initial learning rate used. It controls the step-size
        in updating the weights.

    power_t : double, optional, default 0.5
        The exponent for inverse scaling learning rate.
        It is used for updating learning_rate_ when it
        is set to 'invscaling'.

    verbose : bool, optional, default False
        Whether to print progress messages to stdout.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.

    Attributes
    ----------
    `cost_` : float
        The current cost value computed by the loss function.

    `layers_coef_` : python list, length n_layers - 1
        The ith element in the list represents the weight matrix corresponding
        to layer i.

    `layers_intercept_` : python list, length n_layers - 1
        The ith element in the list represents the bias vector corresponding to
        layer i + 1.

    `learning_rate_` : float
        The current learning rate.

    n_iter_ : int,
        The current number of iterations the algorithm has ran.

    n_layers_ : int
        Number of layers.

    `n_outputs_` : int
        Number of outputs.

    `out_activation_` : string
        Name of the output activation function.

    References
    ----------
    Hinton, Geoffrey E.
        "Connectionist learning procedures." Artificial intelligence 40.1
        (1989): 185-234.

    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """
    def __init__(self, n_hidden=[100], activation="tanh",
                 algorithm='l-bfgs', alpha=0.00001,
                 batch_size=200, learning_rate="constant",
                 learning_rate_init=0.1,
                 power_t=0.5, max_iter=100, shuffle=False,
                 random_state=None, tol=1e-5,
                 verbose=False, warm_start=False):

        sup = super(MultilayerPerceptronRegressor, self)
        sup.__init__(n_hidden=n_hidden, activation=activation,
                     algorithm=algorithm, alpha=alpha, batch_size=batch_size,
                     learning_rate=learning_rate,
                     learning_rate_init=learning_rate_init, power_t=power_t,
                     max_iter=max_iter, loss='squared_loss', shuffle=shuffle,
                     random_state=random_state, tol=tol,
                     verbose=verbose, warm_start=warm_start)

    def predict(self, X):
        """Predict using the multi-layer perceptron model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Data, where n_samples is the number of samples
            and n_features is the number of features.

        Returns
        -------
        y : array-like, shape (n_samples, n_outputs)
            The predicted classes, or the predict values.
        """
        return self._decision_scores(X)
