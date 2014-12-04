"""Multi-layer Perceptron
"""

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# Licence: BSD 3 clause

import numpy as np

from abc import ABCMeta, abstractmethod
from scipy.optimize import fmin_l_bfgs_b
import warnings

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from .base import logistic, softmax
from .base import ACTIVATIONS, DERIVATIVES, LOSS_FUNCTIONS
from ..externals import six
from ..preprocessing import LabelBinarizer
from ..utils import gen_batches, check_random_state
from ..utils import shuffle
from ..utils import check_array, check_X_y, column_or_1d
from ..utils import ConvergenceWarning
from ..utils.extmath import safe_sparse_dot


def _pack(layers_coef_, layers_intercept_):
    """Pack the parameters into a single vector."""
    return np.hstack([l.ravel() for l in layers_coef_ + layers_intercept_])


class BaseMultilayerPerceptron(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for MLP classification and regression.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """

    @abstractmethod
    def __init__(self, hidden_layer_sizes, activation, algorithm,
                 alpha, batch_size, learning_rate, learning_rate_init, power_t,
                 max_iter, loss, shuffle, random_state, tol, verbose,
                 warm_start, momentum):
        self.activation = activation
        self.algorithm = algorithm
        self.alpha = alpha
        self.batch_size = batch_size
        self.learning_rate = learning_rate
        self.learning_rate_init = learning_rate_init
        self.power_t = power_t
        self.max_iter = max_iter
        self.loss = loss
        self.hidden_layer_sizes = hidden_layer_sizes
        self.shuffle = shuffle
        self.random_state = random_state
        self.tol = tol
        self.verbose = verbose
        self.warm_start = warm_start
        self.momentum = momentum

        self.layers_coef_ = None
        self.layers_intercept_ = None
        self.cost_ = None
        self.n_iter_ = None
        self.learning_rate_ = None
        self.classes_ = None

        # iteration count for learning rate schedule
        # must not be int (e.g. if ``learning_rate=='optimal'``)
        self.t_ = None

    def _unpack(self, packed_parameters):
        """Extract the coefficients and intercepts from packed_parameters."""
        for i in range(self.n_layers_ - 1):
            start, end, shape = self._coef_indptr[i]
            self.layers_coef_[i] = np.reshape(packed_parameters[start:end],
                                              shape)

            start, end = self._intercept_indptr[i]
            self.layers_intercept_[i] = packed_parameters[start:end]

    def _forward_pass(self, activations, with_output_activation=True):
        """Perform a forward pass on the network by computing the values
        of the neurons in the hidden layers and the output layer.

        Parameters
        ----------
        activations: list, length = n_layers - 1
            The ith index of the list holds the values of the ith layer.

        with_output_activation : bool, default True
            If True, the output passes through the output activation
            function, which is either the softmax function or the
            logistic function
        """
        # Iterate over the hidden layers
        for i in range(self.n_layers_ - 1):
            activations[i + 1] = safe_sparse_dot(activations[i],
                                                 self.layers_coef_[i])
            activations[i + 1] += self.layers_intercept_[i]

            # For the hidden layers
            if i + 1 != self.n_layers_ - 1:
                hidden_activation = ACTIVATIONS[self.activation]
                activations[i + 1] = hidden_activation(activations[i + 1])

        # For the last layer
        if with_output_activation:
            output_activation = ACTIVATIONS[self.out_activation_]
            activations[i + 1] = output_activation(activations[i + 1])

        return activations

    def _compute_cost_grad(self, layer, n_samples, activations, deltas,
                           coef_grads, intercept_grads):
        """Compute the cost gradient for the layer."""
        coef_grads[layer] = safe_sparse_dot(activations[layer].T,
                                            deltas[layer])
        coef_grads[layer] += (self.alpha * self.layers_coef_[layer])
        coef_grads[layer] /= n_samples

        intercept_grads[layer] = np.mean(deltas[layer], 0)

        return coef_grads, intercept_grads

    def _cost_grad_lbfgs(self, packed_coef_inter, X, y, activations, deltas,
                         coef_grads, intercept_grads):
        """Compute the MLP cost  function and its corresponding derivatives
        with respect to the different parameters given in the initialization.

        Parameters
        ----------
        packed_parameters : array-like
            A vector comprising the  flattened coefficients and intercepts.

        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            The target values.

        activations: list, length = n_layers - 1
            The ith index of the list holds the values of the ith layer.

        deltas : list, length = n_layers - 1
            The ith index of the list holds the difference between the
            activations of the i + 1 layer and the backpropagated error.

        coef_grad : list, length = n_layers - 1
            The ith index contains the amount of change used to update the
            coefficient parameters of the ith layer in an iteration.

        intercept_grads : list, length = n_layers - 1
            The ith index contains the amount of change used to update the
            intercept parameters of the ith layer in an iteration.

        Returns
        -------
        cost : float
        grad : array-like, shape (number of nodes of all layers,)

        """
        self._unpack(packed_coef_inter)
        cost, coef_grads, intercept_grads = self._backprop(X, y, activations,
                                                           deltas, coef_grads,
                                                           intercept_grads)
        self.n_iter_ += 1
        grad = _pack(coef_grads, intercept_grads)
        return cost, grad

    def _backprop(self, X, y, activations, deltas, coef_grads,
                  intercept_grads):
        """Compute the MLP cost function and its corresponding derivatives
        with respect to each parameter: weights and bias vectors.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            The target values.

        activations: list, length = n_layers - 1
             The ith index of the list holds the values of the ith layer.

        deltas : list, length = n_layers - 1
            The ith index of the list holds the difference between the
            activations of the i + 1 layer and the backpropagated error.

        coef_grad : list, length = n_layers - 1
            The ith index contains the amount of change used to update the
            coefficient parameters of the ith layer in an iteration.

        intercept_grads : list, length = n_layers - 1
            The ith index contains the amount of change used to update the
            intercept parameters of the ith layer in an iteration.

        Returns
        -------
        cost : float
        """
        n_samples = X.shape[0]

        # Step (1/3): Forward propagate
        activations = self._forward_pass(activations)

        # Step (2/3): Get cost
        cost = LOSS_FUNCTIONS[self.loss](y, activations[-1])
        # Add L2 regularization term to cost
        values = np.sum(np.array([np.sum(s ** 2) for s in self.layers_coef_]))
        cost += (0.5 * self.alpha) * values / n_samples

        # Step (3/3): Backward propagate
        last = self.n_layers_ - 2

        diff = y - activations[-1]
        deltas[last] = -diff

        # Compute gradient for the last layer
        coef_grads, intercept_grads = self._compute_cost_grad(last, n_samples,
                                                              activations,
                                                              deltas,
                                                              coef_grads,
                                                              intercept_grads)

        # Iterate over the hidden layers
        for i in range(self.n_layers_ - 2, 0, -1):
            deltas[i - 1] = safe_sparse_dot(deltas[i],
                                            self.layers_coef_[i].T)
            derivative = DERIVATIVES[self.activation]
            deltas[i - 1] *= derivative(activations[i])

            coef_grads, \
                intercept_grads = self._compute_cost_grad(i - 1,
                                                          n_samples,
                                                          activations,
                                                          deltas,
                                                          coef_grads,
                                                          intercept_grads)

        return cost, coef_grads, intercept_grads

    def _fit(self, X, y, incremental=False):
        # Make sure self.hidden_layer_sizes is a list
        hidden_layer_sizes = self.hidden_layer_sizes
        if not hasattr(hidden_layer_sizes, "__iter__"):
            hidden_layer_sizes = [hidden_layer_sizes]
        hidden_layer_sizes = list(hidden_layer_sizes)

        # Validate input parameters.
        if np.any(np.array(hidden_layer_sizes) <= 0):
            raise ValueError("hidden_layer_sizes must be > 0, got %s." %
                             hidden_layer_sizes)
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

            if self.classes_ is None or not incremental:
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

        layer_units = ([n_features] + hidden_layer_sizes +
                       [self.n_outputs_])

        # First time training the model
        if self.layers_coef_ is None or (not self.warm_start and
                                         not incremental):
            # Initialize parameters
            self.n_iter_ = 0
            self.t_ = 0
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
                rng = check_random_state(self.random_state)

                n_fan_in = layer_units[i]
                n_fan_out = layer_units[i + 1]

                # Use the initialization method recommended by
                # Glorot et al.
                weight_init_bound = np.sqrt(6. / (n_fan_in + n_fan_out))

                self.layers_coef_.append(rng.uniform(-weight_init_bound,
                                                     weight_init_bound,
                                                     (n_fan_in, n_fan_out)))
                self.layers_intercept_.append(rng.uniform(-weight_init_bound,
                                                          weight_init_bound,
                                                          n_fan_out))

        if self.shuffle:
            X, y = shuffle(X, y, random_state=self.random_state)

        # l-bfgs does not support mini-batches
        if self.algorithm == 'l-bfgs':
            batch_size = n_samples
        else:
            batch_size = np.clip(self.batch_size, 1, n_samples)

        # Initialize lists
        activations = [X]
        activations.extend(np.empty((batch_size, n_fan_out))
                           for n_fan_out in layer_units[1:])
        deltas = [np.empty_like(a_layer) for a_layer in activations]

        coef_grads = [np.empty((n_fan_in, n_fan_out)) for n_fan_in,
                      n_fan_out in zip(layer_units[:-1],
                                       layer_units[1:])]

        intercept_grads = [np.empty(n_fan_out) for n_fan_out in
                           layer_units[1:]]

        # Run the Stochastic Gradient Descent algorithm
        if self.algorithm == 'sgd':
            self._fit_sgd(X, y, activations, deltas, coef_grads,
                          intercept_grads, layer_units, incremental)

        # Run the LBFGS algorithm
        elif self.algorithm == 'l-bfgs':
            self._fit_lbfgs(X, y, activations, deltas, coef_grads,
                            intercept_grads, layer_units)
        return self

    def _fit_lbfgs(self, X, y, activations, deltas, coef_grads, intercept_grads,
                   layer_units):
        # Store meta information for the parameters
        self._coef_indptr = []
        self._intercept_indptr = []
        start = 0

        # Save sizes and indices of coefficients for faster unpacking
        for i in range(self.n_layers_ - 1):
            n_fan_in, n_fan_out = layer_units[i], layer_units[i + 1]

            end = start + (n_fan_in * n_fan_out)
            self._coef_indptr.append((start, end, (n_fan_in, n_fan_out)))
            start = end

        # Save sizes and indices of intercepts for faster unpacking
        for i in range(self.n_layers_ - 1):
            end = start + layer_units[i + 1]
            self._intercept_indptr.append((start, end))
            start = end

        # Run LBFGS
        packed_coef_inter = _pack(self.layers_coef_,
                                  self.layers_intercept_)

        if self.verbose is True or self.verbose >= 1:
            iprint = 1
        else:
            iprint = -1

        optimal_parameters, self.cost_, d = fmin_l_bfgs_b(
            x0=packed_coef_inter,
            func=self._cost_grad_lbfgs,
            maxfun=self.max_iter,
            iprint=iprint,
            pgtol=self.tol,
            args=(X, y, activations, deltas, coef_grads, intercept_grads))

        self._unpack(optimal_parameters)

    def _fit_sgd(self, X, y, activations, deltas, coef_grads, intercept_grads,
                 layer_units, incremental):
        prev_cost = np.inf
        cost_increase_count = 0
        n_samples = X.shape[0]
        batch_size = np.clip(self.batch_size, 1, n_samples)
        intercept_update_prev = [np.zeros_like(grads) for grads in intercept_grads]
        coef_update_prev = [np.zeros_like(grads) for grads in coef_grads]

        try:
            for i in range(self.max_iter):
                for batch_slice in gen_batches(n_samples, batch_size):
                    activations[0] = X[batch_slice]
                    self.cost_, coef_grads, intercept_grads = self._backprop(
                        X[batch_slice], y[batch_slice],
                        activations, deltas, coef_grads,
                        intercept_grads)

                    # update weights
                    for i in range(self.n_layers_ - 1):
                        coef_update_prev[i] = ((1 - self.momentum) * coef_grads[i]
                                               + self.momentum * coef_update_prev[i])
                        self.layers_coef_[i] -= self.learning_rate_ * coef_update_prev[i]

                        intercept_update_prev[i] = ((1 - self.momentum) * intercept_grads[i]
                                                    + self.momentum * intercept_update_prev[i])
                        self.layers_intercept_[i] -= self.learning_rate_ * intercept_update_prev[i]

                    if self.learning_rate == 'invscaling':
                        self.learning_rate_ = self.learning_rate_init / \
                            (self.t_ + 1) ** self.power_t

                if incremental is False:
                    self.n_iter_ += 1

                self.t_ += n_samples

                if self.verbose:
                    print("Iteration %d, cost = %.8f" % (self.n_iter_,
                                                         self.cost_))

                if self.cost_ > prev_cost:
                    cost_increase_count += 1
                    if cost_increase_count == 0.2 * self.max_iter:
                        warnings.warn('Cost is increasing for more than 20%%'
                                      ' of the iterations. Consider reducing'
                                      ' learning_rate_init and preprocessing'
                                      ' your data with StandardScaler or '
                                      ' MinMaxScaler.'
                                      % self.cost_, ConvergenceWarning)

                elif prev_cost - self.cost_ < self.tol or incremental:
                    break

                prev_cost = self.cost_

                if self.n_iter_ == self.max_iter:
                    warnings.warn('SGD: Maximum iterations have reached and'
                                  ' the optimization hasn\'t converged yet.'
                                  % (), ConvergenceWarning)
        except KeyboardInterrupt:
            pass

    def fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            The target values.

        Returns
        -------
        self : returns a trained MLP model.
        """
        return self._fit(X, y, incremental=False)

    def partial_fit(self, X, y):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            The target values.

        Returns
        -------
        self : returns a trained MLP model.
        """
        if self.algorithm != 'sgd':
            raise ValueError("only SGD algorithm supports partial fit")

        return self._fit(X, y, incremental=True)

    def _decision_scores(self, X):
        """Predict using the trained model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y_pred : array-like, shape (n_samples,) or (n_samples, n_outputs)
            The predicted values.
        """
        X = check_array(X, accept_sparse=['csr', 'csc', 'coo'])

        # Make sure self.hidden_layer_sizes is a list
        hidden_layer_sizes = self.hidden_layer_sizes
        if not hasattr(hidden_layer_sizes, "__iter__"):
            hidden_layer_sizes = [hidden_layer_sizes]
        hidden_layer_sizes = list(hidden_layer_sizes)

        layer_units = [X.shape[1]] + hidden_layer_sizes + \
            [self.n_outputs_]

        # Initialize layers
        activations = []
        activations.append(X)

        for i in range(self.n_layers_ - 1):
            activations.append(np.empty((X.shape[0],
                                         layer_units[i + 1])))
        # forward propagate
        self._forward_pass(activations, with_output_activation=False)
        y_pred = activations[-1]

        return y_pred


class MultilayerPerceptronClassifier(BaseMultilayerPerceptron,
                                     ClassifierMixin):
    """Multi-layer Perceptron classifier.

    This algorithm optimizes the logistic loss function using l-bfgs or
    gradient descent.

    Parameters
    ----------
    hidden_layer_sizes : tuple, length = n_layers - 2, default (100,)
        The ith index in list contains the number of neurons in the ith
        hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'relu'
        Activation function for the hidden layer.

        -  'logistic', the logistic sigmoid function,
            returns f(x) = 1 / (1 + exp(x)).

        - 'tanh', the hyperbolic tan function,
           returns f(x) = tanh(x).

        - 'relu', the rectified linear unit function,
           returns f(x) = max(0, x)

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

         Only used when algorithm='sgd'.

    max_iter : int, optional, default 200
        Maximum number of iterations. The algorithm
        iterates until convergence (determined by 'tol') or
        this number of iterations.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    shuffle : bool, optional, default True
        Whether to shuffle samples in each iteration before extracting
        minibatches.

    tol : float, optional, default 1e-5
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached and the algorithm exits.

    learning_rate_init : double, optional, default 0.5
        The initial learning rate used. It controls the step-size
        in updating the weights. Only used when algorithm='sgd'.

    power_t : double, optional, default 0.5
        The exponent for inverse scaling learning rate.
        It is used in updating learning_rate_init when the learning_rate
        is set to 'invscaling'. Only used when algorithm='sgd'.

    verbose : bool, optional, default False
        Whether to print progress messages to stdout.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.

    momentum : float, default 0
        Momentum for gradient descent update.  Should be between 0 and 1. Only
        used when algorithm='sgd'.

    Attributes
    ----------
    `classes_` : array or list of array of shape (n_classes,)
        Class labels for each output.

    `cost_` : float
        The current cost value computed by the loss function.

    `label_binarizer_` : LabelBinarizer
        A LabelBinarizer object trained on the training set.

    `layers_coef_` : list, length n_layers - 1
        The ith element in the list represents the weight matrix corresponding
        to layer i.

    `layers_intercept_` : list, length n_layers - 1
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

    Notes
    -----
    MultilayerPerceptronClassifier trains iteratively since at each time step
    the partial derivatives of the loss function with respect to the model
    parameters are computed to update the parameters.

    It can also use regularizer as a penalty term added to the loss function
    that shrinks model parameters towards zero.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    References
    ----------
    Hinton, Geoffrey E.
        "Connectionist learning procedures." Artificial intelligence 40.1
        (1989): 185-234.

    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """
    def __init__(self, hidden_layer_sizes=(100,), activation="relu",
                 algorithm='l-bfgs', alpha=0.00001,
                 batch_size=200, learning_rate="constant",
                 learning_rate_init=0.5, power_t=0.5, max_iter=200,
                 shuffle=True, random_state=None, tol=1e-5,
                 verbose=False, warm_start=False, momentum=0):

        sup = super(MultilayerPerceptronClassifier, self)
        sup.__init__(hidden_layer_sizes=hidden_layer_sizes,
                     activation=activation, algorithm=algorithm, alpha=alpha,
                     batch_size=batch_size, learning_rate=learning_rate,
                     learning_rate_init=learning_rate_init, power_t=power_t,
                     max_iter=max_iter, loss='log_loss', shuffle=shuffle,
                     random_state=random_state, tol=tol,
                     verbose=verbose, warm_start=warm_start, momentum=momentum)

        self.label_binarizer_ = LabelBinarizer()

    def decision_function(self, X):
        """Decision function of the elm model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : array-like, shape (n_samples,) or (n_samples, n_classes)
            The predicted values.
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
            The input data.

        Returns
        -------
        y : array-like, shape (n_samples,) or (n_samples, n_classes)
            The predicted classes, or the predicted values.
        """
        y_scores = self.decision_function(X)
        y_scores = ACTIVATIONS[self.out_activation_](y_scores)

        return self.label_binarizer_.inverse_transform(y_scores)

    def partial_fit(self, X, y, classes=None):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            The predicted values.

        classes : array, shape (n_classes)
            Classes across all calls to partial_fit.
            Can be obtained by via `np.unique(y_all)`, where y_all is the
            target vector of the entire dataset.
            This argument is required for the first call to partial_fit
            and can be omitted in the subsequent calls.
            Note that y doesn't need to contain all labels in `classes`.

        Returns
        -------
        self : returns a trained MLP model.
        """
        self.classes_ = classes

        super(MultilayerPerceptronClassifier, self).partial_fit(X, y)

        return self

    def predict_log_proba(self, X):
        """Return the log of probability estimates.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            The input data.

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
            The input data.

        Returns
        -------
        y_prob : array-like, shape (n_samples, n_classes)
            The predicted probability of the sample for each class in the
            model, where classes are ordered as they are in `self.classes_`.
        """
        y_scores = self.decision_function(X)

        if y_scores.ndim == 1:
            y_scores = logistic(y_scores)
            return np.vstack([1 - y_scores, y_scores]).T
        else:
            return softmax(y_scores)


class MultilayerPerceptronRegressor(BaseMultilayerPerceptron, RegressorMixin):
    """Multi-layer Perceptron regressor.

    This algorithm optimizes the squared-err loss function using l-bfgs or
    gradient descent.

    Parameters
    ----------
    hidden_layer_sizes : tuple, length = n_layers - 2, default (100,)
        The ith index in list contains the number of neurons in the ith
        hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'relu'
        Activation function for the hidden layer.

        - 'logistic', the logistic sigmoid function,
           returns f(x) = 1 / (1 + exp(x)).

        - 'tanh', the hyperbolic tan function,
           returns f(x) = tanh(x).

        - 'relu', the rectified linear unit function,
           returns f(x) = max(0, x)

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

         Only used when algorithm='sgd'.

    max_iter : int, optional, default 200
        Maximum number of iterations. The algorithm
        iterates until convergence (determined by 'tol') or
        this number of iterations.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    shuffle : bool, optional, default True
        Whether to shuffle samples in each iteration before extracting
        minibatches. Only used when algorithm='sgd'.

    tol : float, optional, default 1e-5
        Tolerance for the optimization. When the loss at iteration i+1 differs
        less than this amount from that at iteration i, convergence is
        considered to be reached and the algorithm exits.

    learning_rate_init : double, optional, default 0.5
        The initial learning rate used. It controls the step-size
        in updating the weights. Only used when algorithm='sgd'.

    power_t : double, optional, default 0.5
        The exponent for inverse scaling learning rate.
        It is used for updating learning_rate_ when it
        is set to 'invscaling'. Only used when algorithm='sgd'.

    verbose : bool, optional, default False
        Whether to print progress messages to stdout.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.

    momentum : float, default 0
        Momentum for gradient descent update.  Should be between 0 and 1. Only
        used when algorithm='sgd'.

    Attributes
    ----------
    `cost_` : float
        The current cost value computed by the loss function.

    `layers_coef_` : list, length n_layers - 1
        The ith element in the list represents the weight matrix corresponding
        to layer i.

    `layers_intercept_` : list, length n_layers - 1
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

    Notes
    -----
    MultilayerPerceptronRegressor trains iteratively since at each time step
    the partial derivatives of the loss function with respect to the model
    parameters are computed to update the parameters.

    It can also use regularizer as a penalty term added to the loss function
    that shrinks model parameters towards zero.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    References
    ----------
    Hinton, Geoffrey E.
        "Connectionist learning procedures." Artificial intelligence 40.1
        (1989): 185-234.

    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """
    def __init__(self, hidden_layer_sizes=(100,), activation="relu",
                 algorithm='l-bfgs', alpha=0.00001,
                 batch_size=200, learning_rate="constant",
                 learning_rate_init=0.1,
                 power_t=0.5, max_iter=100, shuffle=True,
                 random_state=None, tol=1e-5,
                 verbose=False, warm_start=False, momentum=0):

        sup = super(MultilayerPerceptronRegressor, self)
        sup.__init__(hidden_layer_sizes=hidden_layer_sizes,
                     activation=activation, algorithm=algorithm, alpha=alpha,
                     batch_size=batch_size, learning_rate=learning_rate,
                     learning_rate_init=learning_rate_init, power_t=power_t,
                     max_iter=max_iter, loss='squared_loss', shuffle=shuffle,
                     random_state=random_state, tol=tol,
                     verbose=verbose, warm_start=warm_start, momentum=momentum)

    def predict(self, X):
        """Predict using the multi-layer perceptron model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : array-like, shape (n_samples, n_outputs)
            The predicted values.
        """
        return self._decision_scores(X)
