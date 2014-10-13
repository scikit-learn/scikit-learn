"""Extreme Learning Machines
"""

# Author: Issam H. Laradji <issam.laradji@gmail.com>
# Licence: BSD 3 clause


from abc import ABCMeta, abstractmethod

import numpy as np
from scipy import linalg

from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from .base import logistic, softmax, ACTIVATIONS
from ..externals import six
from ..preprocessing import LabelBinarizer
from ..metrics import mean_squared_error
from ..linear_model.ridge import ridge_regression
from ..utils import gen_batches, check_random_state
from ..utils import check_array, check_X_y, column_or_1d
from ..utils.extmath import safe_sparse_dot
from ..utils.class_weight import compute_sample_weight


def _multiply_weights(X, sample_weight):
    """Return W*X if sample_weight is not None."""
    if sample_weight is None:
        return X
    else:
        return X * sample_weight[:, np.newaxis]


class BaseELM(six.with_metaclass(ABCMeta, BaseEstimator)):
    """Base class for ELM classification and regression.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """
    @abstractmethod
    def __init__(self, n_hidden, activation, C, class_weight,
                 weight_scale, batch_size, verbose, warm_start,
                 random_state):
        self.C = C
        self.activation = activation
        self.class_weight = class_weight
        self.weight_scale = weight_scale
        self.batch_size = batch_size
        self.n_hidden = n_hidden
        self.verbose = verbose
        self.warm_start = warm_start
        self.random_state = random_state

    def _init_weights(self, n_features):
        """Initialize the parameter weights."""
        rng = check_random_state(self.random_state)

        # Use the initialization method recommended by Glorot et al.
        weight_init_bound = np.sqrt(6. / (n_features + self.n_hidden))

        self.coef_hidden_ = rng.uniform(-weight_init_bound,
                                        weight_init_bound, (n_features,
                                                            self.n_hidden))
        self.intercept_hidden_ = rng.uniform(-weight_init_bound,
                                             weight_init_bound,
                                             self.n_hidden)
        if self.weight_scale != 1:
            self.coef_hidden_ *= self.weight_scale
            self.intercept_hidden_ *= self.weight_scale

    def _compute_hidden_activations(self, X):
        """Compute the hidden activations."""

        hidden_activations = safe_sparse_dot(X, self.coef_hidden_)
        hidden_activations += self.intercept_hidden_

        # Apply the activation method
        activation = ACTIVATIONS[self.activation]
        hidden_activations = activation(hidden_activations)

        return hidden_activations

    def _fit(self, X, y, sample_weight=None, incremental=False):
        """Fit the model to the data X and target y."""
        # Validate input params
        if self.n_hidden <= 0:
            raise ValueError("n_hidden must be > 0, got %s." % self.n_hidden)
        if self.C <= 0.0:
            raise ValueError("C must be > 0, got %s." % self.C)
        if self.activation not in ACTIVATIONS:
            raise ValueError("The activation %s is not supported. Supported "
                             "activation are %s." % (self.activation,
                                                     ACTIVATIONS))

        # Initialize public attributes
        if not hasattr(self, 'classes_'):
            self.classes_ = None
        if not hasattr(self, 'coef_hidden_'):
            self.coef_hidden_ = None

        # Initialize private attributes
        if not hasattr(self, '_HT_H_accumulated'):
            self._HT_H_accumulated = None

        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc', 'coo'],
                         dtype=np.float64, order="C", multi_output=True)

        # This outputs a warning when a 1d array is expected
        if y.ndim == 2 and y.shape[1] == 1:
            y = column_or_1d(y, warn=True)

        # Classification
        if isinstance(self, ClassifierMixin):
            self.label_binarizer_.fit(y)

            if self.classes_ is None or not incremental:
                self.classes_ = self.label_binarizer_.classes_
                if sample_weight is None:
                    sample_weight = compute_sample_weight(self.class_weight,
                                                          self.classes_, y)
            else:
                classes = self.label_binarizer_.classes_
                if not np.all(np.in1d(classes, self.classes_)):
                    raise ValueError("`y` has classes not in `self.classes_`."
                                     " `self.classes_` has %s. 'y' has %s." %
                                     (self.classes_, classes))

            y = self.label_binarizer_.transform(y)

        # Ensure y is 2D
        if y.ndim == 1:
            y = np.reshape(y, (-1, 1))

        n_samples, n_features = X.shape
        self.n_outputs_ = y.shape[1]

        # Step (1/2): Compute the hidden layer coefficients
        if (self.coef_hidden_ is None or (not incremental and
                                          not self.warm_start)):
            # Randomize and scale the input-to-hidden coefficients
            self._init_weights(n_features)

        # Step (2/2): Compute hidden-to-output coefficients
        if self.batch_size is None:
            # Run the least-square algorithm on the whole dataset
            batch_size = n_samples
        else:
            # Run the recursive least-square algorithm on mini-batches
            batch_size = self.batch_size

        batches = gen_batches(n_samples, batch_size)

        # (First time call) Run the least-square algorithm on batch 0
        if not incremental or self._HT_H_accumulated is None:
            batch_slice = next(batches)
            H_batch = self._compute_hidden_activations(X[batch_slice])

            # Get sample weights for the batch
            if sample_weight is None:
                sw = None
            else:
                sw = sample_weight[batch_slice]

            # beta_{0} = inv(H_{0}^T H_{0} + (1. / C) * I) * H_{0}.T y_{0}
            self.coef_output_ = ridge_regression(H_batch, y[batch_slice],
                                                 1. / self.C,
                                                 sample_weight=sw).T

            # Initialize K if this is batch based or partial_fit
            if self.batch_size is not None or incremental:
                # K_{0} = H_{0}^T * W * H_{0}
                weighted_H_batch = _multiply_weights(H_batch, sw)
                self._HT_H_accumulated = safe_sparse_dot(H_batch.T,
                                                         weighted_H_batch)

            if self.verbose:
                y_scores = self._decision_scores(X[batch_slice])

                if self.batch_size is None:
                    verbose_string = "Training mean squared error ="
                else:
                    verbose_string = "Batch 0, Training mean squared error ="

                print("%s %f" % (verbose_string,
                                 mean_squared_error(y[batch_slice], y_scores,
                                                    sample_weight=sw)))

        # Run the least-square algorithm on batch 1, 2, ..., n
        for batch, batch_slice in enumerate(batches):
            # Compute hidden activations H_{i} for batch i
            H_batch = self._compute_hidden_activations(X[batch_slice])

            # Get sample weights (sw) for the batch
            if sample_weight is None:
                sw = None
            else:
                sw = sample_weight[batch_slice]

            weighted_H_batch = _multiply_weights(H_batch, sw)

            # Update K_{i+1} by H_{i}^T * W * H_{i}
            self._HT_H_accumulated += safe_sparse_dot(H_batch.T,
                                                      weighted_H_batch)

            # Update beta_{i+1} by
            # K_{i+1}^{-1} * H_{i+1}^T * W * (y_{i+1} - H_{i+1} * beta_{i})
            y_batch = y[batch_slice] - safe_sparse_dot(H_batch,
                                                       self.coef_output_)

            weighted_y_batch = _multiply_weights(y_batch, sw)
            Hy_batch = safe_sparse_dot(H_batch.T, weighted_y_batch)

            # Update hidden-to-output coefficients
            regularized_HT_H = self._HT_H_accumulated.copy()
            regularized_HT_H.flat[::self.n_hidden + 1] += 1. / self.C

            # It is safe to use linalg.solve (instead of linalg.lstsq
            # which is slow) since it is highly unlikely that
            # regularized_HT_H is singular due to the random
            # projection of the first layer and 'C' regularization being
            # not dangerously large.
            self.coef_output_ += linalg.solve(regularized_HT_H, Hy_batch,
                                              sym_pos=True, overwrite_a=True,
                                              overwrite_b=True)
            if self.verbose:
                y_scores = self._decision_scores(X[batch_slice])
                print("Batch %d, Training mean squared error = %f" %
                      (batch + 1, mean_squared_error(y[batch_slice], y_scores,
                                                     sample_weight=sw)))
        return self

    def fit(self, X, y, sample_weight=None):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            Target values.

        sample_weight : array-like, shape (n_samples,)
            Per-sample weights. Rescale C per sample. Higher weights
            force the classifier to put more emphasis on these points.

        Returns
        -------
        self : returns a trained ELM usable for prediction.
        """
        return self._fit(X, y, sample_weight=sample_weight, incremental=False)

    def partial_fit(self, X, y, sample_weight=None):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            Subset of training data.

        y : array-like, shape (n_samples,)
            Subset of target values.

        sample_weight : array-like, shape (n_samples,)
            Per-sample weights. Rescale C per sample. Higher weights
            force the classifier to put more emphasis on these points.

        Returns
        -------
        self : returns a trained ELM usable for prediction.
        """
        self._fit(X, y, sample_weight=sample_weight, incremental=True)

        return self

    def _decision_scores(self, X):
        """Predict using the ELM model

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

        if self.batch_size is None:
            hidden_activations = self._compute_hidden_activations(X)
            y_pred = safe_sparse_dot(hidden_activations, self.coef_output_)
        else:
            n_samples = X.shape[0]
            batches = gen_batches(n_samples, self.batch_size)

            y_pred = np.zeros((n_samples, self.n_outputs_))
            for batch in batches:
                h_batch = self._compute_hidden_activations(X[batch])
                y_pred[batch] = safe_sparse_dot(h_batch, self.coef_output_)

        return y_pred


class ELMClassifier(BaseELM, ClassifierMixin):
    """Extreme learning machine classifier.

    The algorithm trains a single-hidden layer feedforward network by computing
    the hidden layer values using randomized parameters, then solving
    for the output weights using least-square solutions. For prediction,
    after computing the forward pass, the continuous output values pass
    through a gate function converting them to integers that represent classes.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    Parameters
    ----------
    C : float, optional, default 100
        A regularization term that controls the linearity of the decision
        function. Smaller value of C makes the decision boundary more linear.

    class_weight : {dict, 'auto', None}, default None
        If 'auto', class weights will be given inversely proportional
        to the frequency of the class in the data.
        If a dictionary is given, keys are the class labels and the
        corresponding values are the class weights.
        If None is given, then no class weights will be applied.

    weight_scale : float, default 1.
        Initializes and scales the input-to-hidden weights.
        The weight values will range between plus and minus
        'sqrt(weight_scale * 6. / (n_features + n_hidden))' based on the
        uniform distribution.

    n_hidden : int, default 100
        The number of units in the hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'relu'
        Activation function for the hidden layer.

         - 'logistic', the logistic sigmoid function,
            returns f(x) = 1 / (1 + exp(x)).

         - 'tanh', the hyperbolic tan function,
            returns f(x) = tanh(x).

         - 'relu', the rectified linear unit function,
            returns f(x) = max(0, x).

    batch_size : int, optional, default None
        If None is given, batch_size is set as the number of samples.
        Otherwise, it will be set as the given integer.

    verbose : bool, optional, default False
        Whether to print the training score.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    Attributes
    ----------
    `classes_` : array-list, shape (n_classes,)
        Class labels for each output.

    `n_outputs_` : int
        Number of output neurons.

    `coef_hidden_` : array-like, shape (n_features, n_hidden)
        The input-to-hidden weights.

    `intercept_hidden_` : array-like, shape (n_hidden,)
        The bias added to the hidden layer neurons.

    `coef_output_` : array-like, shape (n_hidden, n_outputs_)
        The hidden-to-output weights.

    `label_binarizer_` : LabelBinarizer
        A LabelBinarizer object trained on the training set.

    References
    ----------
    Liang, Nan-Ying, et al.
        "A fast and accurate online sequential learning algorithm for
        feedforward networks." Neural Networks, IEEE Transactions on
        17.6 (2006): 1411-1423.
        http://www.ntu.edu.sg/home/egbhuang/pdf/OS-ELM-TNN.pdf

    Zong, Weiwei, Guang-Bin Huang, and Yiqiang Chen.
        "Weighted extreme learning machine for imbalance learning."
        Neurocomputing 101 (2013): 229-242.

    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """
    def __init__(self, n_hidden=100, activation='relu', C=1,
                 class_weight=None, weight_scale=1.0, batch_size=None,
                 verbose=False, warm_start=False, random_state=None):
        super(ELMClassifier, self).__init__(n_hidden=n_hidden,
                                            activation=activation,
                                            C=C, class_weight=class_weight,
                                            weight_scale=weight_scale,
                                            batch_size=batch_size,
                                            verbose=verbose,
                                            warm_start=warm_start,
                                            random_state=random_state)

        self.label_binarizer_ = LabelBinarizer(-1, 1)

    def partial_fit(self, X, y, classes=None, sample_weight=None):
        """Fit the model to the data X and target y.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        y : array-like, shape (n_samples,)
            Subset of the target values.

        classes : array-like, shape (n_classes,)
            List of all the classes that can possibly appear in the y vector.

            Must be provided at the first call to partial_fit, can be omitted
            in subsequent calls.

        sample_weight : array-like, shape (n_samples,)
            Per-sample weights. Rescale C per sample. Higher weights
            force the classifier to put more emphasis on these points.

        Returns
        -------
        self : returns a trained elm usable for prediction.
        """
        self.classes_ = classes

        super(ELMClassifier, self).partial_fit(X, y, sample_weight)

        return self

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
        """Predict using the ELM model

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : array-like, shape (n_samples,) or (n_samples, n_classes)
            The predicted classes, or the predicted values.
        """
        y_scores = self._decision_scores(X)

        return self.label_binarizer_.inverse_transform(y_scores)

    def predict_proba(self, X):
        """Probability estimates.

        Warning: the estimates aren't callibrated since the model optimizes a
        penalized least squares objective function based on the One Vs Rest
        binary encoding of the class membership.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y_prob : array-like, shape (n_samples, n_classes)
            The predicted probability of the sample for each class in the
            model, where classes are ordered as they are in
            `self.classes_`.
        """
        y_scores = self._decision_scores(X)

        if len(self.classes_) == 2:
            y_scores = logistic(y_scores)
            return np.hstack([1 - y_scores, y_scores])
        else:
            return softmax(y_scores)


class ELMRegressor(BaseELM, RegressorMixin):
    """Extreme learning machine regressor.

    The algorithm trains a single-hidden layer feedforward network by computing
    the hidden layer values using randomized parameters, then solving
    for the output weights using least-square solutions. For prediction,
    ELMRegressor computes the forward pass resulting in continuous output
    values.

    This implementation works with data represented as dense and sparse numpy
    arrays of floating point values for the features.

    Parameters
    ----------
    C : float, optional, default 100
        A regularization term that controls the linearity of the decision
        function. Smaller value of C makes the decision boundary more linear.

    weight_scale : float, default 1.
        Initializes and scales the input-to-hidden weights.
        The weight values will range between plus and minus
        'sqrt(weight_scale * 6. / (n_features + n_hidden))' based on the
        uniform distribution.

    n_hidden : int, default 100
        The number of units in the hidden layer.

    activation : {'logistic', 'tanh', 'relu'}, default 'relu'
        Activation function for the hidden layer.

         - 'logistic', the logistic sigmoid function,
            returns f(x) = 1 / (1 + exp(x)).

         - 'tanh', the hyperbolic tan function,
            returns f(x) = tanh(x).

         - 'relu', the rectified linear unit function,
            returns f(x) = max(0, x).

    batch_size : int, optional, default None
        If None is given, batch_size is set as the number of samples.
        Otherwise, it will be set as the given integer.

    verbose : bool, optional, default False
        Whether to print the training score.

    warm_start : bool, optional, default False
        When set to True, reuse the solution of the previous
        call to fit as initialization, otherwise, just erase the
        previous solution.

    random_state : int or RandomState, optional, default None
        State of or seed for random number generator.

    Attributes
    ----------
    `classes_` : array-list, shape (n_classes,)
        Class labels for each output.

    `n_outputs_` : int
        Number of output neurons.

    `coef_hidden_` : array-like, shape (n_features, n_hidden)
        The input-to-hidden weights.

    `intercept_hidden_` : array-like, shape (n_hidden,)
        The bias added to the hidden layer neurons.

    `coef_output_` : array-like, shape (n_hidden, n_outputs_)
        The hidden-to-output weights.

    References
    ----------
    Liang, Nan-Ying, et al.
        "A fast and accurate online sequential learning algorithm for
        feedforward networks." Neural Networks, IEEE Transactions on
        17.6 (2006): 1411-1423.
        http://www.ntu.edu.sg/home/egbhuang/pdf/OS-ELM-TNN.pdf

    Zong, Weiwei, Guang-Bin Huang, and Yiqiang Chen.
        "Weighted extreme learning machine for imbalance learning."
        Neurocomputing 101 (2013): 229-242.

    Glorot, Xavier, and Yoshua Bengio. "Understanding the difficulty of
        training deep feedforward neural networks." International Conference
        on Artificial Intelligence and Statistics. 2010.
    """
    def __init__(self, n_hidden=100, activation='relu', weight_scale=1.0,
                 batch_size=None, C=1, verbose=False, warm_start=False,
                 random_state=None):
        super(ELMRegressor, self).__init__(n_hidden=n_hidden,
                                           activation=activation,
                                           C=C, class_weight=None,
                                           weight_scale=weight_scale,
                                           batch_size=batch_size,
                                           verbose=verbose,
                                           warm_start=warm_start,
                                           random_state=random_state)

    def predict(self, X):
        """Predict using the ELM model.

        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The input data.

        Returns
        -------
        y : array-like, shape (n_samples,) or (n_samples, n_outputs)
            The predicted values.
        """
        y_pred = self._decision_scores(X)

        if self.n_outputs_ == 1:
            return y_pred.ravel()
        else:
            return y_pred
