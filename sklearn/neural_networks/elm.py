# -*- coding: utf8
# Author: David C. Lambert <dcl@panix.com>
# License: Simple BSD

###################
# docstring chunks
###################
_ELM_blurb = """An Extreme Learning Machine (ELM) is a single layer feedforward
network with a random hidden layer components and least-squares fitting
of the hidden->output weights by default. [1][2]
"""

_ELM_refs = """References
----------
.. [1] http://www.extreme-learning-machines.org
.. [2] G.-B. Huang, Q.-Y. Zhu and C.-K. Siew, "Extreme Learning Machine:
          Theory and Applications", Neurocomputing, vol. 70, pp. 489-501,
          2006.
"""

############################
# module docstring assembly
############################
__doc__ = """The :mod:`sklearn.neural_networks.elm` module implements the
Extreme Learning Machine Classifiers and Regressors (ELMClassifier,
ELMRegressor, SimpleELMRegressor, SimpleELMClassifier).

{}
{}""".format(_ELM_blurb, _ELM_refs)

from abc import ABCMeta, abstractmethod

import numpy as np
from scipy.linalg import pinv2

from ..utils import as_float_array
from ..utils.extmath import safe_sparse_dot
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..preprocessing import LabelBinarizer

from .random_hidden_layer import SimpleRandomHiddenLayer

__all__ = ["ELMRegressor",
           "ELMClassifier",
           "SimpleELMRegressor",
           "SimpleELMClassifier"]


########################
# more docstring chunks
########################
_ELM_params = """Parameters
----------
`hidden_layer` : random_hidden_layer instance, optional
    (default=SimpleRandomHiddenLayer(random_state=0))

`regressor`    : linear_model instance, optional (default=None)
    If provided, this object is used to perform the regression from hidden
    unit activations to the outputs and subsequent predictions.  If not
    present, a simple least squares fit is performed internally.
"""

_ELM_seealso = """See Also
--------
RBFRandomHiddenLayer, SimpleRandomHiddenLayer, ELMRegressor, ELMClassifier
SimpleELMRegressor, SimpleELMClassifier
"""

_SimpleELM_params = """Parameters
----------
`n_hidden` : int, optional (default=20)
    Number of units to generate in the SimpleRandomHiddenLayer

`xfer_func` : {callable, string} optional (default='tanh')
    Function used to transform input activation
    It must be one of 'tanh', 'sine', 'tribas', 'sigmoid', 'hardlim' or
    a callable.  If none is given, 'tanh' will be used. If a callable
    is given, it will be used to compute the hidden unit activations.

`xfer_args` : dictionary, optional (default=None)
    Supplies keyword arguments for a callable xfer_func

`random_state`  : int, RandomState instance or None (default=None)
    Control the pseudo random number generator used to generate the
    hidden unit weights at fit time.
"""


# little function to propagate docstrings
# tinily tweaked from Paul McGuire
def _take_docstring_from(cls):
    def docstring_decorator(fn):
        fn.__doc__ = getattr(cls, fn.__name__).__doc__
        return fn
    return docstring_decorator


# BaseELM class, regressor and hidden_layer attributes
# and provides defaults for docstrings
class BaseELM(BaseEstimator):
    """Base class for ELMs.

    Warning: This class should not be used directly.
    Use derived classes instead.
    """
    __metaclass__ = ABCMeta

    def __init__(self, hidden_layer, regressor):
        self.regressor = regressor
        self.hidden_layer = hidden_layer

    @abstractmethod
    def fit(self, X, y):
        """
        Fit the model using X, y as training data.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like of shape [n_samples, n_outputs]
            Target values (class labels in classification, real numbers in
            regression)

        Returns
        -------
        self : object

            Returns an instance of self.
        """

    @abstractmethod
    def predict(self, X):
        """Predict values using the model

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape [n_samples, n_features]

        Returns
        -------
        C : numpy array of shape [n_samples, n_outputs]
            Predicted values.
        """


###################################
# ELMRegressor doc string assembly
###################################
_ELMRegressor_doc = """
ELMRegressor is a regressor based on the Extreme Learning Machine.

{}
{}
Attributes
----------
`coefs_` : numpy array
    Fitted regression coefficients if no regressor supplied.

`fitted_` : bool
    Flag set when fit has been called already.

`hidden_activations_` : numpy array of shape [n_samples, n_hidden]
    Hidden layer activations for last input.

{}
{}
""".format(_ELM_blurb, _ELM_params, _ELM_seealso, _ELM_refs)


# workhorse class
class ELMRegressor(BaseELM, RegressorMixin):
    __doc__ = _ELMRegressor_doc

    def __init__(self,
                 hidden_layer=SimpleRandomHiddenLayer(random_state=0),
                 regressor=None):

        super(ELMRegressor, self).__init__(hidden_layer, regressor)

        self.coefs_ = None
        self.fitted_ = False
        self.hidden_activations_ = None

    def _fit_regression(self, y):
        """fit regression using internal least squares/supplied regressor"""
        if (self.regressor is None):
            self.coefs_ = safe_sparse_dot(pinv2(self.hidden_activations_), y)
        else:
            self.regressor.fit(self.hidden_activations_, y)

        self.fitted_ = True

    @_take_docstring_from(BaseELM)
    def fit(self, X, y):
        # fit random hidden layer and compute the hidden layer activations
        self.hidden_activations_ = self.hidden_layer.fit_transform(X)

        # solve the regression from hidden activations to outputs
        self._fit_regression(as_float_array(y, copy=True))

        return self

    def _get_predictions(self, X):
        """get predictions using internal least squares/supplied regressor"""
        if (self.regressor is None):
            preds = safe_sparse_dot(self.hidden_activations_, self.coefs_)
        else:
            preds = self.regressor.predict(self.hidden_activations_)

        return preds

    @_take_docstring_from(BaseELM)
    def predict(self, X):
        if (not self.fitted_):
            raise ValueError("ELMRegressor not fitted")

        # compute hidden layer activations
        self.hidden_activations_ = self.hidden_layer.transform(X)

        # compute output predictions for new hidden activations
        predictions = self._get_predictions(X)

        return predictions


####################################
# ELMClassifier doc string assembly
####################################
_ELMClassifier_doc = """
ELMClassifier is a classifier based on the Extreme Learning Machine.

{}
{}
Attributes
----------
`classes_` : numpy array of shape [n_classes]
    Array of class labels

`binarizer_` : LabelBinarizer instance
    Used to transform class labels

`elm_regressor_` : ELMRegressor instance
    Performs actual fit of binarized values

{}
{}
""".format(_ELM_blurb, _ELM_params, _ELM_seealso, _ELM_refs)


class ELMClassifier(BaseELM, ClassifierMixin):
    __doc__ = _ELMClassifier_doc

    def __init__(self,
                 hidden_layer=SimpleRandomHiddenLayer(random_state=0),
                 regressor=None):

        super(ELMClassifier, self).__init__(hidden_layer, regressor)

        self.classes_ = None
        self.binarizer_ = LabelBinarizer(-1, 1)
        self.elm_regressor_ = ELMRegressor(hidden_layer, regressor)

    def decision_function(self, X):
        """
        This function return the decision function values related to each
        class on an array of test vectors X.

        Parameters
        ----------
        X : array-like of shape [n_samples, n_features]

        Returns
        -------
        C : array of shape [n_samples, n_classes] or [n_samples,]
            Decision function values related to each class, per sample.
            In the two-class case, the shape is [n_samples,]
        """
        return self.elm_regressor_.predict(X)

    @_take_docstring_from(BaseELM)
    def fit(self, X, y):
        self.classes_ = np.unique(y)

        y_bin = self.binarizer_.fit_transform(y)

        self.elm_regressor_.fit(X, y_bin)
        return self

    @_take_docstring_from(BaseELM)
    def predict(self, X):
        raw_predictions = self.decision_function(X)
        class_predictions = self.binarizer_.inverse_transform(raw_predictions)

        return class_predictions


#########################################
# SimpleELMRegressor doc string assembly
#########################################
_SimpleELMRegressor_doc = """
SimpleELMRegressor is a regressor based on the Extreme Learning Machine.

{}
SimpleELMRegressor is a wrapper for an ELMRegressor that uses a
SimpleRandomHiddenLayer and passes the __init__ parameters through
to the hidden layer generated by the fit() method.

{}
Attributes
----------
`elm_regressor_` : ELMRegressor object
    Wrapped object that actually performs the fit.

{}
{}
""".format(_ELM_blurb, _SimpleELM_params, _ELM_seealso, _ELM_refs)


# ELMRegressor with default SimpleRandomHiddenLayer
class SimpleELMRegressor(BaseEstimator, RegressorMixin):
    __doc__ = _SimpleELMRegressor_doc

    def __init__(self, n_hidden=20, xfer_func='tanh', xfer_args=None,
                 random_state=None):

        self.n_hidden = n_hidden
        self.xfer_func = xfer_func
        self.xfer_args = xfer_args
        self.random_state = random_state

        self.elm_regressor_ = None

    @_take_docstring_from(BaseELM)
    def fit(self, X, y):
        rhl = SimpleRandomHiddenLayer(n_hidden=self.n_hidden,
                                      xfer_func=self.xfer_func,
                                      xfer_args=self.xfer_args,
                                      random_state=self.random_state)

        self.elm_regressor_ = ELMRegressor(hidden_layer=rhl)
        self.elm_regressor_.fit(X, y)
        return self

    @_take_docstring_from(BaseELM)
    def predict(self, X):
        if (self.elm_regressor_ is None):
            raise ValueError("SimpleELMRegressor not fitted")

        return self.elm_regressor_.predict(X)


##########################################
# SimpleELMClassifier doc string assembly
##########################################
_SimpleELMClassifier_doc = """
SimpleELMClassifier is a classifier based on the Extreme Learning Machine.

{}
SimpleELMClassifier is a wrapper for an ELMClassifier that uses a
SimpleRandomHiddenLayer and passes the __init__ parameters through
to the hidden layer generated by the fit() method.

{}
Attributes
----------
`classes_` : numpy array of shape [n_classes]
    Array of class labels

`elm_classifier_` : ELMClassifier object
    Wrapped object that actually performs the fit

{}
{}
""".format(_ELM_blurb, _SimpleELM_params, _ELM_seealso, _ELM_refs)


# ELMClassifier with default SimpleRandomHiddenLayer
class SimpleELMClassifier(BaseEstimator, ClassifierMixin):
    __doc__ = _SimpleELMClassifier_doc

    def __init__(self, n_hidden=20, xfer_func='tanh', xfer_args=None,
                 random_state=None):

        self.n_hidden = n_hidden
        self.xfer_func = xfer_func
        self.xfer_args = xfer_args
        self.random_state = random_state

        self.elm_classifier_ = None

    @property
    def classes_(self):
        return self.elm_classifier_.classes_

    @_take_docstring_from(ELMClassifier)
    def decision_function(self, X):
        return self.elm_classifier_.decision_function(X)

    @_take_docstring_from(BaseELM)
    def fit(self, X, y):
        rhl = SimpleRandomHiddenLayer(n_hidden=self.n_hidden,
                                      xfer_func=self.xfer_func,
                                      xfer_args=self.xfer_args,
                                      random_state=self.random_state)

        self.elm_classifier_ = ELMClassifier(hidden_layer=rhl)
        self.elm_classifier_.fit(X, y)

        return self

    @_take_docstring_from(BaseELM)
    def predict(self, X):
        if (self.elm_classifier_ is None):
            raise ValueError("SimpleELMClassifier not fitted")

        return self.elm_classifier_.predict(X)
