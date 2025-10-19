"""Base classes for all estimators and various utility functions."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import copy
import functools
import inspect
import platform
import re
import warnings
from collections import defaultdict

import numpy as np

from sklearn import __version__
from sklearn._config import config_context, get_config
from sklearn.exceptions import InconsistentVersionWarning
from sklearn.utils._metadata_requests import _MetadataRequester, _routing_enabled
from sklearn.utils._missing import is_pandas_na, is_scalar_nan
from sklearn.utils._param_validation import validate_parameter_constraints
from sklearn.utils._repr_html.base import ReprHTMLMixin, _HTMLDocumentationLinkMixin
from sklearn.utils._repr_html.estimator import estimator_html_repr
from sklearn.utils._repr_html.params import ParamsDict
from sklearn.utils._set_output import _SetOutputMixin
from sklearn.utils._tags import (
    ClassifierTags,
    RegressorTags,
    Tags,
    TargetTags,
    TransformerTags,
    get_tags,
)
from sklearn.utils.fixes import _IS_32BIT
from sklearn.utils.validation import (
    _check_feature_names_in,
    _generate_get_feature_names_out,
    _is_fitted,
    check_array,
    check_is_fitted,
)


def clone(estimator, *, safe=True):
    """Construct a new unfitted estimator with the same parameters.

    Clone does a deep copy of the model in an estimator
    without actually copying attached data. It returns a new estimator
    with the same parameters that has not been fitted on any data.

    .. versionchanged:: 1.3
        Delegates to `estimator.__sklearn_clone__` if the method exists.

    Parameters
    ----------
    estimator : {list, tuple, set} of estimator instance or a single \
            estimator instance
        The estimator or group of estimators to be cloned.
    safe : bool, default=True
        If safe is False, clone will fall back to a deep copy on objects
        that are not estimators. Ignored if `estimator.__sklearn_clone__`
        exists.

    Returns
    -------
    estimator : object
        The deep copy of the input, an estimator if input is an estimator.

    Notes
    -----
    If the estimator's `random_state` parameter is an integer (or if the
    estimator doesn't have a `random_state` parameter), an *exact clone* is
    returned: the clone and the original estimator will give the exact same
    results. Otherwise, *statistical clone* is returned: the clone might
    return different results from the original estimator. More details can be
    found in :ref:`randomness`.

    Examples
    --------
    >>> from sklearn.base import clone
    >>> from sklearn.linear_model import LogisticRegression
    >>> X = [[-1, 0], [0, 1], [0, -1], [1, 0]]
    >>> y = [0, 0, 1, 1]
    >>> classifier = LogisticRegression().fit(X, y)
    >>> cloned_classifier = clone(classifier)
    >>> hasattr(classifier, "classes_")
    True
    >>> hasattr(cloned_classifier, "classes_")
    False
    >>> classifier is cloned_classifier
    False
    """
    if hasattr(estimator, "__sklearn_clone__") and not inspect.isclass(estimator):
        return estimator.__sklearn_clone__()
    return _clone_parametrized(estimator, safe=safe)


def _clone_parametrized(estimator, *, safe=True):
    """Default implementation of clone. See :func:`sklearn.base.clone` for details."""

    estimator_type = type(estimator)
    if estimator_type is dict:
        return {k: clone(v, safe=safe) for k, v in estimator.items()}
    elif estimator_type in (list, tuple, set, frozenset):
        return estimator_type([clone(e, safe=safe) for e in estimator])
    elif not hasattr(estimator, "get_params") or isinstance(estimator, type):
        if not safe:
            return copy.deepcopy(estimator)
        else:
            if isinstance(estimator, type):
                raise TypeError(
                    "Cannot clone object. "
                    "You should provide an instance of "
                    "scikit-learn estimator instead of a class."
                )
            else:
                raise TypeError(
                    "Cannot clone object '%s' (type %s): "
                    "it does not seem to be a scikit-learn "
                    "estimator as it does not implement a "
                    "'get_params' method." % (repr(estimator), type(estimator))
                )

    klass = estimator.__class__
    new_object_params = estimator.get_params(deep=False)
    for name, param in new_object_params.items():
        new_object_params[name] = clone(param, safe=False)

    new_object = klass(**new_object_params)
    try:
        new_object._metadata_request = copy.deepcopy(estimator._metadata_request)
    except AttributeError:
        pass

    params_set = new_object.get_params(deep=False)

    # quick sanity check of the parameters of the clone
    for name in new_object_params:
        param1 = new_object_params[name]
        param2 = params_set[name]
        if param1 is not param2:
            raise RuntimeError(
                "Cannot clone object %s, as the constructor "
                "either does not set or modifies parameter %s" % (estimator, name)
            )

    # _sklearn_output_config is used by `set_output` to configure the output
    # container of an estimator.
    if hasattr(estimator, "_sklearn_output_config"):
        new_object._sklearn_output_config = copy.deepcopy(
            estimator._sklearn_output_config
        )
    return new_object


class BaseEstimator(ReprHTMLMixin, _HTMLDocumentationLinkMixin, _MetadataRequester):
    """Base class for all estimators in scikit-learn.

    Inheriting from this class provides default implementations of:

    - setting and getting parameters used by `GridSearchCV` and friends;
    - textual and HTML representation displayed in terminals and IDEs;
    - estimator serialization;
    - parameters validation;
    - data validation;
    - feature names validation.

    Read more in the :ref:`User Guide <rolling_your_own_estimator>`.


    Notes
    -----
    All estimators should specify all the parameters that can be set
    at the class level in their ``__init__`` as explicit keyword
    arguments (no ``*args`` or ``**kwargs``).

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator
    >>> class MyEstimator(BaseEstimator):
    ...     def __init__(self, *, param=1):
    ...         self.param = param
    ...     def fit(self, X, y=None):
    ...         self.is_fitted_ = True
    ...         return self
    ...     def predict(self, X):
    ...         return np.full(shape=X.shape[0], fill_value=self.param)
    >>> estimator = MyEstimator(param=2)
    >>> estimator.get_params()
    {'param': 2}
    >>> X = np.array([[1, 2], [2, 3], [3, 4]])
    >>> y = np.array([1, 0, 1])
    >>> estimator.fit(X, y).predict(X)
    array([2, 2, 2])
    >>> estimator.set_params(param=3).fit(X, y).predict(X)
    array([3, 3, 3])
    """

    def __dir__(self):
        # Filters conditional methods that should be hidden based
        # on the `available_if` decorator
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=FutureWarning)
            return [attr for attr in super().__dir__() if hasattr(self, attr)]

    _html_repr = estimator_html_repr

    @classmethod
    def _get_param_names(cls):
        """Get parameter names for the estimator"""
        # fetch the constructor or the original constructor before
        # deprecation wrapping if any
        init = getattr(cls.__init__, "deprecated_original", cls.__init__)
        if init is object.__init__:
            # No explicit constructor to introspect
            return []

        # introspect the constructor arguments to find the model parameters
        # to represent
        init_signature = inspect.signature(init)
        # Consider the constructor parameters excluding 'self'
        parameters = [
            p
            for p in init_signature.parameters.values()
            if p.name != "self" and p.kind != p.VAR_KEYWORD
        ]
        for p in parameters:
            if p.kind == p.VAR_POSITIONAL:
                raise RuntimeError(
                    "scikit-learn estimators should always "
                    "specify their parameters in the signature"
                    " of their __init__ (no varargs)."
                    " %s with constructor %s doesn't "
                    " follow this convention." % (cls, init_signature)
                )
        # Extract and sort argument names excluding 'self'
        return sorted([p.name for p in parameters])

    def get_params(self, deep=True):
        """
        Get parameters for this estimator.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : dict
            Parameter names mapped to their values.
        """
        out = dict()
        for key in self._get_param_names():
            value = getattr(self, key)
            if deep and hasattr(value, "get_params") and not isinstance(value, type):
                deep_items = value.get_params().items()
                out.update((key + "__" + k, val) for k, val in deep_items)
            out[key] = value
        return out

    def _get_params_html(self, deep=True, doc_link=""):
        """
        Get parameters for this estimator with a specific HTML representation.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        doc_link : str
            URL to the estimator documentation.
            Used for linking to the estimator's parameters documentation
            available in HTML displays.

        Returns
        -------
        params : ParamsDict
            Parameter names mapped to their values. We return a `ParamsDict`
            dictionary, which renders a specific HTML representation in table
            form.
        """
        out = self.get_params(deep=deep)

        init_func = getattr(self.__init__, "deprecated_original", self.__init__)
        init_default_params = inspect.signature(init_func).parameters
        init_default_params = {
            name: param.default for name, param in init_default_params.items()
        }

        def is_non_default(param_name, param_value):
            """Finds the parameters that have been set by the user."""
            if param_name not in init_default_params:
                # happens if k is part of a **kwargs
                return True
            if init_default_params[param_name] == inspect._empty:
                # k has no default value
                return True
            # avoid calling repr on nested estimators
            if isinstance(param_value, BaseEstimator) and type(param_value) is not type(
                init_default_params[param_name]
            ):
                return True
            if is_pandas_na(param_value) and not is_pandas_na(
                init_default_params[param_name]
            ):
                return True
            if not np.array_equal(
                param_value, init_default_params[param_name]
            ) and not (
                is_scalar_nan(init_default_params[param_name])
                and is_scalar_nan(param_value)
            ):
                return True

            return False

        # reorder the parameters from `self.get_params` using the `__init__`
        # signature
        remaining_params = [name for name in out if name not in init_default_params]
        ordered_out = {name: out[name] for name in init_default_params if name in out}
        ordered_out.update({name: out[name] for name in remaining_params})

        non_default_ls = tuple(
            [name for name, value in ordered_out.items() if is_non_default(name, value)]
        )

        return ParamsDict(
            params=ordered_out,
            non_default=non_default_ls,
            estimator_class=self.__class__,
            doc_link=doc_link,
        )

    def set_params(self, **params):
        """Set the parameters of this estimator.

        The method works on simple estimators as well as on nested objects
        (such as :class:`~sklearn.pipeline.Pipeline`). The latter have
        parameters of the form ``<component>__<parameter>`` so that it's
        possible to update each component of a nested object.

        Parameters
        ----------
        **params : dict
            Estimator parameters.

        Returns
        -------
        self : estimator instance
            Estimator instance.
        """
        if not params:
            # Simple optimization to gain speed (inspect is slow)
            return self
        valid_params = self.get_params(deep=True)

        nested_params = defaultdict(dict)  # grouped by prefix
        for key, value in params.items():
            key, delim, sub_key = key.partition("__")
            if key not in valid_params:
                local_valid_params = self._get_param_names()
                raise ValueError(
                    f"Invalid parameter {key!r} for estimator {self}. "
                    f"Valid parameters are: {local_valid_params!r}."
                )

            if delim:
                nested_params[key][sub_key] = value
            else:
                setattr(self, key, value)
                valid_params[key] = value

        for key, sub_params in nested_params.items():
            valid_params[key].set_params(**sub_params)

        return self

    def __sklearn_clone__(self):
        return _clone_parametrized(self)

    def __repr__(self, N_CHAR_MAX=700):
        # N_CHAR_MAX is the (approximate) maximum number of non-blank
        # characters to render. We pass it as an optional parameter to ease
        # the tests.

        from sklearn.utils._pprint import _EstimatorPrettyPrinter

        N_MAX_ELEMENTS_TO_SHOW = 30  # number of elements to show in sequences

        # use ellipsis for sequences with a lot of elements
        pp = _EstimatorPrettyPrinter(
            compact=True,
            indent=1,
            indent_at_name=True,
            n_max_elements_to_show=N_MAX_ELEMENTS_TO_SHOW,
        )

        repr_ = pp.pformat(self)

        # Use bruteforce ellipsis when there are a lot of non-blank characters
        n_nonblank = len("".join(repr_.split()))
        if n_nonblank > N_CHAR_MAX:
            lim = N_CHAR_MAX // 2  # apprx number of chars to keep on both ends
            regex = r"^(\s*\S){%d}" % lim
            # The regex '^(\s*\S){%d}' % n
            # matches from the start of the string until the nth non-blank
            # character:
            # - ^ matches the start of string
            # - (pattern){n} matches n repetitions of pattern
            # - \s*\S matches a non-blank char following zero or more blanks
            left_lim = re.match(regex, repr_).end()
            right_lim = re.match(regex, repr_[::-1]).end()

            if "\n" in repr_[left_lim:-right_lim]:
                # The left side and right side aren't on the same line.
                # To avoid weird cuts, e.g.:
                # categoric...ore',
                # we need to start the right side with an appropriate newline
                # character so that it renders properly as:
                # categoric...
                # handle_unknown='ignore',
                # so we add [^\n]*\n which matches until the next \n
                regex += r"[^\n]*\n"
                right_lim = re.match(regex, repr_[::-1]).end()

            ellipsis = "..."
            if left_lim + len(ellipsis) < len(repr_) - right_lim:
                # Only add ellipsis if it results in a shorter repr
                repr_ = repr_[:left_lim] + "..." + repr_[-right_lim:]

        return repr_

    def __getstate__(self):
        if getattr(self, "__slots__", None):
            raise TypeError(
                "You cannot use `__slots__` in objects inheriting from "
                "`sklearn.base.BaseEstimator`."
            )

        try:
            state = super().__getstate__()
            if state is None:
                # For Python 3.11+, empty instance (no `__slots__`,
                # and `__dict__`) will return a state equal to `None`.
                state = self.__dict__.copy()
        except AttributeError:
            # Python < 3.11
            state = self.__dict__.copy()

        if type(self).__module__.startswith("sklearn."):
            return dict(state.items(), _sklearn_version=__version__)
        else:
            return state

    def __setstate__(self, state):
        if type(self).__module__.startswith("sklearn."):
            pickle_version = state.pop("_sklearn_version", "pre-0.18")
            if pickle_version != __version__:
                warnings.warn(
                    InconsistentVersionWarning(
                        estimator_name=self.__class__.__name__,
                        current_sklearn_version=__version__,
                        original_sklearn_version=pickle_version,
                    ),
                )
        try:
            super().__setstate__(state)
        except AttributeError:
            self.__dict__.update(state)

    def __sklearn_tags__(self):
        return Tags(
            estimator_type=None,
            target_tags=TargetTags(required=False),
            transformer_tags=None,
            regressor_tags=None,
            classifier_tags=None,
        )

    def _validate_params(self):
        """Validate types and values of constructor parameters

        The expected type and values must be defined in the `_parameter_constraints`
        class attribute, which is a dictionary `param_name: list of constraints`. See
        the docstring of `validate_parameter_constraints` for a description of the
        accepted constraints.
        """
        validate_parameter_constraints(
            self._parameter_constraints,
            self.get_params(deep=False),
            caller_name=self.__class__.__name__,
        )


class ClassifierMixin:
    """Mixin class for all classifiers in scikit-learn.

    This mixin defines the following functionality:

    - set estimator type to `"classifier"` through the `estimator_type` tag;
    - `score` method that default to :func:`~sklearn.metrics.accuracy_score`.
    - enforce that `fit` requires `y` to be passed through the `requires_y` tag,
      which is done by setting the classifier type tag.

    Read more in the :ref:`User Guide <rolling_your_own_estimator>`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator, ClassifierMixin
    >>> # Mixin classes should always be on the left-hand side for a correct MRO
    >>> class MyEstimator(ClassifierMixin, BaseEstimator):
    ...     def __init__(self, *, param=1):
    ...         self.param = param
    ...     def fit(self, X, y=None):
    ...         self.is_fitted_ = True
    ...         return self
    ...     def predict(self, X):
    ...         return np.full(shape=X.shape[0], fill_value=self.param)
    >>> estimator = MyEstimator(param=1)
    >>> X = np.array([[1, 2], [2, 3], [3, 4]])
    >>> y = np.array([1, 0, 1])
    >>> estimator.fit(X, y).predict(X)
    array([1, 1, 1])
    >>> estimator.score(X, y)
    0.66...
    """

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.estimator_type = "classifier"
        tags.classifier_tags = ClassifierTags()
        tags.target_tags.required = True
        return tags

    def score(self, X, y, sample_weight=None):
        """
        Return :ref:`accuracy <accuracy_score>` on provided data and labels.

        In multi-label classification, this is the subset accuracy
        which is a harsh metric since you require for each sample that
        each label set be correctly predicted.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test samples.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            True labels for `X`.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        score : float
            Mean accuracy of ``self.predict(X)`` w.r.t. `y`.
        """
        from sklearn.metrics import accuracy_score

        return accuracy_score(y, self.predict(X), sample_weight=sample_weight)


class RegressorMixin:
    """Mixin class for all regression estimators in scikit-learn.

    This mixin defines the following functionality:

    - set estimator type to `"regressor"` through the `estimator_type` tag;
    - `score` method that default to :func:`~sklearn.metrics.r2_score`.
    - enforce that `fit` requires `y` to be passed through the `requires_y` tag,
      which is done by setting the regressor type tag.

    Read more in the :ref:`User Guide <rolling_your_own_estimator>`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator, RegressorMixin
    >>> # Mixin classes should always be on the left-hand side for a correct MRO
    >>> class MyEstimator(RegressorMixin, BaseEstimator):
    ...     def __init__(self, *, param=1):
    ...         self.param = param
    ...     def fit(self, X, y=None):
    ...         self.is_fitted_ = True
    ...         return self
    ...     def predict(self, X):
    ...         return np.full(shape=X.shape[0], fill_value=self.param)
    >>> estimator = MyEstimator(param=0)
    >>> X = np.array([[1, 2], [2, 3], [3, 4]])
    >>> y = np.array([-1, 0, 1])
    >>> estimator.fit(X, y).predict(X)
    array([0, 0, 0])
    >>> estimator.score(X, y)
    0.0
    """

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.estimator_type = "regressor"
        tags.regressor_tags = RegressorTags()
        tags.target_tags.required = True
        return tags

    def score(self, X, y, sample_weight=None):
        """Return :ref:`coefficient of determination <r2_score>` on test data.

        The coefficient of determination, :math:`R^2`, is defined as
        :math:`(1 - \\frac{u}{v})`, where :math:`u` is the residual
        sum of squares ``((y_true - y_pred)** 2).sum()`` and :math:`v`
        is the total sum of squares ``((y_true - y_true.mean()) ** 2).sum()``.
        The best possible score is 1.0 and it can be negative (because the
        model can be arbitrarily worse). A constant model that always predicts
        the expected value of `y`, disregarding the input features, would get
        a :math:`R^2` score of 0.0.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test samples. For some estimators this may be a precomputed
            kernel matrix or a list of generic objects instead with shape
            ``(n_samples, n_samples_fitted)``, where ``n_samples_fitted``
            is the number of samples used in the fitting for the estimator.

        y : array-like of shape (n_samples,) or (n_samples, n_outputs)
            True values for `X`.

        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.

        Returns
        -------
        score : float
            :math:`R^2` of ``self.predict(X)`` w.r.t. `y`.

        Notes
        -----
        The :math:`R^2` score used when calling ``score`` on a regressor uses
        ``multioutput='uniform_average'`` from version 0.23 to keep consistent
        with default value of :func:`~sklearn.metrics.r2_score`.
        This influences the ``score`` method of all the multioutput
        regressors (except for
        :class:`~sklearn.multioutput.MultiOutputRegressor`).
        """

        from sklearn.metrics import r2_score

        y_pred = self.predict(X)
        return r2_score(y, y_pred, sample_weight=sample_weight)


class ClusterMixin:
    """Mixin class for all cluster estimators in scikit-learn.

    - set estimator type to `"clusterer"` through the `estimator_type` tag;
    - `fit_predict` method returning the cluster labels associated to each sample.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator, ClusterMixin
    >>> class MyClusterer(ClusterMixin, BaseEstimator):
    ...     def fit(self, X, y=None):
    ...         self.labels_ = np.ones(shape=(len(X),), dtype=np.int64)
    ...         return self
    >>> X = [[1, 2], [2, 3], [3, 4]]
    >>> MyClusterer().fit_predict(X)
    array([1, 1, 1])
    """

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.estimator_type = "clusterer"
        if tags.transformer_tags is not None:
            tags.transformer_tags.preserves_dtype = []
        return tags

    def fit_predict(self, X, y=None, **kwargs):
        """
        Perform clustering on `X` and returns cluster labels.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        y : Ignored
            Not used, present for API consistency by convention.

        **kwargs : dict
            Arguments to be passed to ``fit``.

            .. versionadded:: 1.4

        Returns
        -------
        labels : ndarray of shape (n_samples,), dtype=np.int64
            Cluster labels.
        """
        # non-optimized default implementation; override when a better
        # method is possible for a given clustering algorithm
        self.fit(X, **kwargs)
        return self.labels_


class BiclusterMixin:
    """Mixin class for all bicluster estimators in scikit-learn.

    This mixin defines the following functionality:

    - `biclusters_` property that returns the row and column indicators;
    - `get_indices` method that returns the row and column indices of a bicluster;
    - `get_shape` method that returns the shape of a bicluster;
    - `get_submatrix` method that returns the submatrix corresponding to a bicluster.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator, BiclusterMixin
    >>> class DummyBiClustering(BiclusterMixin, BaseEstimator):
    ...     def fit(self, X, y=None):
    ...         self.rows_ = np.ones(shape=(1, X.shape[0]), dtype=bool)
    ...         self.columns_ = np.ones(shape=(1, X.shape[1]), dtype=bool)
    ...         return self
    >>> X = np.array([[1, 1], [2, 1], [1, 0],
    ...               [4, 7], [3, 5], [3, 6]])
    >>> bicluster = DummyBiClustering().fit(X)
    >>> hasattr(bicluster, "biclusters_")
    True
    >>> bicluster.get_indices(0)
    (array([0, 1, 2, 3, 4, 5]), array([0, 1]))
    """

    @property
    def biclusters_(self):
        """Convenient way to get row and column indicators together.

        Returns the ``rows_`` and ``columns_`` members.
        """
        return self.rows_, self.columns_

    def get_indices(self, i):
        """Row and column indices of the `i`'th bicluster.

        Only works if ``rows_`` and ``columns_`` attributes exist.

        Parameters
        ----------
        i : int
            The index of the cluster.

        Returns
        -------
        row_ind : ndarray, dtype=np.intp
            Indices of rows in the dataset that belong to the bicluster.
        col_ind : ndarray, dtype=np.intp
            Indices of columns in the dataset that belong to the bicluster.
        """
        rows = self.rows_[i]
        columns = self.columns_[i]
        return np.nonzero(rows)[0], np.nonzero(columns)[0]

    def get_shape(self, i):
        """Shape of the `i`'th bicluster.

        Parameters
        ----------
        i : int
            The index of the cluster.

        Returns
        -------
        n_rows : int
            Number of rows in the bicluster.

        n_cols : int
            Number of columns in the bicluster.
        """
        indices = self.get_indices(i)
        return tuple(len(i) for i in indices)

    def get_submatrix(self, i, data):
        """Return the submatrix corresponding to bicluster `i`.

        Parameters
        ----------
        i : int
            The index of the cluster.
        data : array-like of shape (n_samples, n_features)
            The data.

        Returns
        -------
        submatrix : ndarray of shape (n_rows, n_cols)
            The submatrix corresponding to bicluster `i`.

        Notes
        -----
        Works with sparse matrices. Only works if ``rows_`` and
        ``columns_`` attributes exist.
        """

        data = check_array(data, accept_sparse="csr")
        row_ind, col_ind = self.get_indices(i)
        return data[row_ind[:, np.newaxis], col_ind]


class TransformerMixin(_SetOutputMixin):
    """Mixin class for all transformers in scikit-learn.

    This mixin defines the following functionality:

    - a `fit_transform` method that delegates to `fit` and `transform`;
    - a `set_output` method to output `X` as a specific container type.

    If :term:`get_feature_names_out` is defined, then :class:`BaseEstimator` will
    automatically wrap `transform` and `fit_transform` to follow the `set_output`
    API. See the :ref:`developer_api_set_output` for details.

    :class:`OneToOneFeatureMixin` and
    :class:`ClassNamePrefixFeaturesOutMixin` are helpful mixins for
    defining :term:`get_feature_names_out`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator, TransformerMixin
    >>> class MyTransformer(TransformerMixin, BaseEstimator):
    ...     def __init__(self, *, param=1):
    ...         self.param = param
    ...     def fit(self, X, y=None):
    ...         return self
    ...     def transform(self, X):
    ...         return np.full(shape=len(X), fill_value=self.param)
    >>> transformer = MyTransformer()
    >>> X = [[1, 2], [2, 3], [3, 4]]
    >>> transformer.fit_transform(X)
    array([1, 1, 1])
    """

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.transformer_tags = TransformerTags()
        return tags

    def fit_transform(self, X, y=None, **fit_params):
        """
        Fit to data, then transform it.

        Fits transformer to `X` and `y` with optional parameters `fit_params`
        and returns a transformed version of `X`.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input samples.

        y :  array-like of shape (n_samples,) or (n_samples, n_outputs), \
                default=None
            Target values (None for unsupervised transformations).

        **fit_params : dict
            Additional fit parameters.
            Pass only if the estimator accepts additional params in its `fit` method.

        Returns
        -------
        X_new : ndarray array of shape (n_samples, n_features_new)
            Transformed array.
        """
        # non-optimized default implementation; override when a better
        # method is possible for a given clustering algorithm

        # we do not route parameters here, since consumers don't route. But
        # since it's possible for a `transform` method to also consume
        # metadata, we check if that's the case, and we raise a warning telling
        # users that they should implement a custom `fit_transform` method
        # to forward metadata to `transform` as well.
        #
        # For that, we calculate routing and check if anything would be routed
        # to `transform` if we were to route them.
        if _routing_enabled():
            transform_params = self.get_metadata_routing().consumes(
                method="transform", params=fit_params.keys()
            )
            if transform_params:
                warnings.warn(
                    (
                        f"This object ({self.__class__.__name__}) has a `transform`"
                        " method which consumes metadata, but `fit_transform` does not"
                        " forward metadata to `transform`. Please implement a custom"
                        " `fit_transform` method to forward metadata to `transform` as"
                        " well. Alternatively, you can explicitly do"
                        " `set_transform_request`and set all values to `False` to"
                        " disable metadata routed to `transform`, if that's an option."
                    ),
                    UserWarning,
                )

        if y is None:
            # fit method of arity 1 (unsupervised transformation)
            return self.fit(X, **fit_params).transform(X)
        else:
            # fit method of arity 2 (supervised transformation)
            return self.fit(X, y, **fit_params).transform(X)


class OneToOneFeatureMixin:
    """Provides `get_feature_names_out` for simple transformers.

    This mixin assumes there's a 1-to-1 correspondence between input features
    and output features, such as :class:`~sklearn.preprocessing.StandardScaler`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import OneToOneFeatureMixin, BaseEstimator
    >>> class MyEstimator(OneToOneFeatureMixin, BaseEstimator):
    ...     def fit(self, X, y=None):
    ...         self.n_features_in_ = X.shape[1]
    ...         return self
    >>> X = np.array([[1, 2], [3, 4]])
    >>> MyEstimator().fit(X).get_feature_names_out()
    array(['x0', 'x1'], dtype=object)
    """

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Input features.

            - If `input_features` is `None`, then `feature_names_in_` is
              used as feature names in. If `feature_names_in_` is not defined,
              then the following input feature names are generated:
              `["x0", "x1", ..., "x(n_features_in_ - 1)"]`.
            - If `input_features` is an array-like, then `input_features` must
              match `feature_names_in_` if `feature_names_in_` is defined.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Same as input features.
        """
        # Note that passing attributes="n_features_in_" forces check_is_fitted
        # to check if the attribute is present. Otherwise it will pass on
        # stateless estimators (requires_fit=False)
        check_is_fitted(self, attributes="n_features_in_")
        return _check_feature_names_in(self, input_features)


class ClassNamePrefixFeaturesOutMixin:
    """Mixin class for transformers that generate their own names by prefixing.

    This mixin is useful when the transformer needs to generate its own feature
    names out, such as :class:`~sklearn.decomposition.PCA`. For example, if
    :class:`~sklearn.decomposition.PCA` outputs 3 features, then the generated feature
    names out are: `["pca0", "pca1", "pca2"]`.

    This mixin assumes that a `_n_features_out` attribute is defined when the
    transformer is fitted. `_n_features_out` is the number of output features
    that the transformer will return in `transform` of `fit_transform`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import ClassNamePrefixFeaturesOutMixin, BaseEstimator
    >>> class MyEstimator(ClassNamePrefixFeaturesOutMixin, BaseEstimator):
    ...     def fit(self, X, y=None):
    ...         self._n_features_out = X.shape[1]
    ...         return self
    >>> X = np.array([[1, 2], [3, 4]])
    >>> MyEstimator().fit(X).get_feature_names_out()
    array(['myestimator0', 'myestimator1'], dtype=object)
    """

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        The feature names out will prefixed by the lowercased class name. For
        example, if the transformer outputs 3 features, then the feature names
        out are: `["class_name0", "class_name1", "class_name2"]`.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Only used to validate feature names with the names seen in `fit`.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names.
        """
        check_is_fitted(self, "_n_features_out")
        return _generate_get_feature_names_out(
            self, self._n_features_out, input_features=input_features
        )


class DensityMixin:
    """Mixin class for all density estimators in scikit-learn.

    This mixin defines the following functionality:

    - sets estimator type to `"density_estimator"` through the `estimator_type` tag;
    - `score` method that default that do no-op.

    Examples
    --------
    >>> from sklearn.base import DensityMixin
    >>> class MyEstimator(DensityMixin):
    ...     def fit(self, X, y=None):
    ...         self.is_fitted_ = True
    ...         return self
    >>> estimator = MyEstimator()
    >>> hasattr(estimator, "score")
    True
    """

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.estimator_type = "density_estimator"
        return tags

    def score(self, X, y=None):
        """Return the score of the model on the data `X`.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Test samples.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        score : float
        """
        pass


class OutlierMixin:
    """Mixin class for all outlier detection estimators in scikit-learn.

    This mixin defines the following functionality:

    - set estimator type to `"outlier_detector"` through the `estimator_type` tag;
    - `fit_predict` method that default to `fit` and `predict`.

    Examples
    --------
    >>> import numpy as np
    >>> from sklearn.base import BaseEstimator, OutlierMixin
    >>> class MyEstimator(OutlierMixin):
    ...     def fit(self, X, y=None):
    ...         self.is_fitted_ = True
    ...         return self
    ...     def predict(self, X):
    ...         return np.ones(shape=len(X))
    >>> estimator = MyEstimator()
    >>> X = np.array([[1, 2], [2, 3], [3, 4]])
    >>> estimator.fit_predict(X)
    array([1., 1., 1.])
    """

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.estimator_type = "outlier_detector"
        return tags

    def fit_predict(self, X, y=None, **kwargs):
        """Perform fit on X and returns labels for X.

        Returns -1 for outliers and 1 for inliers.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            The input samples.

        y : Ignored
            Not used, present for API consistency by convention.

        **kwargs : dict
            Arguments to be passed to ``fit``.

            .. versionadded:: 1.4

        Returns
        -------
        y : ndarray of shape (n_samples,)
            1 for inliers, -1 for outliers.
        """
        # we do not route parameters here, since consumers don't route. But
        # since it's possible for a `predict` method to also consume
        # metadata, we check if that's the case, and we raise a warning telling
        # users that they should implement a custom `fit_predict` method
        # to forward metadata to `predict` as well.
        #
        # For that, we calculate routing and check if anything would be routed
        # to `predict` if we were to route them.
        if _routing_enabled():
            transform_params = self.get_metadata_routing().consumes(
                method="predict", params=kwargs.keys()
            )
            if transform_params:
                warnings.warn(
                    (
                        f"This object ({self.__class__.__name__}) has a `predict` "
                        "method which consumes metadata, but `fit_predict` does not "
                        "forward metadata to `predict`. Please implement a custom "
                        "`fit_predict` method to forward metadata to `predict` as well."
                        "Alternatively, you can explicitly do `set_predict_request`"
                        "and set all values to `False` to disable metadata routed to "
                        "`predict`, if that's an option."
                    ),
                    UserWarning,
                )

        # override for transductive outlier detectors like LocalOulierFactor
        return self.fit(X, **kwargs).predict(X)


class MetaEstimatorMixin:
    """Mixin class for all meta estimators in scikit-learn.

    This mixin is empty, and only exists to indicate that the estimator is a
    meta-estimator.

    .. versionchanged:: 1.6
        The `_required_parameters` is now removed and is unnecessary since tests are
        refactored and don't use this anymore.

    Examples
    --------
    >>> from sklearn.base import MetaEstimatorMixin
    >>> from sklearn.datasets import load_iris
    >>> from sklearn.linear_model import LogisticRegression
    >>> class MyEstimator(MetaEstimatorMixin):
    ...     def __init__(self, *, estimator=None):
    ...         self.estimator = estimator
    ...     def fit(self, X, y=None):
    ...         if self.estimator is None:
    ...             self.estimator_ = LogisticRegression()
    ...         else:
    ...             self.estimator_ = self.estimator
    ...         return self
    >>> X, y = load_iris(return_X_y=True)
    >>> estimator = MyEstimator().fit(X, y)
    >>> estimator.estimator_
    LogisticRegression()
    """


class MultiOutputMixin:
    """Mixin to mark estimators that support multioutput."""

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.target_tags.multi_output = True
        return tags


class _UnstableArchMixin:
    """Mark estimators that are non-determinstic on 32bit or PowerPC"""

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        tags.non_deterministic = _IS_32BIT or platform.machine().startswith(
            ("ppc", "powerpc")
        )
        return tags


def is_classifier(estimator):
    """Return True if the given estimator is (probably) a classifier.

    Parameters
    ----------
    estimator : estimator instance
        Estimator object to test.

    Returns
    -------
    out : bool
        True if estimator is a classifier and False otherwise.

    Examples
    --------
    >>> from sklearn.base import is_classifier
    >>> from sklearn.cluster import KMeans
    >>> from sklearn.svm import SVC, SVR
    >>> classifier = SVC()
    >>> regressor = SVR()
    >>> kmeans = KMeans()
    >>> is_classifier(classifier)
    True
    >>> is_classifier(regressor)
    False
    >>> is_classifier(kmeans)
    False
    """
    return get_tags(estimator).estimator_type == "classifier"


def is_regressor(estimator):
    """Return True if the given estimator is (probably) a regressor.

    Parameters
    ----------
    estimator : estimator instance
        Estimator object to test.

    Returns
    -------
    out : bool
        True if estimator is a regressor and False otherwise.

    Examples
    --------
    >>> from sklearn.base import is_regressor
    >>> from sklearn.cluster import KMeans
    >>> from sklearn.svm import SVC, SVR
    >>> classifier = SVC()
    >>> regressor = SVR()
    >>> kmeans = KMeans()
    >>> is_regressor(classifier)
    False
    >>> is_regressor(regressor)
    True
    >>> is_regressor(kmeans)
    False
    """
    return get_tags(estimator).estimator_type == "regressor"


def is_clusterer(estimator):
    """Return True if the given estimator is (probably) a clusterer.

    .. versionadded:: 1.6

    Parameters
    ----------
    estimator : estimator instance
        Estimator object to test.

    Returns
    -------
    out : bool
        True if estimator is a clusterer and False otherwise.

    Examples
    --------
    >>> from sklearn.base import is_clusterer
    >>> from sklearn.cluster import KMeans
    >>> from sklearn.svm import SVC, SVR
    >>> classifier = SVC()
    >>> regressor = SVR()
    >>> kmeans = KMeans()
    >>> is_clusterer(classifier)
    False
    >>> is_clusterer(regressor)
    False
    >>> is_clusterer(kmeans)
    True
    """
    return get_tags(estimator).estimator_type == "clusterer"


def is_outlier_detector(estimator):
    """Return True if the given estimator is (probably) an outlier detector.

    Parameters
    ----------
    estimator : estimator instance
        Estimator object to test.

    Returns
    -------
    out : bool
        True if estimator is an outlier detector and False otherwise.
    """
    return get_tags(estimator).estimator_type == "outlier_detector"


def _fit_context(*, prefer_skip_nested_validation):
    """Decorator to run the fit methods of estimators within context managers.

    Parameters
    ----------
    prefer_skip_nested_validation : bool
        If True, the validation of parameters of inner estimators or functions
        called during fit will be skipped.

        This is useful to avoid validating many times the parameters passed by the
        user from the public facing API. It's also useful to avoid validating
        parameters that we pass internally to inner functions that are guaranteed to
        be valid by the test suite.

        It should be set to True for most estimators, except for those that receive
        non-validated objects as parameters, such as meta-estimators that are given
        estimator objects.

    Returns
    -------
    decorated_fit : method
        The decorated fit method.
    """

    def decorator(fit_method):
        @functools.wraps(fit_method)
        def wrapper(estimator, *args, **kwargs):
            global_skip_validation = get_config()["skip_parameter_validation"]

            # we don't want to validate again for each call to partial_fit
            partial_fit_and_fitted = (
                fit_method.__name__ == "partial_fit" and _is_fitted(estimator)
            )

            if not global_skip_validation and not partial_fit_and_fitted:
                estimator._validate_params()

            with config_context(
                skip_parameter_validation=(
                    prefer_skip_nested_validation or global_skip_validation
                )
            ):
                return fit_method(estimator, *args, **kwargs)

        return wrapper

    return decorator
