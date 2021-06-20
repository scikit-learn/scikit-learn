"""Base classes for all estimators."""

# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# License: BSD 3 clause

import copy
import warnings
from collections import defaultdict
import platform
import inspect
import re
import numpy as np
from typing import Union, Optional

from . import __version__
from .utils import _IS_32BIT
from .utils import _standardize_metadata_request
from ._config import get_config
from .utils._tags import (
    _DEFAULT_TAGS,
    _safe_tags,
)
from .utils.validation import check_X_y
from .utils.validation import check_array
from .utils._estimator_html_repr import estimator_html_repr
from .utils.validation import _deprecate_positional_args
from .utils.metadata_requests import MetadataRequest
from .utils.metadata_requests import metadata_request_factory
from .utils.metadata_requests import UNCHANGED
from .utils.metadata_requests import RequestType


@_deprecate_positional_args
def clone(estimator, *, safe=True):
    """Constructs a new unfitted estimator with the same parameters.

    Clone does a deep copy of the model in an estimator
    without actually copying attached data. It yields a new estimator
    with the same parameters that has not been fitted on any data.

    If the estimator's `random_state` parameter is an integer (or if the
    estimator doesn't have a `random_state` parameter), an *exact clone* is
    returned: the clone and the original estimator will give the exact same
    results. Otherwise, *statistical clone* is returned: the clone might
    yield different results from the original estimator. More details can be
    found in :ref:`randomness`.

    Parameters
    ----------
    estimator : {list, tuple, set} of estimator instance or a single \
            estimator instance
        The estimator or group of estimators to be cloned.

    safe : bool, default=True
        If safe is False, clone will fall back to a deep copy on objects
        that are not estimators.

    """
    estimator_type = type(estimator)
    # XXX: not handling dictionaries
    if estimator_type in (list, tuple, set, frozenset):
        return estimator_type([clone(e, safe=safe) for e in estimator])
    elif not hasattr(estimator, 'get_params') or isinstance(estimator, type):
        if not safe:
            return copy.deepcopy(estimator)
        else:
            if isinstance(estimator, type):
                raise TypeError("Cannot clone object. " +
                                "You should provide an instance of " +
                                "scikit-learn estimator instead of a class.")
            else:
                raise TypeError("Cannot clone object '%s' (type %s): "
                                "it does not seem to be a scikit-learn "
                                "estimator as it does not implement a "
                                "'get_params' method."
                                % (repr(estimator), type(estimator)))

    klass = estimator.__class__
    new_object_params = estimator.get_params(deep=False)
    for name, param in new_object_params.items():
        new_object_params[name] = clone(param, safe=False)

    new_object = klass(**new_object_params)
    try:
        new_object._metadata_request = copy.deepcopy(
            estimator._metadata_request)
    except AttributeError:
        pass

    params_set = new_object.get_params(deep=False)

    # quick sanity check of the parameters of the clone
    for name in new_object_params:
        param1 = new_object_params[name]
        param2 = params_set[name]
        if param1 is not param2:
            raise RuntimeError('Cannot clone object %s, as the constructor '
                               'either does not set or modifies parameter %s' %
                               (estimator, name))
    return new_object


def _pprint(params, offset=0, printer=repr):
    """Pretty print the dictionary 'params'

    Parameters
    ----------
    params : dict
        The dictionary to pretty print

    offset : int, default=0
        The offset in characters to add at the begin of each line.

    printer : callable, default=repr
        The function to convert entries to strings, typically
        the builtin str or repr

    """
    # Do a multi-line justified repr:
    options = np.get_printoptions()
    np.set_printoptions(precision=5, threshold=64, edgeitems=2)
    params_list = list()
    this_line_length = offset
    line_sep = ',\n' + (1 + offset // 2) * ' '
    for i, (k, v) in enumerate(sorted(params.items())):
        if type(v) is float:
            # use str for representing floating point numbers
            # this way we get consistent representation across
            # architectures and versions.
            this_repr = '%s=%s' % (k, str(v))
        else:
            # use repr of the rest
            this_repr = '%s=%s' % (k, printer(v))
        if len(this_repr) > 500:
            this_repr = this_repr[:300] + '...' + this_repr[-100:]
        if i > 0:
            if (this_line_length + len(this_repr) >= 75 or '\n' in this_repr):
                params_list.append(line_sep)
                this_line_length = len(line_sep)
            else:
                params_list.append(', ')
                this_line_length += 2
        params_list.append(this_repr)
        this_line_length += len(this_repr)

    np.set_printoptions(**options)
    lines = ''.join(params_list)
    # Strip trailing space to avoid nightmare in doctests
    lines = '\n'.join(l.rstrip(' ') for l in lines.split('\n'))
    return lines


class RequestMethod:
    """
    A descriptor to create a function to be set as request methods.

    Parameters
    ----------
    name : str
        The name of the method for which the request function should be
        created, e.g. ``"fit"`` would create a ``fit_requests`` function.

    keys : list of str
        A list of strings which are accepted parameters by the created
        function, e.g. ``["sample_weight"]`` if the corresponding method
        accepts it as a metadata.

    Notes
    -----
    This class is a descriptor [1]_ and uses PEP-362 to set the signature of
    the returned function [2]_.

    References
    ----------
    .. [1] https://docs.python.org/3/howto/descriptor.html

    .. [2] https://www.python.org/dev/peps/pep-0362/
    """
    def __init__(self, name, keys):
        self.name = name
        self.keys = keys

    def __get__(self, instance, owner):
        # we would want to have a method which accepts only the expected args
        def func(**kw):
            if set(kw) - set(self.keys):
                raise TypeError(f"Unexpected args: {set(kw) - set(self.keys)}")

            requests = metadata_request_factory(instance)

            try:
                method_metadata_request = getattr(requests, self.name)
            except AttributeError:
                raise ValueError(f"{self.name} is not a supported method.")

            for prop, alias in kw.items():
                if alias is not UNCHANGED:
                    method_metadata_request.add_request(prop=prop, alias=alias)
            instance._metadata_request = requests.to_dict()

            return instance

        # Now we set the relevant attributes of the function so that it seems
        # like a normal method to the end user, with known expected arguments.
        func.__name__ = f"{self.name}_requests"
        func.__signature__ = inspect.Signature([
            inspect.Parameter(k,
                              inspect.Parameter.KEYWORD_ONLY,
                              default=UNCHANGED,
                              annotation=Optional[Union[RequestType, str]])
            for k in self.keys
        ], return_annotation=type(instance))
        func.__doc__ = "Something useful"
        return func


class _MetadataConsumer:
    def __init_subclass__(cls, **kwargs):
        """Set the ``{method}_requests`` methods.

        This uses PEP-487 [1]_ to set the ``{method}_requests`` methods. It
        looks for the information available in the set default values which are
        set using ``_metadata_request__*`` class attributes.

        References
        ----------
        .. [1] https://www.python.org/dev/peps/pep-0487
        """
        try:
            requests = cls._get_default_requests().to_dict()
        except Exception:
            # if there are any issues in the default values, it will be raised
            # when ``get_metadata_request`` is called. Here we are going to
            # ignore all the issues such as bad defaults etc.`
            super().__init_subclass__(**kwargs)
            return

        for request_method, request_keys in requests.items():
            # set ``{method}_requests``` methods
            if not len(request_keys):
                continue
            setattr(
                cls,
                f"{request_method}_requests",
                RequestMethod(request_method, sorted(request_keys))
            )
        super().__init_subclass__(**kwargs)

    @classmethod
    def _get_default_requests(cls):
        """Collect default request values.

        This method combines the information present in ``metadata_request__*``
        class attributes.
        """
        requests = MetadataRequest()

        # need to go through the MRO since this is a class attribute and
        # ``vars`` doesn't report the parent class attributes. We go through
        # the reverse of the MRO since cls is the first in the tuple and object
        # is the last.
        defaults = defaultdict()
        for klass in reversed(inspect.getmro(cls)):
            klass_defaults = {
                attr: value for attr, value in vars(klass).items()
                if attr.startswith("_metadata_request__")
            }
            defaults.update(klass_defaults)
        defaults = dict(sorted(defaults.items()))
        for attr, value in defaults.items():
            requests.add_requests(
                value,
                expected_metadata="__".join(attr.split("__")[1:])
            )
        return requests

    def get_metadata_request(self, output="dict"):
        """Get requested data properties.

        Parameters
        ----------
        output : {"dict", "MetadataRequest}
            Whether the output should be a MetadataRequest instance, or a dict
            representing that instance.

        Returns
        -------
        request : MetadataRequest, or dict
            If dict, it will be a deserialized version of the underlying
            MetadataRequest object: dict of dict of str->value. The key to the
            first dict is the name of the method, and the key to the second
            dict is the name of the argument requested by the method.
        """
        if output not in {"dict", "MetadataRequest"}:
            raise ValueError(
                "output can be one of {'dict', 'MetadataRequest'}."
            )
        if hasattr(self, "_metadata_request"):
            requests = MetadataRequest(self._metadata_request)
        else:
            requests = self._get_default_requests()
            self._metadata_request = requests

        if output == "dict":
            return requests.to_dict()
        else:
            return requests


class MetadataConsumer:
    def _remove_param(self, *, method_params, param):
        user_provides_list = list(method_params.keys())
        for user_provides in user_provides_list:
            method_params[user_provides] -= {param}
            if not len(method_params[user_provides]):
                del method_params[user_provides]

    def _request_key_for_method(self, *, method, param, user_provides):
        if user_provides is None:
            return
        if not hasattr(self, "_metadata_request"):
            self._metadata_request = self.get_metadata_request()
        self._metadata_request = _standardize_metadata_request(
            self._metadata_request
        )

        if user_provides == param:
            user_provides = True

        method_dict = self._metadata_request[method]
        self._remove_param(method_params=method_dict, param=param)

        if user_provides is True:
            method_dict[param] = {param}
        elif user_provides is False:
            # the parameter is already removed, nothing to do here
            pass
        elif isinstance(user_provides, str):
            routing = method_dict.get(user_provides, set())
            method_dict[user_provides] = routing.union({param})
        else:
            raise TypeError

    def _set_metadata_request(self, props):
        if props is None:
            try:
                del self._metadata_request
            except AttributeError:
                pass
            return self

        self._metadata_request = _standardize_metadata_request(props)
        return self


class SampleWeightConsumer:
    _metadata_request__sample_weight = {
        "fit": {"sample_weight": None},
        "score": {"sample_weight": None},
    }

    def request_sample_weight(self, *, fit=None, score=None):
        """Define how to receive sample_weight from a parent meta-estimator

        Parameters
        ----------
        fit : string or bool, default=None
            The fit parameter name that a meta-estimator should pass to this
            estimator as sample_weight. If true, the name will be
            'sample_weight'.
            If False, no fit parameter will be passed.
            If None, no change in routing.

        score : string or bool, default=None
            The parameter name that a meta-estimator should pass to this
            estimator as sample_weight when calling its `score` method.
            If true, the name will be 'sample_weight'.
            If False, no score parameter will be passed.
            If None, no change in routing.

        Returns
        -------
        self

        Examples
        --------
        >>> from sklearn.linear_model import LogisticRegression
        >>> from sklearn.model_selection import GridSearchCV
        >>> from sklearn.datasets import load_iris
        >>> X, y = load_iris(return_X_y=True)
        >>> sample_weight = np.random.rand(len(X))
        >>> clf = LogisticRegression()
        >>> gs = GridSearchCV(clf, {"C": [.1, 1, 10]})
        >>> # Unweighted fitting and scoring
        >>> gs.fit(X, y)
        GridSearchCV(...)
        >>> # Weighted fitting, unweighted scoring
        >>> clf.request_sample_weight(fit=True, score=False)
        LogisticRegression()
        >>> gs.fit(X, y, sample_weight=sample_weight)
        GridSearchCV(...)
        >>> # Weighted fitting and scoring
        >>> clf.request_sample_weight(fit=True, score=True)
        LogisticRegression()
        >>> gs.fit(X, y, sample_weight=sample_weight)
        GridSearchCV(...)
        >>> # Weighted scoring only
        >>> clf.request_sample_weight(fit=False, score=True)
        LogisticRegression()
        >>> gs.fit(X, y, sample_weight=sample_weight)
        GridSearchCV(...)
        >>> # Distinct weights for fit and score
        >>> score_sample_weight = np.random.rand(len(X))
        >>> clf.request_sample_weight(fit='fit_sample_weight',
        ...     score='score_sample_weight')
        LogisticRegression()
        >>> gs.fit(X, y, fit_sample_weight=sample_weight,
        ...        score_sample_weight=score_sample_weight)
        GridSearchCV(...)
        """
        self._request_key_for_method(
            method="fit", param="sample_weight", user_provides=fit
        )
        self._request_key_for_method(
            method="score", param="sample_weight", user_provides=score
        )
        return self


class BaseEstimator(_MetadataConsumer):
    """Base class for all estimators in scikit-learn

    Notes
    -----
    All estimators should specify all the parameters that can be set
    at the class level in their ``__init__`` as explicit keyword
    arguments (no ``*args`` or ``**kwargs``).
    """

    @classmethod
    def _get_param_names(cls):
        """Get parameter names for the estimator"""
        # fetch the constructor or the original constructor before
        # deprecation wrapping if any
        init = getattr(cls.__init__, 'deprecated_original', cls.__init__)
        if init is object.__init__:
            # No explicit constructor to introspect
            return []

        # introspect the constructor arguments to find the model parameters
        # to represent
        init_signature = inspect.signature(init)
        # Consider the constructor parameters excluding 'self'
        parameters = [p for p in init_signature.parameters.values()
                      if p.name != 'self' and p.kind != p.VAR_KEYWORD]
        for p in parameters:
            if p.kind == p.VAR_POSITIONAL:
                raise RuntimeError("scikit-learn estimators should always "
                                   "specify their parameters in the signature"
                                   " of their __init__ (no varargs)."
                                   " %s with constructor %s doesn't "
                                   " follow this convention."
                                   % (cls, init_signature))
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
            if deep and hasattr(value, 'get_params'):
                deep_items = value.get_params().items()
                out.update((key + '__' + k, val) for k, val in deep_items)
            out[key] = value
        return out

    def set_params(self, **params):
        """
        Set the parameters of this estimator.

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
            key, delim, sub_key = key.partition('__')
            if key not in valid_params:
                raise ValueError('Invalid parameter %s for estimator %s. '
                                 'Check the list of available parameters '
                                 'with `estimator.get_params().keys()`.' %
                                 (key, self))

            if delim:
                nested_params[key][sub_key] = value
            else:
                setattr(self, key, value)
                valid_params[key] = value

        for key, sub_params in nested_params.items():
            valid_params[key].set_params(**sub_params)

        return self

    def __repr__(self, N_CHAR_MAX=700):
        # N_CHAR_MAX is the (approximate) maximum number of non-blank
        # characters to render. We pass it as an optional parameter to ease
        # the tests.

        from .utils._pprint import _EstimatorPrettyPrinter

        N_MAX_ELEMENTS_TO_SHOW = 30  # number of elements to show in sequences

        # use ellipsis for sequences with a lot of elements
        pp = _EstimatorPrettyPrinter(
            compact=True, indent=1, indent_at_name=True,
            n_max_elements_to_show=N_MAX_ELEMENTS_TO_SHOW)

        repr_ = pp.pformat(self)

        # Use bruteforce ellipsis when there are a lot of non-blank characters
        n_nonblank = len(''.join(repr_.split()))
        if n_nonblank > N_CHAR_MAX:
            lim = N_CHAR_MAX // 2  # apprx number of chars to keep on both ends
            regex = r'^(\s*\S){%d}' % lim
            # The regex '^(\s*\S){%d}' % n
            # matches from the start of the string until the nth non-blank
            # character:
            # - ^ matches the start of string
            # - (pattern){n} matches n repetitions of pattern
            # - \s*\S matches a non-blank char following zero or more blanks
            left_lim = re.match(regex, repr_).end()
            right_lim = re.match(regex, repr_[::-1]).end()

            if '\n' in repr_[left_lim:-right_lim]:
                # The left side and right side aren't on the same line.
                # To avoid weird cuts, e.g.:
                # categoric...ore',
                # we need to start the right side with an appropriate newline
                # character so that it renders properly as:
                # categoric...
                # handle_unknown='ignore',
                # so we add [^\n]*\n which matches until the next \n
                regex += r'[^\n]*\n'
                right_lim = re.match(regex, repr_[::-1]).end()

            ellipsis = '...'
            if left_lim + len(ellipsis) < len(repr_) - right_lim:
                # Only add ellipsis if it results in a shorter repr
                repr_ = repr_[:left_lim] + '...' + repr_[-right_lim:]

        return repr_

    def __getstate__(self):
        try:
            state = super().__getstate__()
        except AttributeError:
            state = self.__dict__.copy()

        if type(self).__module__.startswith('sklearn.'):
            return dict(state.items(), _sklearn_version=__version__)
        else:
            return state

    def __setstate__(self, state):
        if type(self).__module__.startswith('sklearn.'):
            pickle_version = state.pop("_sklearn_version", "pre-0.18")
            if pickle_version != __version__:
                warnings.warn(
                    "Trying to unpickle estimator {0} from version {1} when "
                    "using version {2}. This might lead to breaking code or "
                    "invalid results. Use at your own risk.".format(
                        self.__class__.__name__, pickle_version, __version__),
                    UserWarning)
        try:
            super().__setstate__(state)
        except AttributeError:
            self.__dict__.update(state)

    def _more_tags(self):
        return _DEFAULT_TAGS

    def _get_tags(self):
        collected_tags = {}
        for base_class in reversed(inspect.getmro(self.__class__)):
            if hasattr(base_class, '_more_tags'):
                # need the if because mixins might not have _more_tags
                # but might do redundant work in estimators
                # (i.e. calling more tags on BaseEstimator multiple times)
                more_tags = base_class._more_tags(self)
                collected_tags.update(more_tags)
        return collected_tags

    def _check_n_features(self, X, reset):
        """Set the `n_features_in_` attribute, or check against it.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input samples.
        reset : bool
            If True, the `n_features_in_` attribute is set to `X.shape[1]`.
            If False and the attribute exists, then check that it is equal to
            `X.shape[1]`. If False and the attribute does *not* exist, then
            the check is skipped.
            .. note::
               It is recommended to call reset=True in `fit` and in the first
               call to `partial_fit`. All other methods that validate `X`
               should set `reset=False`.
        """
        n_features = X.shape[1]

        if reset:
            self.n_features_in_ = n_features
            return

        if not hasattr(self, "n_features_in_"):
            # Skip this check if the expected number of expected input features
            # was not recorded by calling fit first. This is typically the case
            # for stateless transformers.
            return

        if n_features != self.n_features_in_:
            raise ValueError(
                f"X has {n_features} features, but {self.__class__.__name__} "
                f"is expecting {self.n_features_in_} features as input.")

    def _validate_data(self, X, y='no_validation', reset=True,
                       validate_separately=False, **check_params):
        """Validate input data and set or check the `n_features_in_` attribute.

        Parameters
        ----------
        X : {array-like, sparse matrix, dataframe} of shape \
                (n_samples, n_features)
            The input samples.
        y : array-like of shape (n_samples,), default='no_validation'
            The targets.

            - If `None`, `check_array` is called on `X`. If the estimator's
              requires_y tag is True, then an error will be raised.
            - If `'no_validation'`, `check_array` is called on `X` and the
              estimator's requires_y tag is ignored. This is a default
              placeholder and is never meant to be explicitly set.
            - Otherwise, both `X` and `y` are checked with either `check_array`
              or `check_X_y` depending on `validate_separately`.

        reset : bool, default=True
            Whether to reset the `n_features_in_` attribute.
            If False, the input will be checked for consistency with data
            provided when reset was last True.
            .. note::
               It is recommended to call reset=True in `fit` and in the first
               call to `partial_fit`. All other methods that validate `X`
               should set `reset=False`.
        validate_separately : False or tuple of dicts, default=False
            Only used if y is not None.
            If False, call validate_X_y(). Else, it must be a tuple of kwargs
            to be used for calling check_array() on X and y respectively.
        **check_params : kwargs
            Parameters passed to :func:`sklearn.utils.check_array` or
            :func:`sklearn.utils.check_X_y`. Ignored if validate_separately
            is not False.

        Returns
        -------
        out : {ndarray, sparse matrix} or tuple of these
            The validated input. A tuple is returned if `y` is not None.
        """

        if y is None:
            if self._get_tags()['requires_y']:
                raise ValueError(
                    f"This {self.__class__.__name__} estimator "
                    f"requires y to be passed, but the target y is None."
                )
            X = check_array(X, **check_params)
            out = X
        elif isinstance(y, str) and y == 'no_validation':
            X = check_array(X, **check_params)
            out = X
        else:
            if validate_separately:
                # We need this because some estimators validate X and y
                # separately, and in general, separately calling check_array()
                # on X and y isn't equivalent to just calling check_X_y()
                # :(
                check_X_params, check_y_params = validate_separately
                X = check_array(X, **check_X_params)
                y = check_array(y, **check_y_params)
            else:
                X, y = check_X_y(X, y, **check_params)
            out = X, y

        if check_params.get('ensure_2d', True):
            self._check_n_features(X, reset=reset)

        return out

    @property
    def _repr_html_(self):
        """HTML representation of estimator.

        This is redundant with the logic of `_repr_mimebundle_`. The latter
        should be favorted in the long term, `_repr_html_` is only
        implemented for consumers who do not interpret `_repr_mimbundle_`.
        """
        if get_config()["display"] != 'diagram':
            raise AttributeError("_repr_html_ is only defined when the "
                                 "'display' configuration option is set to "
                                 "'diagram'")
        return self._repr_html_inner

    def _repr_html_inner(self):
        """This function is returned by the @property `_repr_html_` to make
        `hasattr(estimator, "_repr_html_") return `True` or `False` depending
        on `get_config()["display"]`.
        """
        return estimator_html_repr(self)

    def _repr_mimebundle_(self, **kwargs):
        """Mime bundle used by jupyter kernels to display estimator"""
        output = {"text/plain": repr(self)}
        if get_config()["display"] == 'diagram':
            output["text/html"] = estimator_html_repr(self)
        return output


class ClassifierMixin(SampleWeightConsumer):
    """Mixin class for all classifiers in scikit-learn."""

    _estimator_type = "classifier"

    def score(self, X, y, sample_weight=None):
        """
        Return the mean accuracy on the given test data and labels.

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
            Mean accuracy of ``self.predict(X)`` wrt. `y`.
        """
        from .metrics import accuracy_score
        return accuracy_score(y, self.predict(X), sample_weight=sample_weight)

    def _more_tags(self):
        return {'requires_y': True}


class RegressorMixin(SampleWeightConsumer):
    """Mixin class for all regression estimators in scikit-learn."""
    _estimator_type = "regressor"

    def score(self, X, y, sample_weight=None):
        """Return the coefficient of determination :math:`R^2` of the
        prediction.

        The coefficient :math:`R^2` is defined as :math:`(1 - \\frac{u}{v})`,
        where :math:`u` is the residual sum of squares ``((y_true - y_pred)
        ** 2).sum()`` and :math:`v` is the total sum of squares ``((y_true -
        y_true.mean()) ** 2).sum()``. The best possible score is 1.0 and it
        can be negative (because the model can be arbitrarily worse). A
        constant model that always predicts the expected value of `y`,
        disregarding the input features, would get a :math:`R^2` score of
        0.0.

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
            :math:`R^2` of ``self.predict(X)`` wrt. `y`.

        Notes
        -----
        The :math:`R^2` score used when calling ``score`` on a regressor uses
        ``multioutput='uniform_average'`` from version 0.23 to keep consistent
        with default value of :func:`~sklearn.metrics.r2_score`.
        This influences the ``score`` method of all the multioutput
        regressors (except for
        :class:`~sklearn.multioutput.MultiOutputRegressor`).
        """

        from .metrics import r2_score
        y_pred = self.predict(X)
        return r2_score(y, y_pred, sample_weight=sample_weight)

    def _more_tags(self):
        return {'requires_y': True}


class ClusterMixin:
    """Mixin class for all cluster estimators in scikit-learn."""
    _estimator_type = "clusterer"

    def fit_predict(self, X, y=None):
        """
        Perform clustering on `X` and returns cluster labels.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        labels : ndarray of shape (n_samples,), dtype=np.int64
            Cluster labels.
        """
        # non-optimized default implementation; override when a better
        # method is possible for a given clustering algorithm
        self.fit(X)
        return self.labels_

    def _more_tags(self):
        return {"preserves_dtype": []}


class BiclusterMixin:
    """Mixin class for all bicluster estimators in scikit-learn."""

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
        from .utils.validation import check_array
        data = check_array(data, accept_sparse='csr')
        row_ind, col_ind = self.get_indices(i)
        return data[row_ind[:, np.newaxis], col_ind]


class TransformerMixin:
    """Mixin class for all transformers in scikit-learn."""

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

        Returns
        -------
        X_new : ndarray array of shape (n_samples, n_features_new)
            Transformed array.
        """
        # non-optimized default implementation; override when a better
        # method is possible for a given clustering algorithm
        if y is None:
            # fit method of arity 1 (unsupervised transformation)
            return self.fit(X, **fit_params).transform(X)
        else:
            # fit method of arity 2 (supervised transformation)
            return self.fit(X, y, **fit_params).transform(X)


class DensityMixin:
    """Mixin class for all density estimators in scikit-learn."""
    _estimator_type = "DensityEstimator"

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
    """Mixin class for all outlier detection estimators in scikit-learn."""
    _estimator_type = "outlier_detector"

    def fit_predict(self, X, y=None):
        """Perform fit on X and returns labels for X.

        Returns -1 for outliers and 1 for inliers.

        Parameters
        ----------
        X : {array-like, sparse matrix, dataframe} of shape \
            (n_samples, n_features)

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        y : ndarray of shape (n_samples,)
            1 for inliers, -1 for outliers.
        """
        # override for transductive outlier detectors like LocalOulierFactor
        return self.fit(X).predict(X)


class MetaEstimatorMixin:
    _required_parameters = ["estimator"]
    """Mixin class for all meta estimators in scikit-learn."""


class MultiOutputMixin:
    """Mixin to mark estimators that support multioutput."""
    def _more_tags(self):
        return {'multioutput': True}


class _UnstableArchMixin:
    """Mark estimators that are non-determinstic on 32bit or PowerPC"""
    def _more_tags(self):
        return {'non_deterministic': (
            _IS_32BIT or platform.machine().startswith(('ppc', 'powerpc')))}


def is_classifier(estimator):
    """Return True if the given estimator is (probably) a classifier.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    Returns
    -------
    out : bool
        True if estimator is a classifier and False otherwise.
    """
    return getattr(estimator, "_estimator_type", None) == "classifier"


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
    """
    return getattr(estimator, "_estimator_type", None) == "regressor"


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
    return getattr(estimator, "_estimator_type", None) == "outlier_detector"


def _is_pairwise(estimator):
    """Returns True if estimator is pairwise.

    - If the `_pairwise` attribute and the tag are present and consistent,
      then use the value and not issue a warning.
    - If the `_pairwise` attribute and the tag are present and not
      consistent, use the `_pairwise` value and issue a deprecation
      warning.
    - If only the `_pairwise` attribute is present and it is not False,
      issue a deprecation warning and use the `_pairwise` value.

    Parameters
    ----------
    estimator : object
        Estimator object to test.

    Returns
    -------
    out : bool
        True if the estimator is pairwise and False otherwise.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=FutureWarning)
        has_pairwise_attribute = hasattr(estimator, '_pairwise')
        pairwise_attribute = getattr(estimator, '_pairwise', False)
    pairwise_tag = _safe_tags(estimator, key="pairwise")

    if has_pairwise_attribute:
        if pairwise_attribute != pairwise_tag:
            warnings.warn(
                "_pairwise was deprecated in 0.24 and will be removed in 1.1 "
                "(renaming of 0.26). Set the estimator tags of your estimator "
                "instead",
                FutureWarning
            )
        return pairwise_attribute

    # use pairwise tag when the attribute is not present
    return pairwise_tag
