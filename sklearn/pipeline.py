"""Utilities to build a composite estimator as a chain of transforms and estimators."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings
from collections import Counter, defaultdict
from contextlib import contextmanager
from copy import deepcopy
from itertools import chain, islice

import numpy as np
from scipy import sparse

from .base import TransformerMixin, _fit_context, clone
from .exceptions import NotFittedError
from .preprocessing import FunctionTransformer
from .utils import Bunch
from .utils._metadata_requests import METHODS
from .utils._param_validation import HasMethods, Hidden
from .utils._repr_html.estimator import _VisualBlock
from .utils._set_output import (
    _get_container_adapter,
    _safe_set_output,
)
from .utils._tags import get_tags
from .utils._user_interface import _print_elapsed_time
from .utils.metadata_routing import (
    MetadataRouter,
    MethodMapping,
    _raise_for_params,
    _routing_enabled,
    get_routing_for_object,
    process_routing,
)
from .utils.metaestimators import _BaseComposition, available_if
from .utils.parallel import Parallel, delayed
from .utils.validation import check_is_fitted, check_memory

__all__ = ["FeatureUnion", "Pipeline", "make_pipeline", "make_union"]


@contextmanager
def _raise_or_warn_if_not_fitted(estimator):
    """A context manager to make sure a NotFittedError is raised, if a sub-estimator
    raises the error.

    Otherwise, we raise a warning if the pipeline is not fitted, with the deprecation.

    TODO(1.8): remove this context manager and replace with check_is_fitted.
    """
    try:
        yield
    except NotFittedError as exc:
        raise NotFittedError("Pipeline is not fitted yet.") from exc

    # we only get here if the above didn't raise
    try:
        check_is_fitted(estimator)
    except NotFittedError:
        warnings.warn(
            "This Pipeline instance is not fitted yet. Call 'fit' with "
            "appropriate arguments before using other methods such as transform, "
            "predict, etc. This will raise an error in 1.8 instead of the current "
            "warning.",
            FutureWarning,
        )


def _final_estimator_has(attr):
    """Check that final_estimator has `attr`.

    Used together with `available_if` in `Pipeline`."""

    def check(self):
        # raise original `AttributeError` if `attr` does not exist
        getattr(self._final_estimator, attr)
        return True

    return check


def _cached_transform(
    sub_pipeline, *, cache, param_name, param_value, transform_params
):
    """Transform a parameter value using a sub-pipeline and cache the result.

    Parameters
    ----------
    sub_pipeline : Pipeline
        The sub-pipeline to be used for transformation.
    cache : dict
        The cache dictionary to store the transformed values.
    param_name : str
        The name of the parameter to be transformed.
    param_value : object
        The value of the parameter to be transformed.
    transform_params : dict
        The metadata to be used for transformation. This passed to the
        `transform` method of the sub-pipeline.

    Returns
    -------
    transformed_value : object
        The transformed value of the parameter.
    """
    if param_name not in cache:
        # If the parameter is a tuple, transform each element of the
        # tuple. This is needed to support the pattern present in
        # `lightgbm` and `xgboost` where users can pass multiple
        # validation sets.
        if isinstance(param_value, tuple):
            cache[param_name] = tuple(
                sub_pipeline.transform(element, **transform_params)
                for element in param_value
            )
        else:
            cache[param_name] = sub_pipeline.transform(param_value, **transform_params)

    return cache[param_name]


class Pipeline(_BaseComposition):
    """
    A sequence of data transformers with an optional final predictor.

    `Pipeline` allows you to sequentially apply a list of transformers to
    preprocess the data and, if desired, conclude the sequence with a final
    :term:`predictor` for predictive modeling.

    Intermediate steps of the pipeline must be transformers, that is, they
    must implement `fit` and `transform` methods.
    The final :term:`estimator` only needs to implement `fit`.
    The transformers in the pipeline can be cached using ``memory`` argument.

    The purpose of the pipeline is to assemble several steps that can be
    cross-validated together while setting different parameters. For this, it
    enables setting parameters of the various steps using their names and the
    parameter name separated by a `'__'`, as in the example below. A step's
    estimator may be replaced entirely by setting the parameter with its name
    to another estimator, or a transformer removed by setting it to
    `'passthrough'` or `None`.

    For an example use case of `Pipeline` combined with
    :class:`~sklearn.model_selection.GridSearchCV`, refer to
    :ref:`sphx_glr_auto_examples_compose_plot_compare_reduction.py`. The
    example :ref:`sphx_glr_auto_examples_compose_plot_digits_pipe.py` shows how
    to grid search on a pipeline using `'__'` as a separator in the parameter names.

    Read more in the :ref:`User Guide <pipeline>`.

    .. versionadded:: 0.5

    Parameters
    ----------
    steps : list of tuples
        List of (name of step, estimator) tuples that are to be chained in
        sequential order. To be compatible with the scikit-learn API, all steps
        must define `fit`. All non-last steps must also define `transform`. See
        :ref:`Combining Estimators <combining_estimators>` for more details.

    transform_input : list of str, default=None
        The names of the :term:`metadata` parameters that should be transformed by the
        pipeline before passing it to the step consuming it.

        This enables transforming some input arguments to ``fit`` (other than ``X``)
        to be transformed by the steps of the pipeline up to the step which requires
        them. Requirement is defined via :ref:`metadata routing <metadata_routing>`.
        For instance, this can be used to pass a validation set through the pipeline.

        You can only set this if metadata routing is enabled, which you
        can enable using ``sklearn.set_config(enable_metadata_routing=True)``.

        .. versionadded:: 1.6

    memory : str or object with the joblib.Memory interface, default=None
        Used to cache the fitted transformers of the pipeline. The last step
        will never be cached, even if it is a transformer. By default, no
        caching is performed. If a string is given, it is the path to the
        caching directory. Enabling caching triggers a clone of the transformers
        before fitting. Therefore, the transformer instance given to the
        pipeline cannot be inspected directly. Use the attribute ``named_steps``
        or ``steps`` to inspect estimators within the pipeline. Caching the
        transformers is advantageous when fitting is time consuming. See
        :ref:`sphx_glr_auto_examples_neighbors_plot_caching_nearest_neighbors.py`
        for an example on how to enable caching.

    verbose : bool, default=False
        If True, the time elapsed while fitting each step will be printed as it
        is completed.

    Attributes
    ----------
    named_steps : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.
        Read-only attribute to access any step parameter by user given name.
        Keys are step names and values are steps parameters.

    classes_ : ndarray of shape (n_classes,)
        The classes labels. Only exist if the last step of the pipeline is a
        classifier.

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying first estimator in `steps` exposes such an attribute
        when fit.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Only defined if the
        underlying estimator exposes such an attribute when fit.

        .. versionadded:: 1.0

    See Also
    --------
    make_pipeline : Convenience function for simplified pipeline construction.

    Examples
    --------
    >>> from sklearn.svm import SVC
    >>> from sklearn.preprocessing import StandardScaler
    >>> from sklearn.datasets import make_classification
    >>> from sklearn.model_selection import train_test_split
    >>> from sklearn.pipeline import Pipeline
    >>> X, y = make_classification(random_state=0)
    >>> X_train, X_test, y_train, y_test = train_test_split(X, y,
    ...                                                     random_state=0)
    >>> pipe = Pipeline([('scaler', StandardScaler()), ('svc', SVC())])
    >>> # The pipeline can be used as any other estimator
    >>> # and avoids leaking the test set into the train set
    >>> pipe.fit(X_train, y_train).score(X_test, y_test)
    0.88
    >>> # An estimator's parameter can be set using '__' syntax
    >>> pipe.set_params(svc__C=10).fit(X_train, y_train).score(X_test, y_test)
    0.76
    """

    # BaseEstimator interface
    _parameter_constraints: dict = {
        "steps": [list, Hidden(tuple)],
        "transform_input": [list, None],
        "memory": [None, str, HasMethods(["cache"])],
        "verbose": ["boolean"],
    }

    def __init__(self, steps, *, transform_input=None, memory=None, verbose=False):
        self.steps = steps
        self.transform_input = transform_input
        self.memory = memory
        self.verbose = verbose

    def set_output(self, *, transform=None):
        """Set the output container when `"transform"` and `"fit_transform"` are called.

        Calling `set_output` will set the output of all estimators in `steps`.

        Parameters
        ----------
        transform : {"default", "pandas", "polars"}, default=None
            Configure output of `transform` and `fit_transform`.

            - `"default"`: Default output format of a transformer
            - `"pandas"`: DataFrame output
            - `"polars"`: Polars output
            - `None`: Transform configuration is unchanged

            .. versionadded:: 1.4
                `"polars"` option was added.

        Returns
        -------
        self : estimator instance
            Estimator instance.
        """
        for _, _, step in self._iter():
            _safe_set_output(step, transform=transform)
        return self

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Returns the parameters given in the constructor as well as the
        estimators contained within the `steps` of the `Pipeline`.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        return self._get_params("steps", deep=deep)

    def set_params(self, **kwargs):
        """Set the parameters of this estimator.

        Valid parameter keys can be listed with ``get_params()``. Note that
        you can directly set the parameters of the estimators contained in
        `steps`.

        Parameters
        ----------
        **kwargs : dict
            Parameters of this estimator or parameters of estimators contained
            in `steps`. Parameters of the steps may be set using its name and
            the parameter name separated by a '__'.

        Returns
        -------
        self : object
            Pipeline class instance.
        """
        self._set_params("steps", **kwargs)
        return self

    def _validate_steps(self):
        names, estimators = zip(*self.steps)

        # validate names
        self._validate_names(names)

        # validate estimators
        transformers = estimators[:-1]
        estimator = estimators[-1]

        for t in transformers:
            if t is None or t == "passthrough":
                continue
            if not (hasattr(t, "fit") or hasattr(t, "fit_transform")) or not hasattr(
                t, "transform"
            ):
                raise TypeError(
                    "All intermediate steps should be "
                    "transformers and implement fit and transform "
                    "or be the string 'passthrough' "
                    "'%s' (type %s) doesn't" % (t, type(t))
                )

        # We allow last estimator to be None as an identity transformation
        if (
            estimator is not None
            and estimator != "passthrough"
            and not hasattr(estimator, "fit")
        ):
            raise TypeError(
                "Last step of Pipeline should implement fit "
                "or be the string 'passthrough'. "
                "'%s' (type %s) doesn't" % (estimator, type(estimator))
            )

    def _iter(self, with_final=True, filter_passthrough=True):
        """
        Generate (idx, (name, trans)) tuples from self.steps

        When filter_passthrough is True, 'passthrough' and None transformers
        are filtered out.
        """
        stop = len(self.steps)
        if not with_final:
            stop -= 1

        for idx, (name, trans) in enumerate(islice(self.steps, 0, stop)):
            if not filter_passthrough:
                yield idx, name, trans
            elif trans is not None and trans != "passthrough":
                yield idx, name, trans

    def __len__(self):
        """
        Returns the length of the Pipeline
        """
        return len(self.steps)

    def __getitem__(self, ind):
        """Returns a sub-pipeline or a single estimator in the pipeline

        Indexing with an integer will return an estimator; using a slice
        returns another Pipeline instance which copies a slice of this
        Pipeline. This copy is shallow: modifying (or fitting) estimators in
        the sub-pipeline will affect the larger pipeline and vice-versa.
        However, replacing a value in `step` will not affect a copy.

        See
        :ref:`sphx_glr_auto_examples_feature_selection_plot_feature_selection_pipeline.py`
        for an example of how to use slicing to inspect part of a pipeline.
        """
        if isinstance(ind, slice):
            if ind.step not in (1, None):
                raise ValueError("Pipeline slicing only supports a step of 1")
            return self.__class__(
                self.steps[ind], memory=self.memory, verbose=self.verbose
            )
        try:
            name, est = self.steps[ind]
        except TypeError:
            # Not an int, try get step by name
            return self.named_steps[ind]
        return est

    # TODO(1.8): Remove this property
    @property
    def _estimator_type(self):
        """Return the estimator type of the last step in the pipeline."""

        if not self.steps:
            return None

        return self.steps[-1][1]._estimator_type

    @property
    def named_steps(self):
        """Access the steps by name.

        Read-only attribute to access any step by given name.
        Keys are steps names and values are the steps objects."""
        # Use Bunch object to improve autocomplete
        return Bunch(**dict(self.steps))

    @property
    def _final_estimator(self):
        try:
            estimator = self.steps[-1][1]
            return "passthrough" if estimator is None else estimator
        except (ValueError, AttributeError, TypeError):
            # This condition happens when a call to a method is first calling
            # `_available_if` and `fit` did not validate `steps` yet. We
            # return `None` and an `InvalidParameterError` will be raised
            # right after.
            return None

    def _log_message(self, step_idx):
        if not self.verbose:
            return None
        name, _ = self.steps[step_idx]

        return "(step %d of %d) Processing %s" % (step_idx + 1, len(self.steps), name)

    def _check_method_params(self, method, props, **kwargs):
        if _routing_enabled():
            routed_params = process_routing(self, method, **props, **kwargs)
            return routed_params
        else:
            fit_params_steps = Bunch(
                **{
                    name: Bunch(**{method: {} for method in METHODS})
                    for name, step in self.steps
                    if step is not None
                }
            )
            for pname, pval in props.items():
                if "__" not in pname:
                    raise ValueError(
                        "Pipeline.fit does not accept the {} parameter. "
                        "You can pass parameters to specific steps of your "
                        "pipeline using the stepname__parameter format, e.g. "
                        "`Pipeline.fit(X, y, logisticregression__sample_weight"
                        "=sample_weight)`.".format(pname)
                    )
                step, param = pname.split("__", 1)
                fit_params_steps[step]["fit"][param] = pval
                # without metadata routing, fit_transform and fit_predict
                # get all the same params and pass it to the last fit.
                fit_params_steps[step]["fit_transform"][param] = pval
                fit_params_steps[step]["fit_predict"][param] = pval
            return fit_params_steps

    def _get_metadata_for_step(self, *, step_idx, step_params, all_params):
        """Get params (metadata) for step `name`.

        This transforms the metadata up to this step if required, which is
        indicated by the `transform_input` parameter.

        If a param in `step_params` is included in the `transform_input` list,
        it will be transformed.

        Parameters
        ----------
        step_idx : int
            Index of the step in the pipeline.

        step_params : dict
            Parameters specific to the step. These are routed parameters, e.g.
            `routed_params[name]`. If a parameter name here is included in the
            `pipeline.transform_input`, then it will be transformed. Note that
            these parameters are *after* routing, so the aliases are already
            resolved.

        all_params : dict
            All parameters passed by the user. Here this is used to call
            `transform` on the slice of the pipeline itself.

        Returns
        -------
        dict
            Parameters to be passed to the step. The ones which should be
            transformed are transformed.
        """
        if (
            self.transform_input is None
            or not all_params
            or not step_params
            or step_idx == 0
        ):
            # we only need to process step_params if transform_input is set
            # and metadata is given by the user.
            return step_params

        sub_pipeline = self[:step_idx]
        sub_metadata_routing = get_routing_for_object(sub_pipeline)
        # here we get the metadata required by sub_pipeline.transform
        transform_params = {
            key: value
            for key, value in all_params.items()
            if key
            in sub_metadata_routing.consumes(
                method="transform", params=all_params.keys()
            )
        }
        transformed_params = dict()  # this is to be returned
        transformed_cache = dict()  # used to transform each param once
        # `step_params` is the output of `process_routing`, so it has a dict for each
        # method (e.g. fit, transform, predict), which are the args to be passed to
        # those methods. We need to transform the parameters which are in the
        # `transform_input`, before returning these dicts.
        for method, method_params in step_params.items():
            transformed_params[method] = Bunch()
            for param_name, param_value in method_params.items():
                # An example of `(param_name, param_value)` is
                # `('sample_weight', array([0.5, 0.5, ...]))`
                if param_name in self.transform_input:
                    # This parameter now needs to be transformed by the sub_pipeline, to
                    # this step. We cache these computations to avoid repeating them.
                    transformed_params[method][param_name] = _cached_transform(
                        sub_pipeline,
                        cache=transformed_cache,
                        param_name=param_name,
                        param_value=param_value,
                        transform_params=transform_params,
                    )
                else:
                    transformed_params[method][param_name] = param_value
        return transformed_params

    # Estimator interface

    def _fit(self, X, y=None, routed_params=None, raw_params=None):
        """Fit the pipeline except the last step.

        routed_params is the output of `process_routing`
        raw_params is the parameters passed by the user, used when `transform_input`
            is set by the user, to transform metadata using a sub-pipeline.
        """
        # shallow copy of steps - this should really be steps_
        self.steps = list(self.steps)
        self._validate_steps()
        # Setup the memory
        memory = check_memory(self.memory)

        fit_transform_one_cached = memory.cache(_fit_transform_one)

        for step_idx, name, transformer in self._iter(
            with_final=False, filter_passthrough=False
        ):
            if transformer is None or transformer == "passthrough":
                with _print_elapsed_time("Pipeline", self._log_message(step_idx)):
                    continue

            if hasattr(memory, "location") and memory.location is None:
                # we do not clone when caching is disabled to
                # preserve backward compatibility
                cloned_transformer = transformer
            else:
                cloned_transformer = clone(transformer)
            # Fit or load from cache the current transformer
            step_params = self._get_metadata_for_step(
                step_idx=step_idx,
                step_params=routed_params[name],
                all_params=raw_params,
            )

            X, fitted_transformer = fit_transform_one_cached(
                cloned_transformer,
                X,
                y,
                weight=None,
                message_clsname="Pipeline",
                message=self._log_message(step_idx),
                params=step_params,
            )
            # Replace the transformer of the step with the fitted
            # transformer. This is necessary when loading the transformer
            # from the cache.
            self.steps[step_idx] = (name, fitted_transformer)
        return X

    @_fit_context(
        # estimators in Pipeline.steps are not validated yet
        prefer_skip_nested_validation=False
    )
    def fit(self, X, y=None, **params):
        """Fit the model.

        Fit all the transformers one after the other and sequentially transform the
        data. Finally, fit the transformed data using the final estimator.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of the
            pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps of
            the pipeline.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters passed to the
              ``fit`` method of each step, where each parameter name is prefixed such
              that parameter ``p`` for step ``s`` has key ``s__p``.

            - If `enable_metadata_routing=True`: Parameters requested and accepted by
              steps. Each step must have requested certain metadata for these parameters
              to be forwarded to them.

            .. versionchanged:: 1.4
                Parameters are now passed to the ``transform`` method of the
                intermediate steps as well, if requested, and if
                `enable_metadata_routing=True` is set via
                :func:`~sklearn.set_config`.

            See :ref:`Metadata Routing User Guide <metadata_routing>` for more
            details.

        Returns
        -------
        self : object
            Pipeline with fitted steps.
        """
        if not _routing_enabled() and self.transform_input is not None:
            raise ValueError(
                "The `transform_input` parameter can only be set if metadata "
                "routing is enabled. You can enable metadata routing using "
                "`sklearn.set_config(enable_metadata_routing=True)`."
            )

        routed_params = self._check_method_params(method="fit", props=params)
        Xt = self._fit(X, y, routed_params, raw_params=params)
        with _print_elapsed_time("Pipeline", self._log_message(len(self.steps) - 1)):
            if self._final_estimator != "passthrough":
                last_step_params = self._get_metadata_for_step(
                    step_idx=len(self) - 1,
                    step_params=routed_params[self.steps[-1][0]],
                    all_params=params,
                )
                self._final_estimator.fit(Xt, y, **last_step_params["fit"])

        return self

    def _can_fit_transform(self):
        return (
            self._final_estimator == "passthrough"
            or hasattr(self._final_estimator, "transform")
            or hasattr(self._final_estimator, "fit_transform")
        )

    @available_if(_can_fit_transform)
    @_fit_context(
        # estimators in Pipeline.steps are not validated yet
        prefer_skip_nested_validation=False
    )
    def fit_transform(self, X, y=None, **params):
        """Fit the model and transform with the final estimator.

        Fit all the transformers one after the other and sequentially transform
        the data. Only valid if the final estimator either implements
        `fit_transform` or `fit` and `transform`.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of the
            pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps of
            the pipeline.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters passed to the
              ``fit`` method of each step, where each parameter name is prefixed such
              that parameter ``p`` for step ``s`` has key ``s__p``.

            - If `enable_metadata_routing=True`: Parameters requested and accepted by
              steps. Each step must have requested certain metadata for these parameters
              to be forwarded to them.

            .. versionchanged:: 1.4
                Parameters are now passed to the ``transform`` method of the
                intermediate steps as well, if requested, and if
                `enable_metadata_routing=True`.

            See :ref:`Metadata Routing User Guide <metadata_routing>` for more
            details.

        Returns
        -------
        Xt : ndarray of shape (n_samples, n_transformed_features)
            Transformed samples.
        """
        routed_params = self._check_method_params(method="fit_transform", props=params)
        Xt = self._fit(X, y, routed_params)

        last_step = self._final_estimator
        with _print_elapsed_time("Pipeline", self._log_message(len(self.steps) - 1)):
            if last_step == "passthrough":
                return Xt
            last_step_params = self._get_metadata_for_step(
                step_idx=len(self) - 1,
                step_params=routed_params[self.steps[-1][0]],
                all_params=params,
            )
            if hasattr(last_step, "fit_transform"):
                return last_step.fit_transform(
                    Xt, y, **last_step_params["fit_transform"]
                )
            else:
                return last_step.fit(Xt, y, **last_step_params["fit"]).transform(
                    Xt, **last_step_params["transform"]
                )

    @available_if(_final_estimator_has("predict"))
    def predict(self, X, **params):
        """Transform the data, and apply `predict` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls `predict`
        method. Only valid if the final estimator implements `predict`.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters to the
              ``predict`` called at the end of all transformations in the pipeline.

            - If `enable_metadata_routing=True`: Parameters requested and accepted by
              steps. Each step must have requested certain metadata for these parameters
              to be forwarded to them.

            .. versionadded:: 0.20

            .. versionchanged:: 1.4
                Parameters are now passed to the ``transform`` method of the
                intermediate steps as well, if requested, and if
                `enable_metadata_routing=True` is set via
                :func:`~sklearn.set_config`.

            See :ref:`Metadata Routing User Guide <metadata_routing>` for more
            details.

            Note that while this may be used to return uncertainties from some
            models with ``return_std`` or ``return_cov``, uncertainties that are
            generated by the transformations in the pipeline are not propagated
            to the final estimator.

        Returns
        -------
        y_pred : ndarray
            Result of calling `predict` on the final estimator.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            Xt = X

            if not _routing_enabled():
                for _, name, transform in self._iter(with_final=False):
                    Xt = transform.transform(Xt)
                return self.steps[-1][1].predict(Xt, **params)

            # metadata routing enabled
            routed_params = process_routing(self, "predict", **params)
            for _, name, transform in self._iter(with_final=False):
                Xt = transform.transform(Xt, **routed_params[name].transform)
            return self.steps[-1][1].predict(
                Xt, **routed_params[self.steps[-1][0]].predict
            )

    @available_if(_final_estimator_has("fit_predict"))
    @_fit_context(
        # estimators in Pipeline.steps are not validated yet
        prefer_skip_nested_validation=False
    )
    def fit_predict(self, X, y=None, **params):
        """Transform the data, and apply `fit_predict` with the final estimator.

        Call `fit_transform` of each transformer in the pipeline. The
        transformed data are finally passed to the final estimator that calls
        `fit_predict` method. Only valid if the final estimator implements
        `fit_predict`.

        Parameters
        ----------
        X : iterable
            Training data. Must fulfill input requirements of first step of
            the pipeline.

        y : iterable, default=None
            Training targets. Must fulfill label requirements for all steps
            of the pipeline.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters to the
              ``predict`` called at the end of all transformations in the pipeline.

            - If `enable_metadata_routing=True`: Parameters requested and accepted by
              steps. Each step must have requested certain metadata for these parameters
              to be forwarded to them.

            .. versionadded:: 0.20

            .. versionchanged:: 1.4
                Parameters are now passed to the ``transform`` method of the
                intermediate steps as well, if requested, and if
                `enable_metadata_routing=True`.

            See :ref:`Metadata Routing User Guide <metadata_routing>` for more
            details.

            Note that while this may be used to return uncertainties from some
            models with ``return_std`` or ``return_cov``, uncertainties that are
            generated by the transformations in the pipeline are not propagated
            to the final estimator.

        Returns
        -------
        y_pred : ndarray
            Result of calling `fit_predict` on the final estimator.
        """
        routed_params = self._check_method_params(method="fit_predict", props=params)
        Xt = self._fit(X, y, routed_params)

        params_last_step = routed_params[self.steps[-1][0]]
        with _print_elapsed_time("Pipeline", self._log_message(len(self.steps) - 1)):
            y_pred = self.steps[-1][1].fit_predict(
                Xt, y, **params_last_step.get("fit_predict", {})
            )
        return y_pred

    @available_if(_final_estimator_has("predict_proba"))
    def predict_proba(self, X, **params):
        """Transform the data, and apply `predict_proba` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls
        `predict_proba` method. Only valid if the final estimator implements
        `predict_proba`.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters to the
              `predict_proba` called at the end of all transformations in the pipeline.

            - If `enable_metadata_routing=True`: Parameters requested and accepted by
              steps. Each step must have requested certain metadata for these parameters
              to be forwarded to them.

            .. versionadded:: 0.20

            .. versionchanged:: 1.4
                Parameters are now passed to the ``transform`` method of the
                intermediate steps as well, if requested, and if
                `enable_metadata_routing=True`.

            See :ref:`Metadata Routing User Guide <metadata_routing>` for more
            details.

        Returns
        -------
        y_proba : ndarray of shape (n_samples, n_classes)
            Result of calling `predict_proba` on the final estimator.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            Xt = X

            if not _routing_enabled():
                for _, name, transform in self._iter(with_final=False):
                    Xt = transform.transform(Xt)
                return self.steps[-1][1].predict_proba(Xt, **params)

            # metadata routing enabled
            routed_params = process_routing(self, "predict_proba", **params)
            for _, name, transform in self._iter(with_final=False):
                Xt = transform.transform(Xt, **routed_params[name].transform)
            return self.steps[-1][1].predict_proba(
                Xt, **routed_params[self.steps[-1][0]].predict_proba
            )

    @available_if(_final_estimator_has("decision_function"))
    def decision_function(self, X, **params):
        """Transform the data, and apply `decision_function` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls
        `decision_function` method. Only valid if the final estimator
        implements `decision_function`.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        **params : dict of string -> object
            Parameters requested and accepted by steps. Each step must have
            requested certain metadata for these parameters to be forwarded to
            them.

            .. versionadded:: 1.4
                Only available if `enable_metadata_routing=True`. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        y_score : ndarray of shape (n_samples, n_classes)
            Result of calling `decision_function` on the final estimator.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            _raise_for_params(params, self, "decision_function")

            # not branching here since params is only available if
            # enable_metadata_routing=True
            routed_params = process_routing(self, "decision_function", **params)

            Xt = X
            for _, name, transform in self._iter(with_final=False):
                Xt = transform.transform(
                    Xt, **routed_params.get(name, {}).get("transform", {})
                )
            return self.steps[-1][1].decision_function(
                Xt,
                **routed_params.get(self.steps[-1][0], {}).get("decision_function", {}),
            )

    @available_if(_final_estimator_has("score_samples"))
    def score_samples(self, X):
        """Transform the data, and apply `score_samples` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls
        `score_samples` method. Only valid if the final estimator implements
        `score_samples`.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        Returns
        -------
        y_score : ndarray of shape (n_samples,)
            Result of calling `score_samples` on the final estimator.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            Xt = X
            for _, _, transformer in self._iter(with_final=False):
                Xt = transformer.transform(Xt)
            return self.steps[-1][1].score_samples(Xt)

    @available_if(_final_estimator_has("predict_log_proba"))
    def predict_log_proba(self, X, **params):
        """Transform the data, and apply `predict_log_proba` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls
        `predict_log_proba` method. Only valid if the final estimator
        implements `predict_log_proba`.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        **params : dict of str -> object
            - If `enable_metadata_routing=False` (default): Parameters to the
              `predict_log_proba` called at the end of all transformations in the
              pipeline.

            - If `enable_metadata_routing=True`: Parameters requested and accepted by
              steps. Each step must have requested certain metadata for these parameters
              to be forwarded to them.

            .. versionadded:: 0.20

            .. versionchanged:: 1.4
                Parameters are now passed to the ``transform`` method of the
                intermediate steps as well, if requested, and if
                `enable_metadata_routing=True`.

            See :ref:`Metadata Routing User Guide <metadata_routing>` for more
            details.

        Returns
        -------
        y_log_proba : ndarray of shape (n_samples, n_classes)
            Result of calling `predict_log_proba` on the final estimator.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            Xt = X

            if not _routing_enabled():
                for _, name, transform in self._iter(with_final=False):
                    Xt = transform.transform(Xt)
                return self.steps[-1][1].predict_log_proba(Xt, **params)

            # metadata routing enabled
            routed_params = process_routing(self, "predict_log_proba", **params)
            for _, name, transform in self._iter(with_final=False):
                Xt = transform.transform(Xt, **routed_params[name].transform)
            return self.steps[-1][1].predict_log_proba(
                Xt, **routed_params[self.steps[-1][0]].predict_log_proba
            )

    def _can_transform(self):
        return self._final_estimator == "passthrough" or hasattr(
            self._final_estimator, "transform"
        )

    @available_if(_can_transform)
    def transform(self, X, **params):
        """Transform the data, and apply `transform` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls
        `transform` method. Only valid if the final estimator
        implements `transform`.

        This also works where final estimator is `None` in which case all prior
        transformations are applied.

        Parameters
        ----------
        X : iterable
            Data to transform. Must fulfill input requirements of first step
            of the pipeline.

        **params : dict of str -> object
            Parameters requested and accepted by steps. Each step must have
            requested certain metadata for these parameters to be forwarded to
            them.

            .. versionadded:: 1.4
                Only available if `enable_metadata_routing=True`. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        Xt : ndarray of shape (n_samples, n_transformed_features)
            Transformed data.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            _raise_for_params(params, self, "transform")

            # not branching here since params is only available if
            # enable_metadata_routing=True
            routed_params = process_routing(self, "transform", **params)
            Xt = X
            for _, name, transform in self._iter():
                Xt = transform.transform(Xt, **routed_params[name].transform)
            return Xt

    def _can_inverse_transform(self):
        return all(hasattr(t, "inverse_transform") for _, _, t in self._iter())

    @available_if(_can_inverse_transform)
    def inverse_transform(self, X, **params):
        """Apply `inverse_transform` for each step in a reverse order.

        All estimators in the pipeline must support `inverse_transform`.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_transformed_features)
            Data samples, where ``n_samples`` is the number of samples and
            ``n_features`` is the number of features. Must fulfill
            input requirements of last step of pipeline's
            ``inverse_transform`` method.

        **params : dict of str -> object
            Parameters requested and accepted by steps. Each step must have
            requested certain metadata for these parameters to be forwarded to
            them.

            .. versionadded:: 1.4
                Only available if `enable_metadata_routing=True`. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Inverse transformed data, that is, data in the original feature
            space.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            _raise_for_params(params, self, "inverse_transform")

            # we don't have to branch here, since params is only non-empty if
            # enable_metadata_routing=True.
            routed_params = process_routing(self, "inverse_transform", **params)
            reverse_iter = reversed(list(self._iter()))
            for _, name, transform in reverse_iter:
                X = transform.inverse_transform(
                    X, **routed_params[name].inverse_transform
                )
            return X

    @available_if(_final_estimator_has("score"))
    def score(self, X, y=None, sample_weight=None, **params):
        """Transform the data, and apply `score` with the final estimator.

        Call `transform` of each transformer in the pipeline. The transformed
        data are finally passed to the final estimator that calls
        `score` method. Only valid if the final estimator implements `score`.

        Parameters
        ----------
        X : iterable
            Data to predict on. Must fulfill input requirements of first step
            of the pipeline.

        y : iterable, default=None
            Targets used for scoring. Must fulfill label requirements for all
            steps of the pipeline.

        sample_weight : array-like, default=None
            If not None, this argument is passed as ``sample_weight`` keyword
            argument to the ``score`` method of the final estimator.

        **params : dict of str -> object
            Parameters requested and accepted by steps. Each step must have
            requested certain metadata for these parameters to be forwarded to
            them.

            .. versionadded:: 1.4
                Only available if `enable_metadata_routing=True`. See
                :ref:`Metadata Routing User Guide <metadata_routing>` for more
                details.

        Returns
        -------
        score : float
            Result of calling `score` on the final estimator.
        """
        # TODO(1.8): Remove the context manager and use check_is_fitted(self)
        with _raise_or_warn_if_not_fitted(self):
            Xt = X
            if not _routing_enabled():
                for _, name, transform in self._iter(with_final=False):
                    Xt = transform.transform(Xt)
                score_params = {}
                if sample_weight is not None:
                    score_params["sample_weight"] = sample_weight
                return self.steps[-1][1].score(Xt, y, **score_params)

            # metadata routing is enabled.
            routed_params = process_routing(
                self, "score", sample_weight=sample_weight, **params
            )

            Xt = X
            for _, name, transform in self._iter(with_final=False):
                Xt = transform.transform(Xt, **routed_params[name].transform)
            return self.steps[-1][1].score(
                Xt, y, **routed_params[self.steps[-1][0]].score
            )

    @property
    def classes_(self):
        """The classes labels. Only exist if the last step is a classifier."""
        return self.steps[-1][1].classes_

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()

        if not self.steps:
            return tags

        try:
            if self.steps[0][1] is not None and self.steps[0][1] != "passthrough":
                tags.input_tags.pairwise = get_tags(
                    self.steps[0][1]
                ).input_tags.pairwise
            # WARNING: the sparse tag can be incorrect.
            # Some Pipelines accepting sparse data are wrongly tagged sparse=False.
            # For example Pipeline([PCA(), estimator]) accepts sparse data
            # even if the estimator doesn't as PCA outputs a dense array.
            tags.input_tags.sparse = all(
                get_tags(step).input_tags.sparse
                for name, step in self.steps
                if step is not None and step != "passthrough"
            )
        except (ValueError, AttributeError, TypeError):
            # This happens when the `steps` is not a list of (name, estimator)
            # tuples and `fit` is not called yet to validate the steps.
            pass

        try:
            if self.steps[-1][1] is not None and self.steps[-1][1] != "passthrough":
                last_step_tags = get_tags(self.steps[-1][1])
                tags.estimator_type = last_step_tags.estimator_type
                tags.target_tags.multi_output = last_step_tags.target_tags.multi_output
                tags.classifier_tags = deepcopy(last_step_tags.classifier_tags)
                tags.regressor_tags = deepcopy(last_step_tags.regressor_tags)
                tags.transformer_tags = deepcopy(last_step_tags.transformer_tags)
        except (ValueError, AttributeError, TypeError):
            # This happens when the `steps` is not a list of (name, estimator)
            # tuples and `fit` is not called yet to validate the steps.
            pass

        return tags

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Transform input features using the pipeline.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Input features.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names.
        """
        feature_names_out = input_features
        for _, name, transform in self._iter():
            if not hasattr(transform, "get_feature_names_out"):
                raise AttributeError(
                    "Estimator {} does not provide get_feature_names_out. "
                    "Did you mean to call pipeline[:-1].get_feature_names_out"
                    "()?".format(name)
                )
            feature_names_out = transform.get_feature_names_out(feature_names_out)
        return feature_names_out

    @property
    def n_features_in_(self):
        """Number of features seen during first step `fit` method."""
        # delegate to first step (which will call check_is_fitted)
        return self.steps[0][1].n_features_in_

    @property
    def feature_names_in_(self):
        """Names of features seen during first step `fit` method."""
        # delegate to first step (which will call check_is_fitted)
        return self.steps[0][1].feature_names_in_

    def __sklearn_is_fitted__(self):
        """Indicate whether pipeline has been fit.

        This is done by checking whether the last non-`passthrough` step of the
        pipeline is fitted.

        An empty pipeline is considered fitted.
        """

        # First find the last step that is not 'passthrough'
        last_step = None
        for _, estimator in reversed(self.steps):
            if estimator != "passthrough":
                last_step = estimator
                break

        if last_step is None:
            # All steps are 'passthrough', so the pipeline is considered fitted
            return True

        try:
            # check if the last step of the pipeline is fitted
            # we only check the last step since if the last step is fit, it
            # means the previous steps should also be fit. This is faster than
            # checking if every step of the pipeline is fit.
            check_is_fitted(last_step)
            return True
        except NotFittedError:
            return False

    def _sk_visual_block_(self):
        _, estimators = zip(*self.steps)

        def _get_name(name, est):
            if est is None or est == "passthrough":
                return f"{name}: passthrough"
            # Is an estimator
            return f"{name}: {est.__class__.__name__}"

        names = [_get_name(name, est) for name, est in self.steps]
        name_details = [str(est) for est in estimators]
        return _VisualBlock(
            "serial",
            estimators,
            names=names,
            name_details=name_details,
            dash_wrapped=False,
        )

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self.__class__.__name__)

        # first we add all steps except the last one
        for _, name, trans in self._iter(with_final=False, filter_passthrough=True):
            method_mapping = MethodMapping()
            # fit, fit_predict, and fit_transform call fit_transform if it
            # exists, or else fit and transform
            if hasattr(trans, "fit_transform"):
                (
                    method_mapping.add(caller="fit", callee="fit_transform")
                    .add(caller="fit_transform", callee="fit_transform")
                    .add(caller="fit_predict", callee="fit_transform")
                )
            else:
                (
                    method_mapping.add(caller="fit", callee="fit")
                    .add(caller="fit", callee="transform")
                    .add(caller="fit_transform", callee="fit")
                    .add(caller="fit_transform", callee="transform")
                    .add(caller="fit_predict", callee="fit")
                    .add(caller="fit_predict", callee="transform")
                )

            (
                method_mapping.add(caller="predict", callee="transform")
                .add(caller="predict", callee="transform")
                .add(caller="predict_proba", callee="transform")
                .add(caller="decision_function", callee="transform")
                .add(caller="predict_log_proba", callee="transform")
                .add(caller="transform", callee="transform")
                .add(caller="inverse_transform", callee="inverse_transform")
                .add(caller="score", callee="transform")
            )

            router.add(method_mapping=method_mapping, **{name: trans})

        final_name, final_est = self.steps[-1]
        if final_est is None or final_est == "passthrough":
            return router

        # then we add the last step
        method_mapping = MethodMapping()
        if hasattr(final_est, "fit_transform"):
            method_mapping.add(caller="fit_transform", callee="fit_transform")
        else:
            method_mapping.add(caller="fit", callee="fit").add(
                caller="fit", callee="transform"
            )
        (
            method_mapping.add(caller="fit", callee="fit")
            .add(caller="predict", callee="predict")
            .add(caller="fit_predict", callee="fit_predict")
            .add(caller="predict_proba", callee="predict_proba")
            .add(caller="decision_function", callee="decision_function")
            .add(caller="predict_log_proba", callee="predict_log_proba")
            .add(caller="transform", callee="transform")
            .add(caller="inverse_transform", callee="inverse_transform")
            .add(caller="score", callee="score")
        )

        router.add(method_mapping=method_mapping, **{final_name: final_est})
        return router


def _name_estimators(estimators):
    """Generate names for estimators."""

    names = [
        estimator if isinstance(estimator, str) else type(estimator).__name__.lower()
        for estimator in estimators
    ]
    namecount = defaultdict(int)
    for est, name in zip(estimators, names):
        namecount[name] += 1

    for k, v in list(namecount.items()):
        if v == 1:
            del namecount[k]

    for i in reversed(range(len(estimators))):
        name = names[i]
        if name in namecount:
            names[i] += "-%d" % namecount[name]
            namecount[name] -= 1

    return list(zip(names, estimators))


def make_pipeline(*steps, memory=None, transform_input=None, verbose=False):
    """Construct a :class:`Pipeline` from the given estimators.

    This is a shorthand for the :class:`Pipeline` constructor; it does not
    require, and does not permit, naming the estimators. Instead, their names
    will be set to the lowercase of their types automatically.

    Parameters
    ----------
    *steps : list of Estimator objects
        List of the scikit-learn estimators that are chained together.

    memory : str or object with the joblib.Memory interface, default=None
        Used to cache the fitted transformers of the pipeline. The last step
        will never be cached, even if it is a transformer. By default, no
        caching is performed. If a string is given, it is the path to the
        caching directory. Enabling caching triggers a clone of the transformers
        before fitting. Therefore, the transformer instance given to the
        pipeline cannot be inspected directly. Use the attribute ``named_steps``
        or ``steps`` to inspect estimators within the pipeline. Caching the
        transformers is advantageous when fitting is time consuming.

    transform_input : list of str, default=None
        This enables transforming some input arguments to ``fit`` (other than ``X``)
        to be transformed by the steps of the pipeline up to the step which requires
        them. Requirement is defined via :ref:`metadata routing <metadata_routing>`.
        This can be used to pass a validation set through the pipeline for instance.

        You can only set this if metadata routing is enabled, which you
        can enable using ``sklearn.set_config(enable_metadata_routing=True)``.

        .. versionadded:: 1.6

    verbose : bool, default=False
        If True, the time elapsed while fitting each step will be printed as it
        is completed.

    Returns
    -------
    p : Pipeline
        Returns a scikit-learn :class:`Pipeline` object.

    See Also
    --------
    Pipeline : Class for creating a pipeline of transforms with a final
        estimator.

    Examples
    --------
    >>> from sklearn.naive_bayes import GaussianNB
    >>> from sklearn.preprocessing import StandardScaler
    >>> from sklearn.pipeline import make_pipeline
    >>> make_pipeline(StandardScaler(), GaussianNB(priors=None))
    Pipeline(steps=[('standardscaler', StandardScaler()),
                    ('gaussiannb', GaussianNB())])
    """
    return Pipeline(
        _name_estimators(steps),
        transform_input=transform_input,
        memory=memory,
        verbose=verbose,
    )


def _transform_one(transformer, X, y, weight, params):
    """Call transform and apply weight to output.

    Parameters
    ----------
    transformer : estimator
        Estimator to be used for transformation.

    X : {array-like, sparse matrix} of shape (n_samples, n_features)
        Input data to be transformed.

    y : ndarray of shape (n_samples,)
        Ignored.

    weight : float
        Weight to be applied to the output of the transformation.

    params : dict
        Parameters to be passed to the transformer's ``transform`` method.

        This should be of the form ``process_routing()["step_name"]``.
    """
    res = transformer.transform(X, **params.transform)
    # if we have a weight for this transformer, multiply output
    if weight is None:
        return res
    return res * weight


def _fit_transform_one(
    transformer, X, y, weight, message_clsname="", message=None, params=None
):
    """
    Fits ``transformer`` to ``X`` and ``y``. The transformed result is returned
    with the fitted transformer. If ``weight`` is not ``None``, the result will
    be multiplied by ``weight``.

    ``params`` needs to be of the form ``process_routing()["step_name"]``.
    """
    params = params or {}
    with _print_elapsed_time(message_clsname, message):
        if hasattr(transformer, "fit_transform"):
            res = transformer.fit_transform(X, y, **params.get("fit_transform", {}))
        else:
            res = transformer.fit(X, y, **params.get("fit", {})).transform(
                X, **params.get("transform", {})
            )

    if weight is None:
        return res, transformer
    return res * weight, transformer


def _fit_one(transformer, X, y, weight, message_clsname="", message=None, params=None):
    """
    Fits ``transformer`` to ``X`` and ``y``.
    """
    with _print_elapsed_time(message_clsname, message):
        return transformer.fit(X, y, **params["fit"])


class FeatureUnion(TransformerMixin, _BaseComposition):
    """Concatenates results of multiple transformer objects.

    This estimator applies a list of transformer objects in parallel to the
    input data, then concatenates the results. This is useful to combine
    several feature extraction mechanisms into a single transformer.

    Parameters of the transformers may be set using its name and the parameter
    name separated by a '__'. A transformer may be replaced entirely by
    setting the parameter with its name to another transformer, removed by
    setting to 'drop' or disabled by setting to 'passthrough' (features are
    passed without transformation).

    Read more in the :ref:`User Guide <feature_union>`.

    .. versionadded:: 0.13

    Parameters
    ----------
    transformer_list : list of (str, transformer) tuples
        List of transformer objects to be applied to the data. The first
        half of each tuple is the name of the transformer. The transformer can
        be 'drop' for it to be ignored or can be 'passthrough' for features to
        be passed unchanged.

        .. versionadded:: 1.1
           Added the option `"passthrough"`.

        .. versionchanged:: 0.22
           Deprecated `None` as a transformer in favor of 'drop'.

    n_jobs : int, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

        .. versionchanged:: v0.20
           `n_jobs` default changed from 1 to None

    transformer_weights : dict, default=None
        Multiplicative weights for features per transformer.
        Keys are transformer names, values the weights.
        Raises ValueError if key not present in ``transformer_list``.

    verbose : bool, default=False
        If True, the time elapsed while fitting each transformer will be
        printed as it is completed.

    verbose_feature_names_out : bool, default=True
        If True, :meth:`get_feature_names_out` will prefix all feature names
        with the name of the transformer that generated that feature.
        If False, :meth:`get_feature_names_out` will not prefix any feature
        names and will error if feature names are not unique.

        .. versionadded:: 1.5

    Attributes
    ----------
    named_transformers : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.
        Read-only attribute to access any transformer parameter by user
        given name. Keys are transformer names and values are
        transformer parameters.

        .. versionadded:: 1.2

    n_features_in_ : int
        Number of features seen during :term:`fit`. Only defined if the
        underlying first transformer in `transformer_list` exposes such an
        attribute when fit.

        .. versionadded:: 0.24

    feature_names_in_ : ndarray of shape (`n_features_in_`,)
        Names of features seen during :term:`fit`. Defined only when
        `X` has feature names that are all strings.

        .. versionadded:: 1.3

    See Also
    --------
    make_union : Convenience function for simplified feature union
        construction.

    Examples
    --------
    >>> from sklearn.pipeline import FeatureUnion
    >>> from sklearn.decomposition import PCA, TruncatedSVD
    >>> union = FeatureUnion([("pca", PCA(n_components=1)),
    ...                       ("svd", TruncatedSVD(n_components=2))])
    >>> X = [[0., 1., 3], [2., 2., 5]]
    >>> union.fit_transform(X)
    array([[-1.5       ,  3.04, -0.872],
           [ 1.5       ,  5.72,  0.463]])
    >>> # An estimator's parameter can be set using '__' syntax
    >>> union.set_params(svd__n_components=1).fit_transform(X)
    array([[-1.5       ,  3.04],
           [ 1.5       ,  5.72]])

    For a more detailed example of usage, see
    :ref:`sphx_glr_auto_examples_compose_plot_feature_union.py`.
    """

    def __init__(
        self,
        transformer_list,
        *,
        n_jobs=None,
        transformer_weights=None,
        verbose=False,
        verbose_feature_names_out=True,
    ):
        self.transformer_list = transformer_list
        self.n_jobs = n_jobs
        self.transformer_weights = transformer_weights
        self.verbose = verbose
        self.verbose_feature_names_out = verbose_feature_names_out

    def set_output(self, *, transform=None):
        """Set the output container when `"transform"` and `"fit_transform"` are called.

        `set_output` will set the output of all estimators in `transformer_list`.

        Parameters
        ----------
        transform : {"default", "pandas", "polars"}, default=None
            Configure output of `transform` and `fit_transform`.

            - `"default"`: Default output format of a transformer
            - `"pandas"`: DataFrame output
            - `"polars"`: Polars output
            - `None`: Transform configuration is unchanged

        Returns
        -------
        self : estimator instance
            Estimator instance.
        """
        super().set_output(transform=transform)
        for _, step, _ in self._iter():
            _safe_set_output(step, transform=transform)
        return self

    @property
    def named_transformers(self):
        # Use Bunch object to improve autocomplete
        return Bunch(**dict(self.transformer_list))

    def get_params(self, deep=True):
        """Get parameters for this estimator.

        Returns the parameters given in the constructor as well as the
        estimators contained within the `transformer_list` of the
        `FeatureUnion`.

        Parameters
        ----------
        deep : bool, default=True
            If True, will return the parameters for this estimator and
            contained subobjects that are estimators.

        Returns
        -------
        params : mapping of string to any
            Parameter names mapped to their values.
        """
        return self._get_params("transformer_list", deep=deep)

    def set_params(self, **kwargs):
        """Set the parameters of this estimator.

        Valid parameter keys can be listed with ``get_params()``. Note that
        you can directly set the parameters of the estimators contained in
        `transformer_list`.

        Parameters
        ----------
        **kwargs : dict
            Parameters of this estimator or parameters of estimators contained
            in `transform_list`. Parameters of the transformers may be set
            using its name and the parameter name separated by a '__'.

        Returns
        -------
        self : object
            FeatureUnion class instance.
        """
        self._set_params("transformer_list", **kwargs)
        return self

    def _validate_transformers(self):
        names, transformers = zip(*self.transformer_list)

        # validate names
        self._validate_names(names)

        # validate estimators
        for t in transformers:
            if t in ("drop", "passthrough"):
                continue
            if not (hasattr(t, "fit") or hasattr(t, "fit_transform")) or not hasattr(
                t, "transform"
            ):
                raise TypeError(
                    "All estimators should implement fit and "
                    "transform. '%s' (type %s) doesn't" % (t, type(t))
                )

    def _validate_transformer_weights(self):
        if not self.transformer_weights:
            return

        transformer_names = set(name for name, _ in self.transformer_list)
        for name in self.transformer_weights:
            if name not in transformer_names:
                raise ValueError(
                    f'Attempting to weight transformer "{name}", '
                    "but it is not present in transformer_list."
                )

    def _iter(self):
        """
        Generate (name, trans, weight) tuples excluding None and
        'drop' transformers.
        """

        get_weight = (self.transformer_weights or {}).get

        for name, trans in self.transformer_list:
            if trans == "drop":
                continue
            if trans == "passthrough":
                trans = FunctionTransformer(feature_names_out="one-to-one")
            yield (name, trans, get_weight(name))

    def get_feature_names_out(self, input_features=None):
        """Get output feature names for transformation.

        Parameters
        ----------
        input_features : array-like of str or None, default=None
            Input features.

        Returns
        -------
        feature_names_out : ndarray of str objects
            Transformed feature names.
        """
        # List of tuples (name, feature_names_out)
        transformer_with_feature_names_out = []
        for name, trans, _ in self._iter():
            if not hasattr(trans, "get_feature_names_out"):
                raise AttributeError(
                    "Transformer %s (type %s) does not provide get_feature_names_out."
                    % (str(name), type(trans).__name__)
                )
            feature_names_out = trans.get_feature_names_out(input_features)
            transformer_with_feature_names_out.append((name, feature_names_out))

        return self._add_prefix_for_feature_names_out(
            transformer_with_feature_names_out
        )

    def _add_prefix_for_feature_names_out(self, transformer_with_feature_names_out):
        """Add prefix for feature names out that includes the transformer names.

        Parameters
        ----------
        transformer_with_feature_names_out : list of tuples of (str, array-like of str)
            The tuple consistent of the transformer's name and its feature names out.

        Returns
        -------
        feature_names_out : ndarray of shape (n_features,), dtype=str
            Transformed feature names.
        """
        if self.verbose_feature_names_out:
            # Prefix the feature names out with the transformers name
            names = list(
                chain.from_iterable(
                    (f"{name}__{i}" for i in feature_names_out)
                    for name, feature_names_out in transformer_with_feature_names_out
                )
            )
            return np.asarray(names, dtype=object)

        # verbose_feature_names_out is False
        # Check that names are all unique without a prefix
        feature_names_count = Counter(
            chain.from_iterable(s for _, s in transformer_with_feature_names_out)
        )
        top_6_overlap = [
            name for name, count in feature_names_count.most_common(6) if count > 1
        ]
        top_6_overlap.sort()
        if top_6_overlap:
            if len(top_6_overlap) == 6:
                # There are more than 5 overlapping names, we only show the 5
                # of the feature names
                names_repr = str(top_6_overlap[:5])[:-1] + ", ...]"
            else:
                names_repr = str(top_6_overlap)
            raise ValueError(
                f"Output feature names: {names_repr} are not unique. Please set "
                "verbose_feature_names_out=True to add prefixes to feature names"
            )

        return np.concatenate(
            [name for _, name in transformer_with_feature_names_out],
        )

    def fit(self, X, y=None, **fit_params):
        """Fit all transformers using X.

        Parameters
        ----------
        X : iterable or array-like, depending on transformers
            Input data, used to fit transformers.

        y : array-like of shape (n_samples, n_outputs), default=None
            Targets for supervised learning.

        **fit_params : dict, default=None
            - If `enable_metadata_routing=False` (default):
              Parameters directly passed to the `fit` methods of the
              sub-transformers.

            - If `enable_metadata_routing=True`:
              Parameters safely routed to the `fit` methods of the
              sub-transformers. See :ref:`Metadata Routing User Guide
              <metadata_routing>` for more details.

            .. versionchanged:: 1.5
                `**fit_params` can be routed via metadata routing API.

        Returns
        -------
        self : object
            FeatureUnion class instance.
        """
        if _routing_enabled():
            routed_params = process_routing(self, "fit", **fit_params)
        else:
            # TODO(SLEP6): remove when metadata routing cannot be disabled.
            routed_params = Bunch()
            for name, _ in self.transformer_list:
                routed_params[name] = Bunch(fit={})
                routed_params[name].fit = fit_params

        transformers = self._parallel_func(X, y, _fit_one, routed_params)

        if not transformers:
            # All transformers are None
            return self

        self._update_transformer_list(transformers)
        return self

    def fit_transform(self, X, y=None, **params):
        """Fit all transformers, transform the data and concatenate results.

        Parameters
        ----------
        X : iterable or array-like, depending on transformers
            Input data to be transformed.

        y : array-like of shape (n_samples, n_outputs), default=None
            Targets for supervised learning.

        **params : dict, default=None
            - If `enable_metadata_routing=False` (default):
              Parameters directly passed to the `fit` methods of the
              sub-transformers.

            - If `enable_metadata_routing=True`:
              Parameters safely routed to the `fit` methods of the
              sub-transformers. See :ref:`Metadata Routing User Guide
              <metadata_routing>` for more details.

            .. versionchanged:: 1.5
                `**params` can now be routed via metadata routing API.

        Returns
        -------
        X_t : array-like or sparse matrix of \
                shape (n_samples, sum_n_components)
            The `hstack` of results of transformers. `sum_n_components` is the
            sum of `n_components` (output dimension) over transformers.
        """
        if _routing_enabled():
            routed_params = process_routing(self, "fit_transform", **params)
        else:
            # TODO(SLEP6): remove when metadata routing cannot be disabled.
            routed_params = Bunch()
            for name, obj in self.transformer_list:
                if hasattr(obj, "fit_transform"):
                    routed_params[name] = Bunch(fit_transform={})
                    routed_params[name].fit_transform = params
                else:
                    routed_params[name] = Bunch(fit={})
                    routed_params[name] = Bunch(transform={})
                    routed_params[name].fit = params

        results = self._parallel_func(X, y, _fit_transform_one, routed_params)
        if not results:
            # All transformers are None
            return np.zeros((X.shape[0], 0))

        Xs, transformers = zip(*results)
        self._update_transformer_list(transformers)

        return self._hstack(Xs)

    def _log_message(self, name, idx, total):
        if not self.verbose:
            return None
        return "(step %d of %d) Processing %s" % (idx, total, name)

    def _parallel_func(self, X, y, func, routed_params):
        """Runs func in parallel on X and y"""
        self.transformer_list = list(self.transformer_list)
        self._validate_transformers()
        self._validate_transformer_weights()
        transformers = list(self._iter())

        return Parallel(n_jobs=self.n_jobs)(
            delayed(func)(
                transformer,
                X,
                y,
                weight,
                message_clsname="FeatureUnion",
                message=self._log_message(name, idx, len(transformers)),
                params=routed_params[name],
            )
            for idx, (name, transformer, weight) in enumerate(transformers, 1)
        )

    def transform(self, X, **params):
        """Transform X separately by each transformer, concatenate results.

        Parameters
        ----------
        X : iterable or array-like, depending on transformers
            Input data to be transformed.

        **params : dict, default=None

            Parameters routed to the `transform` method of the sub-transformers via the
            metadata routing API. See :ref:`Metadata Routing User Guide
            <metadata_routing>` for more details.

            .. versionadded:: 1.5

        Returns
        -------
        X_t : array-like or sparse matrix of shape (n_samples, sum_n_components)
            The `hstack` of results of transformers. `sum_n_components` is the
            sum of `n_components` (output dimension) over transformers.
        """
        _raise_for_params(params, self, "transform")

        if _routing_enabled():
            routed_params = process_routing(self, "transform", **params)
        else:
            # TODO(SLEP6): remove when metadata routing cannot be disabled.
            routed_params = Bunch()
            for name, _ in self.transformer_list:
                routed_params[name] = Bunch(transform={})

        Xs = Parallel(n_jobs=self.n_jobs)(
            delayed(_transform_one)(trans, X, None, weight, params=routed_params[name])
            for name, trans, weight in self._iter()
        )
        if not Xs:
            # All transformers are None
            return np.zeros((X.shape[0], 0))

        return self._hstack(Xs)

    def _hstack(self, Xs):
        adapter = _get_container_adapter("transform", self)
        if adapter and all(adapter.is_supported_container(X) for X in Xs):
            return adapter.hstack(Xs)

        if any(sparse.issparse(f) for f in Xs):
            Xs = sparse.hstack(Xs).tocsr()
        else:
            Xs = np.hstack(Xs)
        return Xs

    def _update_transformer_list(self, transformers):
        transformers = iter(transformers)
        self.transformer_list[:] = [
            (name, old if old == "drop" else next(transformers))
            for name, old in self.transformer_list
        ]

    @property
    def n_features_in_(self):
        """Number of features seen during :term:`fit`."""

        # X is passed to all transformers so we just delegate to the first one
        return self.transformer_list[0][1].n_features_in_

    @property
    def feature_names_in_(self):
        """Names of features seen during :term:`fit`."""
        # X is passed to all transformers -- delegate to the first one
        return self.transformer_list[0][1].feature_names_in_

    def __sklearn_is_fitted__(self):
        # Delegate whether feature union was fitted
        for _, transformer, _ in self._iter():
            check_is_fitted(transformer)
        return True

    def _sk_visual_block_(self):
        names, transformers = zip(*self.transformer_list)
        return _VisualBlock("parallel", transformers, names=names)

    def __getitem__(self, name):
        """Return transformer with name."""
        if not isinstance(name, str):
            raise KeyError("Only string keys are supported")
        return self.named_transformers[name]

    def get_metadata_routing(self):
        """Get metadata routing of this object.

        Please check :ref:`User Guide <metadata_routing>` on how the routing
        mechanism works.

        .. versionadded:: 1.5

        Returns
        -------
        routing : MetadataRouter
            A :class:`~sklearn.utils.metadata_routing.MetadataRouter` encapsulating
            routing information.
        """
        router = MetadataRouter(owner=self.__class__.__name__)

        for name, transformer in self.transformer_list:
            router.add(
                **{name: transformer},
                method_mapping=MethodMapping()
                .add(caller="fit", callee="fit")
                .add(caller="fit_transform", callee="fit_transform")
                .add(caller="fit_transform", callee="fit")
                .add(caller="fit_transform", callee="transform")
                .add(caller="transform", callee="transform"),
            )

        return router

    def __sklearn_tags__(self):
        tags = super().__sklearn_tags__()
        try:
            tags.input_tags.sparse = all(
                get_tags(trans).input_tags.sparse
                for name, trans in self.transformer_list
                if trans not in {"passthrough", "drop"}
            )
        except Exception:
            # If `transformer_list` does not comply with our API (list of tuples)
            # then it will fail. In this case, we assume that `sparse` is False
            # but the parameter validation will raise an error during `fit`.
            pass  # pragma: no cover
        return tags


def make_union(
    *transformers, n_jobs=None, verbose=False, verbose_feature_names_out=True
):
    """Construct a :class:`FeatureUnion` from the given transformers.

    This is a shorthand for the :class:`FeatureUnion` constructor; it does not
    require, and does not permit, naming the transformers. Instead, they will
    be given names automatically based on their types. It also does not allow
    weighting.

    Parameters
    ----------
    *transformers : list of estimators
        One or more estimators.

    n_jobs : int, default=None
        Number of jobs to run in parallel.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

        .. versionchanged:: v0.20
           `n_jobs` default changed from 1 to None.

    verbose : bool, default=False
        If True, the time elapsed while fitting each transformer will be
        printed as it is completed.

    verbose_feature_names_out : bool, default=True
        If True, the feature names generated by `get_feature_names_out` will
        include prefixes derived from the transformer names.

    Returns
    -------
    f : FeatureUnion
        A :class:`FeatureUnion` object for concatenating the results of multiple
        transformer objects.

    See Also
    --------
    FeatureUnion : Class for concatenating the results of multiple transformer
        objects.

    Examples
    --------
    >>> from sklearn.decomposition import PCA, TruncatedSVD
    >>> from sklearn.pipeline import make_union
    >>> make_union(PCA(), TruncatedSVD())
     FeatureUnion(transformer_list=[('pca', PCA()),
                                   ('truncatedsvd', TruncatedSVD())])
    """
    return FeatureUnion(
        _name_estimators(transformers),
        n_jobs=n_jobs,
        verbose=verbose,
        verbose_feature_names_out=verbose_feature_names_out,
    )
