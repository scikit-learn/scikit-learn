.. currentmodule:: sklearn

.. TODO: update doc/conftest.py once document is updated and examples run.

.. _metadata_routing:

Metadata Routing
================

.. note::
  The Metadata Routing API is experimental, and is not yet implemented for all
  estimators. Please refer to the :ref:`list of supported and unsupported
  models <metadata_routing_models>` for more information. It may change without
  the usual deprecation cycle. By default this feature is not enabled. You can
  enable it by setting the ``enable_metadata_routing`` flag to
  ``True``::

    >>> import sklearn
    >>> sklearn.set_config(enable_metadata_routing=True)

  Note that the methods and requirements introduced in this document are only
  relevant if you want to pass :term:`metadata` (e.g. ``sample_weight``) to a method.
  If you're only passing ``X`` and ``y`` and no other parameter / metadata to
  methods such as :term:`fit`, :term:`transform`, etc., then you don't need to set
  anything.

This guide demonstrates how :term:`metadata` can be routed and passed between objects in
scikit-learn. If you are developing a scikit-learn compatible estimator or
meta-estimator, you can check our related developer guide:
:ref:`sphx_glr_auto_examples_miscellaneous_plot_metadata_routing.py`.

Metadata is data that an estimator, scorer, or CV splitter takes into account if the
user explicitly passes it as a parameter. For instance, :class:`~cluster.KMeans` accepts
`sample_weight` in its `fit()` method and considers it to calculate its centroids.
`classes` are consumed by some classifiers and `groups` are used in some splitters, but
any data that is passed into an object's methods apart from X and y can be considered as
metadata. Prior to scikit-learn version 1.3, there was no single API for passing
metadata like that if these objects were used in conjunction with other objects, e.g. a
scorer accepting `sample_weight` inside a :class:`~model_selection.GridSearchCV`.

With the Metadata Routing API, we can transfer metadata to estimators, scorers, and CV
splitters using :term:`meta-estimators` (such as :class:`~pipeline.Pipeline` or
:class:`~model_selection.GridSearchCV`) or functions such as
:func:`~model_selection.cross_validate` which route data to other objects. In order to
pass metadata to a method like ``fit`` or ``score``, the object consuming the metadata,
must *request* it. This is done via `set_{method}_request()` methods, where `{method}`
is substituted by the name of the method that requests the metadata. For instance,
estimators that use the metadata in their `fit()` method would use `set_fit_request()`,
and scorers would use `set_score_request()`. These methods allow us to specify which
metadata to request, for instance `set_fit_request(sample_weight=True)`.

For grouped splitters such as :class:`~model_selection.GroupKFold`, a
``groups`` parameter is requested by default. This is best demonstrated by the
following examples.

Usage Examples
**************
Here we present a few examples to show some common use-cases. Our goal is to pass
`sample_weight` and `groups` through :func:`~model_selection.cross_validate`, which
routes the metadata to :class:`~linear_model.LogisticRegressionCV` and to a custom scorer
made with :func:`~metrics.make_scorer`, both of which *can* use the metadata in their
methods. In these examples we want to individually set whether to use the metadata
within the different :term:`consumers <consumer>`.

The examples in this section require the following imports and data::

  >>> import numpy as np
  >>> from sklearn.metrics import make_scorer, accuracy_score
  >>> from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
  >>> from sklearn.model_selection import cross_validate, GridSearchCV, GroupKFold
  >>> from sklearn.feature_selection import SelectKBest
  >>> from sklearn.pipeline import make_pipeline
  >>> n_samples, n_features = 100, 4
  >>> rng = np.random.RandomState(42)
  >>> X = rng.rand(n_samples, n_features)
  >>> y = rng.randint(0, 2, size=n_samples)
  >>> my_groups = rng.randint(0, 10, size=n_samples)
  >>> my_weights = rng.rand(n_samples)
  >>> my_other_weights = rng.rand(n_samples)

Weighted scoring and fitting
----------------------------

The splitter used internally in :class:`~linear_model.LogisticRegressionCV`,
:class:`~model_selection.GroupKFold`, requests ``groups`` by default. However, we need
to explicitly request `sample_weight` for it and for our custom scorer by specifying
`sample_weight=True` in :class:`~linear_model.LogisticRegressionCV`s `set_fit_request()`
method and in :func:`~metrics.make_scorer`s `set_score_request()` method. Both
:term:`consumers <consumer>` know how to use ``sample_weight`` in their `fit()` or
`score()` methods. We can then pass the metadata in
:func:`~model_selection.cross_validate` which will route it to any active consumers::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(sample_weight=True)
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(),
  ...     scoring=weighted_acc
  ... ).set_fit_request(sample_weight=True)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     params={"sample_weight": my_weights, "groups": my_groups},
  ...     cv=GroupKFold(),
  ...     scoring=weighted_acc,
  ... )

Note that in this example, :func:`~model_selection.cross_validate` routes ``my_weights``
to both the scorer and :class:`~linear_model.LogisticRegressionCV`.

If we would pass `sample_weight` in the params of
:func:`~model_selection.cross_validate`, but not set any object to request it,
`UnsetMetadataPassedError` would be raised, hinting to us that we need to explicitly set
where to route it. The same applies if ``params={"sample_weights": my_weights, ...}``
were passed (note the typo, i.e. ``weights`` instead of ``weight``), since
``sample_weights`` was not requested by any of its underlying objects.

Weighted scoring and unweighted fitting
---------------------------------------

When passing metadata such as ``sample_weight`` into a :term:`router`
(:term:`meta-estimators` or routing function), all ``sample_weight`` :term:`consumers
<consumer>` require weights to be either explicitly requested or explicitly not
requested (i.e. ``True`` or ``False``). Thus, to perform an unweighted fit, we need to
configure :class:`~linear_model.LogisticRegressionCV` to not request sample weights, so
that :func:`~model_selection.cross_validate` does not pass the weights along::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(sample_weight=True)
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).set_fit_request(sample_weight=False)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     params={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

If :meth:`linear_model.LogisticRegressionCV.set_fit_request` had not been called,
:func:`~model_selection.cross_validate` would raise an error because ``sample_weight``
is passed but :class:`~linear_model.LogisticRegressionCV` would not be explicitly
configured to recognize the weights.

Unweighted feature selection
----------------------------

Routing metadata is only possible if the object's method knows how to use the metadata,
which in most cases means they have it as an explicit parameter. Only then we can set
request values for metadata using `set_fit_request(sample_weight=True)`, for instance.
This makes the object a :term:`consumer <consumer>`.

Unlike :class:`~linear_model.LogisticRegressionCV`,
:class:`~feature_selection.SelectKBest` can't consume weights and therefore no request
value for ``sample_weight`` on its instance is set and ``sample_weight`` is not routed
to it::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(sample_weight=True)
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).set_fit_request(sample_weight=True)
  >>> sel = SelectKBest(k=2)
  >>> pipe = make_pipeline(sel, lr)
  >>> cv_results = cross_validate(
  ...     pipe,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     params={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Different scoring and fitting weights
-------------------------------------

Despite :func:`~metrics.make_scorer` and
:class:`~linear_model.LogisticRegressionCV` both expecting the key
``sample_weight``, we can use aliases to pass different weights to different
consumers. In this example, we pass ``scoring_weight`` to the scorer, and
``fitting_weight`` to :class:`~linear_model.LogisticRegressionCV`::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(
  ...    sample_weight="scoring_weight"
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).set_fit_request(sample_weight="fitting_weight")
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     params={
  ...         "scoring_weight": my_weights,
  ...         "fitting_weight": my_other_weights,
  ...         "groups": my_groups,
  ...     },
  ...     scoring=weighted_acc,
  ... )

API Interface
*************

A :term:`consumer` is an object (estimator, meta-estimator, scorer, splitter) which
accepts and uses some :term:`metadata` in at least one of its methods (for instance
``fit``, ``predict``, ``inverse_transform``, ``transform``, ``score``, ``split``).
Meta-estimators which only forward the metadata to other objects (child estimators,
scorers, or splitters) and don't use the metadata themselves are not consumers.
(Meta-)Estimators which route metadata to other objects are :term:`routers <router>`.
A(n) (meta-)estimator can be a :term:`consumer` and a :term:`router` at the same time.
(Meta-)Estimators and splitters expose a `set_{method}_request` method for each method
which accepts at least one metadata. For instance, if an estimator supports
``sample_weight`` in ``fit`` and ``score``, it exposes
``estimator.set_fit_request(sample_weight=value)`` and
``estimator.set_score_request(sample_weight=value)``. Here ``value`` can be:

- ``True``: method requests a ``sample_weight``. This means if the metadata is provided,
  it will be used, otherwise no error is raised.
- ``False``: method does not request a ``sample_weight``.
- ``None``: router will raise an error if ``sample_weight`` is passed. This is in almost
  all cases the default value when an object is instantiated and ensures the user sets
  the metadata requests explicitly when a metadata is passed. The only exception are
  ``Group*Fold`` splitters.
- ``"param_name"``: alias for ``sample_weight`` if we want to pass different weights to
  different consumers. If aliasing is used the meta-estimator should not forward
  ``"param_name"`` to the consumer, but ``sample_weight`` instead, because the consumer
  will expect a param called ``sample_weight``. This means the mapping between the
  metadata required by the object, e.g. ``sample_weight`` and the variable name provided
  by the user, e.g. ``my_weights`` is done at the router level, and not by the consuming
  object itself.

Metadata are requested in the same way for scorers using ``set_score_request``.

If a metadata, e.g. ``sample_weight``, is passed by the user, the metadata request for
all objects which potentially can consume ``sample_weight`` should be set by the user,
otherwise an error is raised by the router object. For example, the following code
raises an error, since it hasn't been explicitly specified whether ``sample_weight``
should be passed to the estimator's scorer or not::

    >>> param_grid = {"C": [0.1, 1]}
    >>> lr = LogisticRegression().set_fit_request(sample_weight=True)
    >>> try:
    ...     GridSearchCV(
    ...         estimator=lr, param_grid=param_grid
    ...     ).fit(X, y, sample_weight=my_weights)
    ... except ValueError as e:
    ...     print(e)
    [sample_weight] are passed but are not explicitly set as requested or not
    requested for LogisticRegression.score, which is used within GridSearchCV.fit.
    Call `LogisticRegression.set_score_request({metadata}=True/False)` for each metadata
    you want to request/ignore.

The issue can be fixed by explicitly setting the request value::

    >>> lr = LogisticRegression().set_fit_request(
    ...     sample_weight=True
    ... ).set_score_request(sample_weight=False)

At the end of the **Usage Examples** section, we disable the configuration flag for
metadata routing::

    >>> sklearn.set_config(enable_metadata_routing=False)

.. _metadata_routing_models:

Metadata Routing Support Status
*******************************
All consumers (i.e. simple estimators which only consume metadata and don't
route them) support metadata routing, meaning they can be used inside
meta-estimators which support metadata routing. However, development of support
for metadata routing for meta-estimators is in progress, and here is a list of
meta-estimators and tools which support and don't yet support metadata routing.


Meta-estimators and functions supporting metadata routing:

- :class:`sklearn.calibration.CalibratedClassifierCV`
- :class:`sklearn.compose.ColumnTransformer`
- :class:`sklearn.compose.TransformedTargetRegressor`
- :class:`sklearn.covariance.GraphicalLassoCV`
- :class:`sklearn.ensemble.StackingClassifier`
- :class:`sklearn.ensemble.StackingRegressor`
- :class:`sklearn.ensemble.VotingClassifier`
- :class:`sklearn.ensemble.VotingRegressor`
- :class:`sklearn.ensemble.BaggingClassifier`
- :class:`sklearn.ensemble.BaggingRegressor`
- :class:`sklearn.feature_selection.RFE`
- :class:`sklearn.feature_selection.RFECV`
- :class:`sklearn.feature_selection.SelectFromModel`
- :class:`sklearn.feature_selection.SequentialFeatureSelector`
- :class:`sklearn.impute.IterativeImputer`
- :class:`sklearn.linear_model.ElasticNetCV`
- :class:`sklearn.linear_model.LarsCV`
- :class:`sklearn.linear_model.LassoCV`
- :class:`sklearn.linear_model.LassoLarsCV`
- :class:`sklearn.linear_model.LogisticRegressionCV`
- :class:`sklearn.linear_model.MultiTaskElasticNetCV`
- :class:`sklearn.linear_model.MultiTaskLassoCV`
- :class:`sklearn.linear_model.OrthogonalMatchingPursuitCV`
- :class:`sklearn.linear_model.RANSACRegressor`
- :class:`sklearn.linear_model.RidgeClassifierCV`
- :class:`sklearn.linear_model.RidgeCV`
- :class:`sklearn.model_selection.GridSearchCV`
- :class:`sklearn.model_selection.HalvingGridSearchCV`
- :class:`sklearn.model_selection.HalvingRandomSearchCV`
- :class:`sklearn.model_selection.RandomizedSearchCV`
- :class:`sklearn.model_selection.permutation_test_score`
- :func:`sklearn.model_selection.cross_validate`
- :func:`sklearn.model_selection.cross_val_score`
- :func:`sklearn.model_selection.cross_val_predict`
- :class:`sklearn.model_selection.learning_curve`
- :class:`sklearn.model_selection.validation_curve`
- :class:`sklearn.multiclass.OneVsOneClassifier`
- :class:`sklearn.multiclass.OneVsRestClassifier`
- :class:`sklearn.multiclass.OutputCodeClassifier`
- :class:`sklearn.multioutput.ClassifierChain`
- :class:`sklearn.multioutput.MultiOutputClassifier`
- :class:`sklearn.multioutput.MultiOutputRegressor`
- :class:`sklearn.multioutput.RegressorChain`
- :class:`sklearn.pipeline.FeatureUnion`
- :class:`sklearn.pipeline.Pipeline`
- :class:`sklearn.semi_supervised.SelfTrainingClassifier`

Meta-estimators and tools not supporting metadata routing yet:

- :class:`sklearn.ensemble.AdaBoostClassifier`
- :class:`sklearn.ensemble.AdaBoostRegressor`
