
.. _metadata_routing:

.. currentmodule:: sklearn

.. TODO: update doc/conftest.py once document is updated and examples run.

Metadata Routing
================

.. note::
  The Metadata Routing API is experimental, and is not implemented yet for many
  estimators. It may change without the usual deprecation cycle. By default
  this feature is not enabled. You can enable this feature  by setting the
  ``enable_metadata_routing`` flag to ``True``:

    >>> import sklearn
    >>> sklearn.set_config(enable_metadata_routing=True)

This guide demonstrates how metadata such as ``sample_weight`` can be routed
and passed along to estimators, scorers, and CV splitters through
meta-estimators such as :class:`~pipeline.Pipeline` and
:class:`~model_selection.GridSearchCV`. In order to pass metadata to a method
such as ``fit`` or ``score``, the object consuming the metadata, must *request*
it. For estimators and splitters, this is done via ``set_*_request`` methods,
e.g. ``set_fit_request(...)``, and for scorers this is done via the
``set_score_request`` method. For grouped splitters such as
:class:`~model_selection.GroupKFold`, a ``groups`` parameter is requested by
default. This is best demonstrated by the following examples.

If you are developing a scikit-learn compatible estimator or meta-estimator,
you can check our related developer guide:
:ref:`sphx_glr_auto_examples_miscellaneous_plot_metadata_routing.py`.

.. note::
  Note that the methods and requirements introduced in this document are only
  relevant if you want to pass :term:`metadata` (e.g. ``sample_weight``) to a method.
  If you're only passing ``X`` and ``y`` and no other parameter / metadata to
  methods such as :term:`fit`, :term:`transform`, etc, then you don't need to set
  anything.

Usage Examples
**************
Here we present a few examples to show different common use-cases. The examples
in this section require the following imports and data::

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

Here :class:`~model_selection.GroupKFold` requests ``groups`` by default. However, we
need to explicitly request weights for our scorer and the internal cross validation of
:class:`~linear_model.LogisticRegressionCV`. Both of these *consumers* know how to use
metadata called ``sample_weight``::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(
  ...     sample_weight=True
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).set_fit_request(sample_weight=True)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     cv=GroupKFold(),
  ...     scoring=weighted_acc,
  ... )

Note that in this example, ``my_weights`` is passed to both the scorer and
:class:`~linear_model.LogisticRegressionCV`.

Error handling: if ``props={"sample_weigh": my_weights, ...}`` were passed
(note the typo), :func:`~model_selection.cross_validate` would raise an error,
since ``sample_weigh`` was not requested by any of its underlying objects.

Weighted scoring and unweighted fitting
---------------------------------------

When passing metadata such as ``sample_weight`` around, all ``sample_weight``
:term:`consumers <consumer>` require weights to be either explicitly requested
or not requested (i.e. ``True`` or ``False``) when used in another
:term:`router` such as a :class:`~pipeline.Pipeline` or a ``*GridSearchCV``. To
perform an unweighted fit, we need to configure
:class:`~linear_model.LogisticRegressionCV` to not request sample weights, so
that :func:`~model_selection.cross_validate` does not pass the weights along::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(
  ...     sample_weight=True
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).set_fit_request(sample_weight=False)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

If :meth:`linear_model.LogisticRegressionCV.set_fit_request` has not
been called, :func:`~model_selection.cross_validate` will raise an
error because ``sample_weight`` is passed in but
:class:`~linear_model.LogisticRegressionCV` would not be explicitly configured
to recognize the weights.

Unweighted feature selection
----------------------------

Setting request values for metadata are only required if the object, e.g. estimator,
scorer, etc., is a consumer of that metadata Unlike
:class:`~linear_model.LogisticRegressionCV`, :class:`~feature_selection.SelectKBest`
doesn't consume weights and therefore no request value for ``sample_weight`` on its
instance is set and ``sample_weight`` is not routed to it::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(
  ...     sample_weight=True
  ... )
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
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Advanced: Different scoring and fitting weights
-----------------------------------------------

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
  ...     props={
  ...         "scoring_weight": my_weights,
  ...         "fitting_weight": my_other_weights,
  ...         "groups": my_groups,
  ...     },
  ...     scoring=weighted_acc,
  ... )

API Interface
*************

A :term:`consumer` is an object (estimator, meta-estimator, scorer, splitter)
which accepts and uses some :term:`metadata` in at least one of its methods
(``fit``, ``predict``, ``inverse_transform``, ``transform``, ``score``,
``split``). Meta-estimators which only forward the metadata to other objects
(the child estimator, scorers, or splitters) and don't use the metadata
themselves are not consumers. (Meta-)Estimators which route metadata to other
objects are :term:`routers <router>`. A(n) (meta-)estimator can be a
:term:`consumer` and a :term:`router` at the same time. (Meta-)Estimators and
splitters expose a ``set_*_request`` method for each method which accepts at
least one metadata. For instance, if an estimator supports ``sample_weight`` in
``fit`` and ``score``, it exposes
``estimator.set_fit_request(sample_weight=value)`` and
``estimator.set_score_request(sample_weight=value)``. Here ``value`` can be:

- ``True``: method requests a ``sample_weight``. This means if the metadata is
  provided, it will be used, otherwise no error is raised.
- ``False``: method does not request a ``sample_weight``.
- ``None``: router will raise an error if ``sample_weight`` is passed. This is
  in almost all cases the default value when an object is instantiated and
  ensures the user sets the metadata requests explicitly when a metadata is
  passed. The only exception are ``Group*Fold`` splitters.
- ``"param_name"``: if this estimator is used in a meta-estimator, the
  meta-estimator should forward ``"param_name"`` as ``sample_weight`` to this
  estimator. This means the mapping between the metadata required by the
  object, e.g. ``sample_weight`` and what is provided by the user, e.g.
  ``my_weights`` is done at the router level, and not by the object, e.g.
  estimator, itself.

Metadata are requested in the same way for scorers using ``set_score_request``.

If a metadata, e.g. ``sample_weight``, is passed by the user, the metadata
request for all objects which potentially can consume ``sample_weight`` should
be set by the user, otherwise an error is raised by the router object. For
example, the following code raises an error, since it hasn't been explicitly
specified whether ``sample_weight`` should be passed to the estimator's scorer
or not::

    >>> param_grid = {"C": [0.1, 1]}
    >>> lr = LogisticRegression().set_fit_request(sample_weight=True)
    >>> try:
    ...     GridSearchCV(
    ...         estimator=lr, param_grid=param_grid
    ...     ).fit(X, y, sample_weight=my_weights)
    ... except ValueError as e:
    ...     print(e)
    [sample_weight] are passed but are not explicitly set as requested or not for
    LogisticRegression.score

The issue can be fixed by explicitly setting the request value::

    >>> lr = LogisticRegression().set_fit_request(
    ...     sample_weight=True
    ... ).set_score_request(sample_weight=False)
