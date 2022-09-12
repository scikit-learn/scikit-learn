
.. _metadata_routing:

.. TODO: update doc/conftest.py once document is updated and examples run.

Metadata Routing
================

This guide demonstrates how metadata such as ``sample_weight`` can be routed
and passed along to estimators, scorers, and CV splitters through
meta-estimators such as ``Pipeline`` and ``GridSearchCV``. In order to pass
metadata to a method such as ``fit`` or ``score``, the object accepting the
metadata, must *request* it. For estimators and splitters this is done via
``set_*_request`` methods, e.g. ``set_fit_request(...)``, and for scorers this
is done via ``set_score_request`` method. For grouped splitters such as
``GroupKFold`` a ``groups`` parameter is requested by default. This is best
demonstrated by the following examples.

If you are developing a scikit-learn compatible estimator or meta-estimator,
you can check our related developer guide:
:ref:`sphx_glr_auto_examples_plot_metadata_routing.py`.

Usage Examples
**************
Here we present a few examples to show different common use-cases. The examples
in this section require the following imports and data::

  >>> import numpy as np
  >>> from sklearn.metrics import make_scorer, accuracy_score
  >>> from sklearn.linear_model import LogisticRegressionCV
  >>> from sklearn.linear_model import LogisticRegression
  >>> from sklearn.model_selection import cross_validate
  >>> from sklearn.model_selection import GridSearchCV
  >>> from sklearn.model_selection import GroupKFold
  >>> from sklearn.feature_selection import SelectKBest
  >>> from sklearn.utils.metadata_requests import RequestType
  >>> from sklearn.pipeline import make_pipeline
  >>> n_samples, n_features = 100, 4
  >>> X = np.random.rand(n_samples, n_features)
  >>> y = np.random.randint(0, 2, size=n_samples)
  >>> my_groups = np.random.randint(0, 10, size=n_samples)
  >>> my_weights = np.random.rand(n_samples)
  >>> my_other_weights = np.random.rand(n_samples)

Weighted scoring and fitting
----------------------------

Here ``GroupKFold`` requests ``groups`` by default. However, we need to
explicitly request weights for our scorer and for ``LogisticRegressionCV``.
Both of these *consumers* know how to use metadata called ``"sample_weight"``::

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
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Note that in this example, ``my_weights`` is passed to both the scorer and
:class:`~linear_model.LogisticRegressionCV`.

Error handling: if ``props={"sample_weigh": my_weights, ...}`` were passed
(note the typo), ``cross_validate`` would raise an error, since
``sample_weigh`` was not requested by any of its children.

Weighted scoring and unweighted fitting
---------------------------------------

All scikit-learn estimators requires weights to be either explicitly requested
or not requested (i.e. ``UNREQUESTED``) when used in another router such as a
``Pipeline`` or a ``*GridSearchCV``. To perform a unweighted fit, we need to
configure :class:`~linear_model.LogisticRegressionCV` to not request sample
weights, so that :func:`~model_selection.cross_validate` does not pass the
weights along::

  >>> weighted_acc = make_scorer(accuracy_score).set_score_request(
  ...     sample_weight=True
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).set_fit_request(sample_weight=RequestType.UNREQUESTED)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Note the usage of ``RequestType`` which in this case is equivalent to
``False``; the type is explained further at the end of this document.

If :class:`~linear_model.LogisticRegressionCV` does not call
``set_fit_request``, :func:`~model_selection.cross_validate` will raise an
error because weights is passed in but
:class:`~linear_model.LogisticRegressionCV` would not be explicitly configured
to recognize the weights.

Unweighted feature selection
----------------------------

Unlike ``LogisticRegressionCV``, ``SelectKBest`` doesn't accept weights and
therefore `"sample_weight"` is not routed to it::

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

Despite ``make_scorer`` and ``LogisticRegressionCV`` both expecting the key
``sample_weight``, we can use aliases to pass different weights to different
consumers. In this example, we pass ``scoring_weight`` to the scorer, and
``fitting_weight`` to ``LogisticRegressionCV``::

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

A *consumer* is an object (estimator, meta-estimator, scorer, splitter) which
accepts and uses some metadata in at least one of its methods (``fit``,
``predict``, ``inverse_transform``, ``transform``, ``score``, ``split``).
Meta-estimators which only forward the metadata to other objects (the child
estimator, scorers, or splitters) and don't use the metadata themselves are not
consumers. (Meta)Estimators which route metadata to other objects are
*routers*. An (meta)estimator can be a consumer and a router at the same time.
(Meta)Estimators and splitters expose a ``set_*_request`` method for each
method which accepts at least one metadata. For instance, if an estimator
supports ``sample_weight`` in ``fit`` and ``score``, it exposes
``estimator.set_fit_request(sample_weight=value)`` and
``estimator.set_score_request(sample_weight=value)``. Here ``value`` can be:

- ``RequestType.REQUESTED`` or ``True``: method requests a ``sample_weight``.
  This means if the metadata is provided, it will be used, otherwise no error
  is raised.
- ``RequestType.UNREQUESTED`` or ``False``: method does not request a
  ``sample_weight``.
- ``RequestType.ERROR_IF_PASSED`` or ``None``: router will raise an error if
  ``sample_weight`` is passed. This is in almost all cases the default value
  when an object is instantiated and ensures the user sets the metadata
  requests explicitly when a metadata is passed. The only exception are
  ``Group*Fold`` splitters.
- ``"param_name"``: if this estimator is used in a meta-estimator, the
  meta-estimator should forward ``"param_name"`` as ``sample_weight`` to this
  estimator. This means the mapping between the metadata required by the
  object, e.g. ``sample_weight`` and what is provided by the user, e.g.
  ``my_weights`` is done at the router level, and not by the object, e.g.
  estimator, itself.

For the scorers, this is done the same way, using ``set_score_request`` method.

If a metadata, e.g. ``sample_weight``, is passed by the user, the metadata
request for all objects which potentially can accept ``sample_weight`` should
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
