
.. _metadata_routing:

Metadata Routing
================

This guide demonstrates how metadata such as ``sample_weight`` can be routed
and passed along to estimators, scorers, and CV splitters through
meta-estimators such as ``Pipeline`` and ``GridSearchCV``. In order to pass
metadata to a method such as ``fit`` or ``score``, the object accepting the
metadata, must *request* it. For estimators and splitters this is done via
``*_requests`` methods, e.g. ``fit_requests(...)``, and for scorers this is
done via ``score_requests`` method of a scorer. For grouped splitters such as
``GroupKFold`` a ``groups`` parameter is requested by default. This is best
demonstrated by the following examples.

Usage Examples
**************
Here we present a few examples to show different common use-cases. The examples
in this section require the following imports and data::

.. TODO: add once implemented
  >>> import numpy as np
  >>> from sklearn.metrics import make_scorer, accuracy_score
  >>> from sklearn.linear_model import LogisticRegressionCV
  >>> from sklearn.linear_model import LogisticRegression
  >>> from sklearn.model_selection import cross_validate
  >>> from sklearn.model_selection import GridSearchCV
  >>> from sklearn.model_selection import GroupKFold
  >>> from sklearn.feature_selection import SelectKBest
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
explicitly request weights in ``make_scorer`` and for ``LogisticRegressionCV``.
Both of these *consumers* understand the meaning of the key
``"sample_weight"``::

.. TODO: add once implemented
  >>> weighted_acc = make_scorer(accuracy_score).score_requests(
  ...     sample_weight=True
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).fit_requests(sample_weight=True)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Error handling: if ``props={'sample_weigh': my_weights, ...}`` were passed
(note the typo), cross_validate would raise an error, since 'sample_weigh' was
not requested by any of its children.

Weighted scoring and unweighted fitting
---------------------------------------

Since ``LogisticRegressionCV``, like all scikit-learn estimators, requires that
weights explicitly be requested, we need to explicitly say that
``sample_weight`` is not used for it, so that ``cross_validate`` doesn't pass
it along.

.. TODO: add once implemented
  >>> weighted_acc = make_scorer(accuracy_score).score_requests(
  ...     sample_weight=True
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).fit_requests(sample_weight=False)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Unweighted feature selection
----------------------------

Unlike ``LogisticRegressionCV``, ``SelectKBest`` doesn't accept weights and
therefore `"sample_weight"` is not routed to it::

.. TODO: add once implemented
  >>> weighted_acc = make_scorer(accuracy_score).score_requests(
  ...     sample_weight=True
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).fit_requests(sample_weight=True)
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

Different scoring and fitting weights
-------------------------------------

Despite ``make_scorer`` and ``LogisticRegressionCV`` both expecting a key
``sample_weight``, we can use aliases to pass different weights to different
consumers. In this example, we pass ``scoring_weight`` to the scorer, and
``fitting_weight`` to ``LogisticRegressionCV``::

.. TODO: add once implemented
  >>> weighted_acc = make_scorer(accuracy_score).score_requests(
  ...    sample_weight="scoring_weight"
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).fit_requests(sample_weight="fitting_weight")
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
accepts and uses some metadata in at least one of their methods (``fit``,
``predict``, ``inverse_transform``, ``transform``, ``score``, ``split``).
Meta-estimators which only forward the metadata other objects (the child
estimator, scorers, or splitters) and don't use the metadata themselves are not
consumers. (Meta)Estimators which route metadata to other objects are routers.
An (meta)estimator can be a consumer and a router at the same time.
(Meta)Estimators and splitters expose a ``*_requests`` method for each method
which accepts at least one metadata. For instance, if an estimator supports
``sample_weight`` in ``fit`` and ``score``, it exposes
``estimator.fit_requests(sample_weight=value)`` and
``estimator.score_requests(sample_weight=value)``. Here ``value`` can be:

- ``RequestType.REQUESTED`` or ``True``: method requests a ``sample_weight``.
  This means if the metadata is provided, it will be used, otherwise no error
  is raised.
- ``RequestType.UNREQUESTED`` or ``False``: method does not request a
  ``sample_weight``.
- ``RequestType.ERROR_IF_PASSED`` or ``None``: router will raise an error if
  ``sample_weight`` is passed. This is in almost all cases the default value
  when an object is instantiated and ensures the user sets the metadata
  requests explicitly when a metadata is passed.
- ``"param_name"``: if this estimator is used in a meta-estimator, the
  meta-estimator should forward ``"param_name"`` as ``sample_weight`` to this
  estimator. This means the mapping between the metadata required by the
  object, e.g. ``sample_weight`` and what is provided by the user, e.g.
  ``my_weights`` is done at the router level, and not by the object, e.g.
  estimator, itself.

For the scorers, this is done the same way, using ``.score_requests`` method.

If a metadata, e.g. ``sample_weight`` is passed by the user, the metadata
request for all objects which potentially can accept ``sample_weight`` should
be set by the user, otherwise an error is raised by the router object. For
example, the following code would raise, since it hasn't been explicitly set
whether ``sample_weight`` should be passed to the estimator's scorer or not::

.. TODO: add once implemented
    >>> param_grid = {"C": [0.1, 1]}
    >>> lr = LogisticRegression().fit_requests(sample_weight=True)
    >>> try:
    ...     GridSearchCV(
    ...         estimator=lr, param_grid=param_grid
    ...     ).fit(X, y, sample_weight=my_weights)
    ... except ValueError as e:
    ...     print(e)
    sample_weight is passed but is not explicitly set as requested or not. In
    method: score
