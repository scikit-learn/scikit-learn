
.. _metadata_routing:

Metadata Routing
================

This guide demonstrates how metadata such as ``sample_weight`` can be routed
and passed along to estimators, scorers, and CV splitters through
meta-estimators such as ``Pipeline`` and ``GridSearchCV``. In order to pass a
metadata to a method such as ``fit`` or ``score``, the object accepting the
metadata, should *request* it. For estimators and splitters this is done via
``request_*`` methods, e.g. ``request_sample_weight(...)``, and for scorers
this is done via passing ``request_props`` to ``make_scorer``. For grouped
splitters such as ``GroupKFold`` a ``groups`` parameter is requested by
default. This is best demonstrated by the following examples.

Usage Examples
**************
Here we present a few examples to show different common usecases. The examples
in this section require the following imports and data:

  >>> import numpy as np
  >>> from sklearn.metrics import make_scorer, accuracy_score
  >>> from sklearn.linear_model import LogisticRegressionCV
  >>> from sklearn.model_selection import cross_validate
  >>> from sklearn.model_selection import GroupKFold
  >>> from sklearn.feature_selection import SelectKBest
  >>> from sklearn.pipeline import make_pipeline
  >>> N, M = 100, 4
  >>> X = np.random.rand(N, M)
  >>> y = np.random.randint(0, 2, size=N)
  >>> my_groups = np.random.randint(0, 10, size=N)
  >>> my_weights = np.random.rand(N)
  >>> my_other_weights = np.random.rand(N)

Weighted scoring and fitting
----------------------------

Here GroupKFold requests ``groups`` by default. We need to explicitly request
weights in ``make_scorer`` and for ``LogisticRegressionCV``. Both of these
consumers understand the meaning of the key ``"sample_weight"``::

  >>> weighted_acc = make_scorer(accuracy_score,
  ...                            request_props=["sample_weight"])
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ...     ).request_sample_weight(fit="sample_weight")
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Error handling: if props={'sample_eight': my_weights, ...} was passed,
cross_validate would raise an error, since 'sample_eight' was not
requested by any of its children.

Weighted scoring and unweighted fitting
---------------------------------------

Since ``LogisticRegressionCV``, like all scikit-learn estimators, requires that
weights explicitly be requested, removing that request means the fitting is
unweighted::

  >>> weighted_acc = make_scorer(accuracy_score,
  ...                            request_props=["sample_weight"])
  >>> lr = LogisticRegressionCV(cv=GroupKFold(), scoring=weighted_acc,)
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

Like ``LogisticRegressionCV``, ``SelectKBest`` needs to request weights
explicitly. Here it does not request them::

  >>> weighted_acc = make_scorer(accuracy_score,
  ...                            request_props=["sample_weight"])
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).request_sample_weight(fit=True)
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

Despite make_scorer and LogisticRegressionCV both expecting a key
sample_weight, we can use aliases to pass different weights to different
consumers::

  >>> weighted_acc = make_scorer(
  ...     accuracy_score, request_props={"scoring_weight": "sample_weight"}
  ... )
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).request_sample_weight(fit="fitting_weight")
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
Meta-estimators which only forward the metadata to the child estimator and
don't use the metadata themselves are not consumers. (Meta)Estimators and
splitters expose a ``request_*`` method for each metadata they accept. For
instance, if an estimator supports ``sample_weight`` in ``fit`` and ``score``,
it exposes ``estimator.request_sample_weight(fit=value, score=value)``. Here
``value`` can be:

- ``True``: method requests a ``sample_weight``.
- ``False``: method does not request a ``sample_weight``.
- ``"param_name"``: if this estimator is used in a meta-estimator, the
  meta-estimator should forward ``"param_name"`` as ``sample_weight`` to this
  estimator.

For the scorers, on the other hand, the user sets the routing via
``make_scorer`` which accepts a ``request_props`` keyword argument, which is
defined as::

    request_props : list of strings, or dict of {str: str}, default=None
        A list of required properties, or a mapping of the form
        ``{provided_metadata: required_metadata}``, or None.
