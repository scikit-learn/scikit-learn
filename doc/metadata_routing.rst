
.. _metadata_routing:

Metadata Routing
================

This guide demonstrates how metadata such as ``sample_weight`` can be routed
and passed along to estimators, scorers, and CV splitters through
meta-estimators such as ``Pipeline`` and ``GridSearchCV``. In order to pass
metadata to a method such as ``fit`` or ``score``, the object accepting the
metadata, must *request* it. For estimators and splitters this is done via
``*_requests`` methods, e.g. ``fit_requests(...)``, and for scorers
this is done via passing ``score_params`` to ``make_scorer``. For grouped
splitters such as ``GroupKFold`` a ``groups`` parameter is requested by
default. This is best demonstrated by the following examples.

Usage Examples
**************
Here we present a few examples to show different common usecases. The examples
in this section require the following imports and data::

  >>> import numpy as np
  >>> from sklearn.metrics import make_scorer, accuracy_score
  >>> from sklearn.linear_model import LogisticRegressionCV
  >>> from sklearn.model_selection import cross_validate
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

Here GroupKFold requests ``groups`` by default. We need to explicitly request
weights in ``make_scorer`` and for ``LogisticRegressionCV``. Both of these
consumers understand the meaning of the key ``"sample_weight"``::

  >>> weighted_acc = make_scorer(accuracy_score,
  ...                            score_params=["sample_weight"])
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ...     ).fit_requests(sample_weight=True)
  >>> cv_results = cross_validate(
  ...     lr,
  ...     X,
  ...     y,
  ...     cv=GroupKFold(),
  ...     props={"sample_weight": my_weights, "groups": my_groups},
  ...     scoring=weighted_acc,
  ... )

Error handling: if props={'sample_weight': my_weights, ...} were passed,
cross_validate would raise an error, since 'sample_weight' was not
requested by any of its children.

Weighted scoring and unweighted fitting
---------------------------------------

Since ``LogisticRegressionCV``, like all scikit-learn estimators, requires that
weights explicitly be requested, we need to explicitly say that ``sample_weight``
is not used for it, so that ``cross_validate`` doesn't pass it along.

  >>> weighted_acc = make_scorer(accuracy_score,
  ...                            score_params=["sample_weight"])
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

Like ``LogisticRegressionCV``, ``SelectKBest`` needs to request weights
explicitly. Here it does not request them::

  >>> weighted_acc = make_scorer(accuracy_score,
  ...                            score_params=["sample_weight"])
  >>> lr = LogisticRegressionCV(
  ...     cv=GroupKFold(), scoring=weighted_acc,
  ... ).fit_requests(sample_weight=True)
  >>> sel = SelectKBest(k=2).fit_requests(sample_weight=False)
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
consumers::

  >>> weighted_acc = make_scorer(
  ...     accuracy_score, score_params={"scoring_weight": "sample_weight"}
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
``make_scorer`` which accepts a ``score_params`` keyword argument, which is
defined as::

    score_params : list of strings, or dict of {str: str}, default=None
        A list of required properties, or a mapping of the form
        ``{"required_metadata": "provided_metadata"}``, or None.
