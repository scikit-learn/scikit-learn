# ruff: noqa: CPY001
"""
========================================
Release Highlights for scikit-learn 1.10
========================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.10! Many bug fixes
and improvements were added, as well as some new key features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_10>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Native categorical support in decision trees and random forests
# ----------------------------------------------------------------
# The following tree-based models now accept a `categorical_features` parameter,
# like :class:`~ensemble.HistGradientBoostingClassifier` already does:
#
# - :class:`~tree.DecisionTreeClassifier`
# - :class:`~tree.DecisionTreeRegressor`
# - :class:`~tree.ExtraTreeClassifier`
# - :class:`~tree.ExtraTreeRegressor`
# - :class:`~ensemble.RandomForestClassifier`
# - :class:`~ensemble.RandomForestRegressor`
# - :class:`~ensemble.ExtraTreesClassifier`
# - :class:`~ensemble.ExtraTreesRegressor`
#
# A split on a categorical feature can send any group of categories to each
# child, instead of being limited to a threshold on an arbitrary encoding of the
# categories (such as ordinal codes).
#
# With `categorical_features="from_dtype"`, columns with a categorical dtype in a
# pandas or polars dataframe are detected automatically. In the toy example below,
# the target is 1 for three of six colors. A single native categorical split
# separates them perfectly, while a split on the ordinal codes cannot.
import numpy as np
import pandas as pd

from sklearn.tree import DecisionTreeClassifier

rng = np.random.default_rng(0)
colors = ["red", "green", "blue", "yellow", "purple", "orange"]
X = pd.DataFrame({"color": pd.Categorical(rng.choice(colors, size=1_000))})
y = X["color"].isin(["red", "blue", "purple"]).astype(int)

native = DecisionTreeClassifier(max_depth=1, categorical_features="from_dtype")
ordinal = DecisionTreeClassifier(max_depth=1)

X_codes = X["color"].cat.codes.to_frame()

native.fit(X, y)
ordinal.fit(X_codes, y)
print(f"Native categorical split accuracy: {native.score(X, y):.2f}")
print(f"Ordinal-encoded split accuracy:    {ordinal.score(X_codes, y):.2f}")

# %%
# Decision trees (with the default `splitter="best"`) and random forests support
# categorical features for binary classification and single-output regression,
# with up to 255 categories per feature. Extra trees use random splits and also
# support multiclass and multi-output targets, with up to ``2**24 - 1`` categories
# per feature.

# %%
# Lasso and Elastic-Net penalties for GLMs and new solvers
# ---------------------------------------------------------
# The generalized linear models :class:`~linear_model.PoissonRegressor`,
# :class:`~linear_model.GammaRegressor` and :class:`~linear_model.TweedieRegressor`
# can now fit :math:`\ell_1` (Lasso) and Elastic-Net penalties with the new `l1_ratio`
# parameter. This brings feature selection to models for counts and for positive,
# skewed targets. On top, you no longer need to pick a solver yourself.
#
# **Two new solvers.** `"newton-cd"` and `"newton-cd-gram"` support :math:`\ell_1` and
# Elastic-Net penalties. They are available for the three models above as well as for
# :class:`~linear_model.LogisticRegression` and
# :class:`~linear_model.LogisticRegressionCV`. Both run Newton iterations around the
# coordinate descent solver of :class:`~linear_model.ElasticNet`, and so benefit from
# its large speed-ups over the last two releases.
#
# - Prefer `"newton-cd-gram"` when `n_samples > n_features`. It builds the full
#   Hessian matrix of size `(n_features, n_features)`.
# - Prefer `"newton-cd"` when `n_features > n_samples`. It avoids building that
#   matrix.
#
# **No need to choose a solver.** The new default `solver="auto"` always selects a
# solver that works with the chosen `l1_ratio` (or `l1_ratios` in
# :class:`~linear_model.LogisticRegressionCV`). Switching from an :math:`\ell_2` to an
# :math:`\ell_1` or Elastic-Net penalty is now a matter of changing `l1_ratio` and
# nothing else. With the default `l1_ratio=0` (pure :math:`\ell_2` penalty), `"auto"`
# still selects `"lbfgs"`, so existing models do not change. The selected solver is
# stored in the new `solver_` attribute. The selection rule may change in future
# versions without deprecation.
#
# **`alpha` instead of `C`.** The penalty strength of
# :class:`~linear_model.LogisticRegression` is now set by `alpha = 1 / C`, as in
# :class:`~linear_model.Ridge`, :class:`~linear_model.Lasso` and almost all other
# penalized linear models: the larger `alpha`, the stronger the penalty. The default
# `alpha=1.0` gives the same model as the previous default `C=1.0`. `C` is deprecated
# and will be removed in 1.14.
#
# In the example below, the :math:`\ell_1` penalty sets the coefficients of the
# features that do not influence the target to exactly zero, for a Poisson regression
# on counts and for a logistic regression on whether the count is positive. Only
# `l1_ratio` is set, the solver is selected automatically.
import numpy as np

from sklearn.linear_model import LogisticRegression, PoissonRegressor

rng = np.random.default_rng(0)
n_samples, n_features = 500, 20
X = rng.normal(size=(n_samples, n_features))
true_coef = np.zeros(n_features)
true_coef[:3] = [0.5, -0.4, 0.3]
y = rng.poisson(np.exp(X @ true_coef))

glm = PoissonRegressor(alpha=0.1, l1_ratio=1.0).fit(X, y)
print(f"PoissonRegressor selected solver: {glm.solver_}")
print(f"Non-zero coefficients: {np.flatnonzero(glm.coef_)}")

clf = LogisticRegression(alpha=30, l1_ratio=1.0).fit(X, y > 0)
print(f"LogisticRegression selected solver: {clf.solver_}")
print(f"Non-zero coefficients: {np.flatnonzero(clf.coef_)}")

# %%
# Easier metadata routing
# ------------------------
# Metadata routing is now much easier to use. It is still experimental and still off by
# default, but once you enable it, scikit-learn estimators auto-request `sample_weight`
# and, where they use a validation set, they also auto-request `X_val`, `y_val`, and
# `sample_weight_val`. You no longer need to call `set_fit_request` or
# `set_score_request` for that metadata.
# Other metadata (except for `groups`) still needs an explicit request.
#
# In the example below, `sample_weight` is auto-requested on estimator methods that
# accept it (typically `fit` and `score`) and on scorers whose metric accepts it. The
# same weights are used to fit :class:`~linear_model.LogisticRegression` and to score
# each candidate, with no explicit request.
import numpy as np

from sklearn import set_config
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV

set_config(enable_metadata_routing=True)

X, y = make_classification(n_samples=150, random_state=0)
sample_weight = np.ones(len(y))
sample_weight[:15] = 2

search = GridSearchCV(LogisticRegression(), param_grid={"alpha": [0.01, 0.1, 1]})
search.fit(X, y, sample_weight=sample_weight)
print(f"Best parameters: {search.best_params_}")
set_config(enable_metadata_routing=False)

# %%
# A request set with `set_{method}_request` still takes precedence
# over the auto-requests.
#
# `set_fit_request(sample_weight=False)` opts out of receiving `sample_weight`.
# `set_config(enable_metadata_auto_requests=False)` restores the previous behaviour:
# passing metadata that was not explicitly requested raises `UnsetMetadataPassedError`.
#
# :class:`~ensemble.HistGradientBoostingClassifier` and
# :class:`~ensemble.HistGradientBoostingRegressor` now also auto-request a validation
# set (`X_val`, `y_val`, and `sample_weight_val`) for early stopping. See
# :ref:`sphx_glr_auto_examples_compose_plot_metadata_routing_validation_set.py`.
# In a :class:`~pipeline.Pipeline`, the default of `transform_input` changed to
# `("X_val",)`, so that a validation set passed to `fit` is always transformed by the
# preceding steps alongside `X`.
#
# For more information about this new feature, see
# :ref:`auto-requested metadata <metadata_routing_auto_request_user>`.
#
# Estimators outside scikit-learn do not follow this policy unless they implement it
# themselves. For information on the developer API we provide for that purpose, see
# :ref:`metadata_routing_auto_request`.

# %%
# Performance update
# ------------------
# TODO
