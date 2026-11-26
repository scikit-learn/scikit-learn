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
# Decision trees, random forests and their extra-trees variants now handle categorical
# features natively, like :class:`~ensemble.HistGradientBoostingClassifier` already
# does. The new `categorical_features` parameter is available in
# :class:`~tree.DecisionTreeClassifier`, :class:`~tree.ExtraTreeClassifier`,
# :class:`~ensemble.RandomForestClassifier`, :class:`~ensemble.ExtraTreesClassifier`
# and their regressor counterparts.
#
# No encoder is needed: with `categorical_features="from_dtype"`, the categorical
# columns of a pandas or polars dataframe are detected automatically, even if they
# contain strings. A split can then send any group of categories to each child,
# instead of being limited to a threshold on an arbitrary ordering of the categories,
# such as ordinal codes.
#
# In the toy example below, the target is 1 for three of six colors. A single native
# categorical split separates them perfectly. On ordinal codes, one split is not
# enough: here the tree needs three.
import numpy as np
import pandas as pd

from sklearn.tree import DecisionTreeClassifier

rng = np.random.default_rng(0)
colors = ["red", "green", "blue", "yellow", "purple", "orange"]
X = pd.DataFrame({"color": pd.Categorical(rng.choice(colors, size=1_000))})
y = X["color"].isin(["red", "blue", "purple"]).astype(int)
X_codes = X["color"].cat.codes.to_frame()

native = DecisionTreeClassifier(max_depth=1, categorical_features="from_dtype")
print(f"One native categorical split: accuracy {native.fit(X, y).score(X, y):.2f}")

ordinal = DecisionTreeClassifier(max_depth=1).fit(X_codes, y)
print(f"One split on ordinal codes:   accuracy {ordinal.score(X_codes, y):.2f}")

ordinal_deep = DecisionTreeClassifier().fit(X_codes, y)
n_splits = ordinal_deep.get_n_leaves() - 1
print(f"Splits needed on ordinal codes for accuracy 1.00: {n_splits}")

# %%
# Decision trees (with the default `splitter="best"`) and random forests support
# categorical features for binary classification and single-output regression, with
# up to 255 categories per feature. Extra trees use random splits and also support
# multiclass and multi-output targets, with up to ``2**24 - 1`` categories per
# feature.

# %%
# Lasso and Elastic-Net penalties for GLMs and new solvers
# ---------------------------------------------------------
# The generalized linear models :class:`~linear_model.PoissonRegressor`,
# :class:`~linear_model.GammaRegressor` and :class:`~linear_model.TweedieRegressor`
# can now fit :math:`\ell_1` (Lasso) and Elastic-Net penalties with the new `l1_ratio`
# parameter. This brings feature selection to models for counts and for positive,
# skewed targets. On top of that, you no longer need to pick a solver yourself.
#
# **Two new solvers.** `"newton-cd"` and `"newton-cd-gram"` support :math:`\ell_1` and
# Elastic-Net penalties. They are available for the three models above as well as for
# :class:`~linear_model.LogisticRegression` and
# :class:`~linear_model.LogisticRegressionCV`. Both run Newton iterations around the
# coordinate descent solver of :class:`~linear_model.ElasticNet`, and so benefit from
# its large speed-ups over the last two releases.
#
# - Prefer `"newton-cd-gram"` when `n_samples > n_features`. It builds the full Hessian
#   matrix, which is square with `n_features` rows.
# - Prefer `"newton-cd"` when `n_features > n_samples`. It never builds that matrix.
#
# For multiclass problems, replace `n_features` by `n_features * n_classes` in the
# two points above (binary classification effectively has `n_classes = 1`).
#
# **Multiclass included.** Both solvers fit the full multinomial loss of
# :class:`~linear_model.LogisticRegression`. Until now, `"saga"` was the only solver
# for a multinomial model with an :math:`\ell_1` or Elastic-Net penalty. Getting
# coordinate descent there was by far the hardest part of this work, because the
# Hessian is not diagonal across classes. For instance, `"newton-cd"` uses an
# analytical :math:`LDL^\top` decomposition of the per-sample Hessian (see
# :doi:`Tanabe & Sagae (1992) <10.1111/J.2517-6161.1992.TB01875.X>`) to avoid
# constructing the full Hessian matrix.
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
# **Penalty strength alpha instead of C.** The penalty strength of
# :class:`~linear_model.LogisticRegression` is now set by `alpha = 1 / C`, as in
# :class:`~linear_model.Ridge`, :class:`~linear_model.Lasso` and almost all other
# penalized linear models: the larger `alpha`, the stronger the penalty. The default
# `alpha=1.0` gives the same model as the previous default `C=1.0`. `C` is deprecated
# and will be removed in 1.14.
#
# In the example below, the :math:`\ell_1` penalty sets the coefficients of the
# features that do not influence the target to exactly zero, for a Poisson regression
# on counts and for a logistic regression on whether the count is positive. Only
# `l1_ratio` is set; the solver is selected automatically.
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
# Easier metadata routing through auto-requests
# ---------------------------------------------
# Metadata routing no longer needs boilerplate. Once it is enabled, scikit-learn
# estimators and scorers auto-request `sample_weight`, and estimators that use a
# validation set also auto-request `X_val`, `y_val` and `sample_weight_val`. Calls to
# `set_fit_request` or `set_score_request` are only needed to opt out or for other
# metadata. Metadata routing itself is still experimental and off by default.
#
# In the example below, 30% of the labels are flipped and those samples get zero
# weight. A single `sample_weight` argument reaches the scaler, the classifier and the
# scorer, so the corrupted samples are ignored everywhere. Without metadata routing,
# the weighted call fails because a :class:`~pipeline.Pipeline` does not know which
# step should receive the weights.
import numpy as np

from sklearn import config_context
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

X, y = make_classification(n_samples=500, random_state=0)
flipped = np.random.default_rng(0).random(len(y)) < 0.3
y_noisy = np.where(flipped, 1 - y, y)
sample_weight = np.where(flipped, 0.0, 1.0)

pipe = make_pipeline(StandardScaler(), LogisticRegression())
search = GridSearchCV(pipe, param_grid={"logisticregression__alpha": [0.1, 1, 10]})

with config_context(enable_metadata_routing=True):
    search.fit(X, y_noisy)
    print(f"Accuracy without weights: {search.best_score_:.2f}")
    search.fit(X, y_noisy, sample_weight=sample_weight)
    print(f"Accuracy with weights:    {search.best_score_:.2f}")

# %%
# An explicit `set_{method}_request` still takes precedence over the auto-requests:
# `set_fit_request(sample_weight=False)` opts out of receiving `sample_weight`.
# `set_config(enable_metadata_auto_requests=False)` restores the previous behaviour,
# where passing metadata that was not explicitly requested raises an
# `UnsetMetadataPassedError`.
#
# A validation set for the early stopping of
# :class:`~ensemble.HistGradientBoostingClassifier` and
# :class:`~ensemble.HistGradientBoostingRegressor` travels the same way: pass `X_val`
# and `y_val` to `fit`. Inside a :class:`~pipeline.Pipeline`, `X_val` is now
# transformed by the preceding steps by default (`transform_input=("X_val",)`). See
# :ref:`sphx_glr_auto_examples_compose_plot_metadata_routing_validation_set.py`.
#
# For more details, see :ref:`auto-requested metadata
# <metadata_routing_auto_request_user>`. Third-party estimators do not auto-request
# metadata unless they implement it; see :ref:`metadata_routing_auto_request`.

# %%
# Performance update
# ------------------
# TODO
