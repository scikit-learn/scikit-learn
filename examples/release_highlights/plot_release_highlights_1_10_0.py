# ruff: noqa: CPY001
"""
========================================
Release Highlights for scikit-learn 1.10
========================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.10! Many bug fixes
and improvements were added, as well as some key new features. Below we
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
# The tree models in :mod:`sklearn.tree`, as well as the random forests and extra
# trees ensembles (:class:`~ensemble.RandomForestClassifier`,
# :class:`~ensemble.RandomForestRegressor`, :class:`~ensemble.ExtraTreesClassifier`
# and :class:`~ensemble.ExtraTreesRegressor`), now accept a `categorical_features`
# parameter, like :class:`~ensemble.HistGradientBoostingClassifier` already does.
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
# Linear models
# -------------
#
# Penalty parameter `alpha` in logistic regression
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# The penalty strength of :class:`~linear_model.LogisticRegression` and
# :class:`~linear_model.LogisticRegressionCV` is now set with `alpha` (`alphas` for
# `LogisticRegressionCV`), as in :class:`~linear_model.Ridge` or
# :class:`~linear_model.ElasticNet`: the larger the value, the stronger the penalty.
# For `LogisticRegression`, `alpha = 1 / C`. The `C` and `Cs` parameters are
# deprecated and will be removed in 1.14. The new defaults give the same models as
# before. Below, `C=10` and `alpha=0.1` give the same model.
import warnings

import numpy as np

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression

X, y = make_classification(random_state=0)

# Before: `C` is the inverse of the penalty strength (deprecated).
with warnings.catch_warnings():
    warnings.simplefilter("ignore", FutureWarning)
    clf_before = LogisticRegression(C=10).fit(X, y)

# After: `alpha` is the penalty strength, with `alpha = 1 / C`.
clf_after = LogisticRegression(alpha=0.1).fit(X, y)

print(f"Same coefficients: {np.allclose(clf_before.coef_, clf_after.coef_)}")

# %%
# New solvers
# ^^^^^^^^^^^
# Two new solvers, `"newton-cd"` and `"newton-cd-gram"`, support :math:`\ell^1` and
# Elastic-Net penalties. They are available for
# :class:`~linear_model.LogisticRegression` and for the generalized linear models
# :class:`~linear_model.PoissonRegressor`, :class:`~linear_model.GammaRegressor` and
# :class:`~linear_model.TweedieRegressor`, which can now fit :math:`\ell^1` and
# Elastic-Net penalties with the new `l1_ratio` parameter.
#
# - The new default `solver="auto"` selects a solver that works with the chosen
#   `l1_ratio`:
#
#   - `l1_ratio=0` (pure :math:`\ell^2` penalty or no penalty): `"lbfgs"`, as before;
#   - `l1_ratio > 0` and `n_samples >= n_features`: `"newton-cd-gram"`, which
#     builds the full Hessian matrix of size `(n_features, n_features)`;
#   - `l1_ratio > 0` and `n_samples < n_features`: `"newton-cd"`, which avoids
#     building that matrix.
#
#   The selected solver is stored in the new `solver_` attribute. The selection rule
#   may change in future versions without deprecation.
# - In logistic regression, `solver="liblinear"`, as well as the `dual` and
#   `intercept_scaling` parameters, are deprecated and will be removed in 1.12. The
#   new solvers are the recommended replacements.
#
# In the example below, the :math:`\ell^1` penalty sets the coefficients of the
# features that do not influence the target to exactly zero.
import numpy as np

from sklearn.linear_model import PoissonRegressor

rng = np.random.default_rng(0)
n_samples, n_features = 500, 20
X = rng.normal(size=(n_samples, n_features))
true_coef = np.zeros(n_features)
true_coef[:3] = [0.5, -0.4, 0.3]
y = rng.poisson(np.exp(X @ true_coef))

glm = PoissonRegressor(alpha=0.1, l1_ratio=1.0, solver="auto").fit(X, y)
print(f"Selected solver: {glm.solver_}")
print(f"Non-zero coefficients: {np.flatnonzero(glm.coef_)}")

# %%
# :class:`~impute.IterativeImputer` is no longer experimental
# -----------------------------------------------------------
# :class:`~impute.IterativeImputer` is now stable and can be imported directly from
# :mod:`sklearn.impute`. Importing `sklearn.experimental.enable_iterative_imputer`
# is no longer needed: it is now a no-op that raises a warning.
#
# `IterativeImputer` also no longer raises a
# :class:`~exceptions.ConvergenceWarning` when `max_iter` is reached before meeting
# the `tol` stopping criterion: the round-robin imputation is not guaranteed to
# converge, and this is expected. See the :ref:`User Guide
# <iterative_imputer_convergence>` for details.
#
# Similarly, `sklearn.experimental.enable_hist_gradient_boosting` is deprecated and
# will be removed in 1.12, since :class:`~ensemble.HistGradientBoostingClassifier` and
# :class:`~ensemble.HistGradientBoostingRegressor` have been stable for a while.
import numpy as np

from sklearn.impute import IterativeImputer

X = np.array([[1, 2], [3, 6], [4, 8], [np.nan, 3], [7, np.nan]])
print(IterativeImputer(random_state=0).fit_transform(X))

# %%
# Easier Metadata Routing
# -----------------------
# Metadata routing is now much easier to use. It is still experimental and still off by
# default, but once you enable it, scikit-learn estimators auto-request `sample_weight`
# and, where they use a validation set, `X_val`, `y_val`, and `sample_weight_val`. You
# no longer need to call `set_fit_request` or `set_score_request` for that metadata.
# Other metadata (except for `groups`) still needs an explicit request.
#
# In the example below, `sample_weight` is auto-requested on estimator methods that
# accept it (typically `fit` and `score`) and on scorers whose metric accepts it. The
# same weights are used to fit :class:`~linear_model.LogisticRegression` and to score
# each candidate, with no explicit request.
from sklearn import set_config
from sklearn.model_selection import GridSearchCV

set_config(enable_metadata_routing=True)

X, y = make_classification(n_samples=150, random_state=0)
sample_weight = np.ones(len(y))
sample_weight[:15] = 2

search = GridSearchCV(LogisticRegression(), param_grid={"alpha": [0.01, 0.1, 1]})
search.fit(X, y, sample_weight=sample_weight)
print(f"Best alpha: {search.best_params_['alpha']}")
set_config(enable_metadata_routing=False)

# %%
# A request set with `set_{method}_request` still takes precedence.
# `set_fit_request(sample_weight=False)` opts out of receiving `sample_weight`.
# `set_config(enable_metadata_auto_requests=False)` restores the previous behaviour:
# passing metadata that was not explicitly requested raises `UnsetMetadataPassedError`.
#
# :class:`~ensemble.HistGradientBoostingClassifier` and
# :class:`~ensemble.HistGradientBoostingRegressor` now also auto-request a validation
# set (`X_val`, `y_val`, and `sample_weight_val`) for early stopping. See
# :ref:`sphx_glr_auto_examples_compose_plot_metadata_routing_validation_set.py`.
#
# For more information about this new feature, see
# :ref:`auto-requested metadata <metadata_routing_auto_request_user>`.
#
# Estimators outside scikit-learn do not follow this policy unless they implement it
# themselves. For information on the developer API we provide for that purpose, see
# :ref:`metadata_routing_auto_request`.

# %%
# Array API support (enables GPU computations)
# --------------------------------------------
# scikit-learn continues to adopt the Python array API standard, which lets
# estimators and functions work directly on PyTorch tensors or CuPy arrays, including
# on GPUs. In this release, array API support was added to:
#
# - :class:`~linear_model.LogisticRegressionCV` (with `solver="lbfgs"` and
#   `solver="newton-cg"`), and :class:`~linear_model.LogisticRegression` with
#   `solver="newton-cg"`;
# - :class:`~covariance.LedoitWolf` and :class:`~covariance.OAS`;
# - :func:`~metrics.matthews_corrcoef` and
#   :func:`~metrics.cluster.contingency_matrix`.
#
# Several estimators, metrics and cross-validation tools now also accept `y` and
# `sample_weight` from a different namespace or device than `X`. For instance, `X`
# can be on a GPU while `y` stays a NumPy array of string labels on the CPU, which
# matters because PyTorch tensors do not support string values:
#
# .. code-block:: python
#
#     import os
#
#     # Must be set before SciPy is imported (scikit-learn imports SciPy).
#     os.environ["SCIPY_ARRAY_API"] = "1"
#
#     import numpy as np
#     import torch
#
#     import sklearn
#     from sklearn.datasets import make_classification
#     from sklearn.linear_model import LogisticRegressionCV
#
#     X, y = make_classification(random_state=0)
#     X_gpu = torch.asarray(X, dtype=torch.float32, device="cuda")
#     y_labels = np.array(["cat", "dog"])[y]  # NumPy string labels on the CPU
#
#     with sklearn.config_context(array_api_dispatch=True):
#         clf = LogisticRegressionCV(solver="lbfgs", scoring="accuracy")
#         clf.fit(X_gpu, y_labels)
#         print(clf.coef_.device, clf.classes_)
#         print(clf.predict(X_gpu[:5]))
#
# Output::
#
#     cuda:0 ['cat' 'dog']
#     ['cat' 'cat' 'dog' 'dog' 'cat']
#
# Array API support is experimental and must be explicitly enabled both in SciPy and
# scikit-learn. See :ref:`array_api` for instructions and :ref:`array_api_supported`
# for the full list of supported estimators and functions.

# %%
# Callbacks
# ---------
# TODO
