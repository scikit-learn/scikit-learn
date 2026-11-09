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
# With the default `splitter="best"` (decision trees and random forests), categorical
# features are supported for binary classification and single-output regression,
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
print(f"Non-zero coefficients: {np.sum(glm.coef_ != 0)} out of {n_features}")

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
# Metadata routing: auto-requests
# -------------------------------
# When :ref:`metadata routing <metadata_routing>` is enabled, estimators, scorers and
# meta-estimators now **auto-request** common metadata such as `sample_weight`. You
# no longer need to call `set_fit_request(sample_weight=True)` or
# `set_score_request(sample_weight=True)` on each step: passing `sample_weight` is
# enough. In the example below, `sample_weight` is routed to both steps of the
# pipeline and to the scorer.
import numpy as np

import sklearn
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_validate
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

rng = np.random.default_rng(0)
X, y = make_classification(random_state=0)
sample_weight = rng.uniform(size=y.shape[0])

pipe = make_pipeline(StandardScaler(), LogisticRegression())

with sklearn.config_context(enable_metadata_routing=True):
    results = cross_validate(
        pipe, X, y, scoring="accuracy", params={"sample_weight": sample_weight}
    )
print(f"Weighted accuracy per fold: {results['test_score'].round(2)}")

# %%
# :class:`~ensemble.HistGradientBoostingClassifier` and
# :class:`~ensemble.HistGradientBoostingRegressor` also auto-request a validation
# set (`X_val`, `y_val` and `sample_weight_val`), see
# :ref:`sphx_glr_auto_examples_compose_plot_metadata_routing_validation_set.py`.
#
# Explicit `set_{method}_request` calls still take precedence, for instance to opt
# out of receiving `sample_weight`. Auto-requests can be disabled globally with
# `sklearn.set_config(enable_metadata_auto_requests=False)`. Third-party estimators
# can opt in with `sklearn.utils.metadata_routing.AutoRequestMixin`. See the
# :ref:`User Guide <metadata_routing_auto_request_user>` for details.

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
# TODO: Callbacks
