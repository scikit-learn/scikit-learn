# ruff: noqa: CPY001
"""
=======================================
Release Highlights for scikit-learn 1.7
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.7! Many bug fixes
and improvements were added, as well as some key new features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_7>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Improved estimator's HTML representation
# ----------------------------------------
# The HTML representation of estimators now includes a section containing the list of
# parameters and their values. Non-default parameters are highlighted in orange. A copy
# button is also available to copy the "fully-qualified" parameter name without the
# need to call the `get_params` method. It is particularly useful when defining a
# parameter grid for a grid-search or a randomized-search with a complex pipeline.
#
# See the example below and click on the different estimator's blocks to see the
# improved HTML representation.

from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

model = make_pipeline(StandardScaler(with_std=False), LogisticRegression(C=2.0))
model

# %%
# Custom validation set for histogram-based Gradient Boosting estimators
# ----------------------------------------------------------------------
# The :class:`ensemble.HistGradientBoostingClassifier` and
# :class:`ensemble.HistGradientBoostingRegressor` now support directly passing a custom
# validation set for early stopping to the `fit` method, using the `X_val`, `y_val`, and
# `sample_weight_val` parameters.
# In a :class:`pipeline.Pipeline`, the validation set `X_val` can be transformed along
# with `X` using the `transform_input` parameter.

import sklearn
from sklearn.datasets import make_classification
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

sklearn.set_config(enable_metadata_routing=True)

X, y = make_classification(random_state=0)
X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=0)

clf = HistGradientBoostingClassifier()
clf.set_fit_request(X_val=True, y_val=True)

model = Pipeline([("sc", StandardScaler()), ("clf", clf)], transform_input=["X_val"])
model.fit(X, y, X_val=X_val, y_val=y_val)

# %%
# Plotting ROC curves from cross-validation results
# -------------------------------------------------
# The class :class:`metrics.RocCurveDisplay` has a new class method `from_cv_results`
# that allows to easily plot multiple ROC curves from the results of
# :func:`model_selection.cross_validate`.

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import RocCurveDisplay
from sklearn.model_selection import cross_validate

X, y = make_classification(n_samples=150, random_state=0)
clf = LogisticRegression(random_state=0)
cv_results = cross_validate(clf, X, y, cv=5, return_estimator=True, return_indices=True)
_ = RocCurveDisplay.from_cv_results(cv_results, X, y)

# %%
# Array API support
# -----------------
# Several functions have been updated to support array API compatible inputs since
# version 1.6, especially metrics from the :mod:`sklearn.metrics` module.
#
# In addition, it is no longer required to install the `array-api-compat` package to use
# the experimental array API support in scikit-learn.
#
# Please refer to the :ref:`array API support<array_api>` page for instructions to use
# scikit-learn with array API compatible libraries such as PyTorch or CuPy.

# %%
# Improved API consistency of Multi-layer Perceptron
# --------------------------------------------------
# The :class:`neural_network.MLPRegressor` has a new parameter `loss` and now supports
# the "poisson" loss in addition to the default "squared_error" loss.
# Moreover, the :class:`neural_network.MLPClassifier` and
# :class:`neural_network.MLPRegressor` estimators now support sample weights.
# These improvements have been made to improve the consistency of these estimators
# with regard to the other estimators in scikit-learn.

# %%
# Migration toward sparse arrays
# ------------------------------
# In order to prepare `SciPy migration from sparse matrices to sparse arrays <https://docs.scipy.org/doc/scipy/reference/sparse.migration_to_sparray.html>`_,
# all scikit-learn estimators that accept sparse matrices as input now also accept
# sparse arrays.
