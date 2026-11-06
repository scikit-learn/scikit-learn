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
# TODO: Linear models (LogisticRegression API changes first, then GLM L1/Elastic-Net
# and the "newton-cd" / "newton-cd-gram" solvers)

# %%
# TODO: IterativeImputer is no longer experimental

# %%
# TODO: Metadata routing (auto-requests)

# %%
# TODO: Array API support

# %%
# TODO: Callbacks
