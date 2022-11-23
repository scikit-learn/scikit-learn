# flake8: noqa
"""
=======================================
Release Highlights for scikit-learn 1.2
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.2! Many bug fixes
and improvements were added, as well as some new key features. We detail
below a few of the major features of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <changes_1_2>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Pandas output with `set_output` API
# -----------------------------------
# scikit-learn's transformers now support pandas output with the `set_output` API.
# To learn more about the `set_output` API see the example:
# :ref:`sphx_glr_auto_examples_miscellaneous_plot_set_output.py` and
# # this `video, pandas DataFrame output for scikit-learn transformers
# (some examples) <https://youtu.be/5bCg8VfX2x8>`__.

import numpy as np
from sklearn.datasets import load_iris
from sklearn.preprocessing import StandardScaler, KBinsDiscretizer
from sklearn.compose import ColumnTransformer

X, y = load_iris(as_frame=True, return_X_y=True)
sepal_cols = ["sepal length (cm)", "sepal width (cm)"]
petal_cols = ["petal length (cm)", "petal width (cm)"]

preprocessor = ColumnTransformer(
    [
        ("scaler", StandardScaler(), sepal_cols),
        ("kbin", KBinsDiscretizer(encode="ordinal"), petal_cols),
    ],
    verbose_feature_names_out=False,
).set_output(transform="pandas")

X_out = preprocessor.fit_transform(X)
X_out.sample(n=5, random_state=0)

# %%
# Interaction constraints in Histogram-based Gradient Boosting Trees
# ------------------------------------------------------------------
# :class:`~ensemble.HistGradientBoostingRegressor` and
# :class:`~ensemble.HistGradientBoostingClassifier` now supports interaction constraints
# with the `interaction_cst` parameter. For details, see the
# :ref:`User Guide <interaction_cst_hgbt>`. In the following example, features are not
# allowed to interact.
from sklearn.datasets import load_diabetes
from sklearn.ensemble import HistGradientBoostingRegressor

X, y = load_diabetes(return_X_y=True, as_frame=True)

hist_no_interact = HistGradientBoostingRegressor(
    interaction_cst=[[i] for i in range(X.shape[1])], random_state=0
)
hist_no_interact.fit(X, y)

# %%
# Experimental Array API support in :class:`~discriminant_analysis.LinearDiscriminantAnalysis`
# --------------------------------------------------------------------------------------------
# Experimental support for the `Array API <https://data-apis.org/array-api/latest/>`_
# specification was added to :class:`~discriminant_analysis.LinearDiscriminantAnalysis`.
# The estimator can now run on any Array API compliant libraries such as
# `CuPy <https://docs.cupy.dev/en/stable/overview.html>`__, a GPU-accelerated array
# library. For details, see the :ref:`User Guide <array_api>`.

# %%
# Improved efficiency of many estimators
# --------------------------------------
# In version 1.1 the efficiency of many estimators relying on the computation of
# pairwise distances was greatly improved for float64 dense input. In version 1.2,
# the efficiency of these estimators was further improved for all combinations of
# float32/float64 and dense/sparse input. It concerns essentially clustering, manifold
# learning and neighbor search algorithms. A detailed list of the impacted estimators
# can be found in the :ref:`changelog <changes_1_2>`. The main benefits are a reduced
# memory footprint and a much better scalability on multi-core machines.
