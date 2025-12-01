# ruff: noqa: CPY001
"""
=======================================
Release Highlights for scikit-learn 1.8
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.8! Many bug fixes
and improvements were added, as well as some key new features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_8>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# Array API stuff
# ---------------
# TODO copy and paste from 1.7 highlights needs tweaking
# Several functions have been updated to support array API compatible inputs since
# version 1.7, especially TODO from the :mod:`sklearn.metrics` module.
#
# Please refer to the :ref:`array API support<array_api>` page for instructions to use
# scikit-learn with array API compatible libraries such as PyTorch or CuPy.
#
# TODO are there more "important ones"? could be in the example see point below
# TODO Only show highlighted code without executing it since we don't have a
# GPU in the doc build? We could also show snippet PyTorch CPU with
# commented out device='cuda' if you want to run on GPU you only have to
# uncomment it. Alternative idea link to Colab notebook?

# %%
# Free-threaded CPython 3.14 support
# ----------------------------------
#
# scikit-learn has support for free-threaded CPython, in particular
# free-threaded wheels are available for all of our supported platforms on Python
# 3.14.
#
# Free-threaded (also known as nogil) CPython is a version of CPython that aims at
# enabling efficient multi-threaded use cases by removing the Global Interpreter
# Lock (GIL).
#
# If you want to try out free-threaded Python, the recommendation is to use
# Python 3.14, that has fixed a number of issues compared to Python 3.13. Feel
# free to try free-threaded on your use case and report any issues!
#
# For more details about free-threaded CPython see `py-free-threading doc <https://py-free-threading.github.io>`_,
# in particular `how to install a free-threaded CPython <https://py-free-threading.github.io/installing_cpython/>`_
# and `Ecosystem compatibility tracking <https://py-free-threading.github.io/tracking/>`_.

# %%
# Temperature scaling in `CalibratedClassifierCV`
# -----------------------------------------------
# TODO

# %%
# Linear models improvements
# --------------------------
# There are two main developments going on for linear models: efficiency and API
# changes.
#
# Efficiency of squared error based models with L1 penalty
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# The first one is a massive improvement of efficiency, in particular a reduced fit
# time, for squared error based estimators with L1 penalty: `ElasticNet`, `Lasso`,
# `MultiTaskElasticNet`, `MultiTaskLasso` and their CV variants. The fit time
# improvement is mainly achieved by **gap safe screening rules**. They enable the
# coordinate descent solver to set feature coefficients early to 0 and not look at them
# again. The stronger the L1 penalty the earlier features can be excluded from further
# updates.

from time import time

from sklearn.datasets import make_regression
from sklearn.linear_model import ElasticNetCV

X, y = make_regression(n_features=5000)
model = ElasticNetCV()
tic = time()
model.fit(X, y)
toc = time()
print(f"Fitting ElasticNetCV took {toc - tic:.3} seconds.")

# %%
# API changes in logistic regression
# ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# After the deprecation (version 1.5) and removal (now in 1.8) of the `multi_class`
# parameter in `LogisticRegression` and `LogisticRegressionCV`, the bookkeeping goes
# on. The goal for `LogisticRegressionCV` is to fix types and shapes of a few fitted
# attributes: `C_`, `l1_ratio_`, `coefs_paths_`, `scores_`, `n_iter_`. To ease the
# transition, you can choose between new and old types and shapes with the new and
# temporary parameter `use_legacy_attributes`.

from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegressionCV

X, y = make_classification(n_classes=3, n_informative=4)
model = LogisticRegressionCV(
    use_legacy_attributes=True, l1_ratios=(0,), solver="newton-cholesky"
).fit(X, y)
model.C_  # ndarray of shape (3,), 3 times the same value
model = LogisticRegressionCV(
    use_legacy_attributes=False, l1_ratios=(0,), solver="newton-cholesky"
).fit(X, y)
model.C_  # single float

# %%
# A further deprecation is going on for the `penalty` parameter in both
# `LogisticRegression` and `LogisticRegressionCV`. It is redundant because `C` together
# with `l1_ratio` for `LogisticRegression` and `l1_ratios` for `LogisticRegressionCV`
# contains the same information. Removing `penalty` can ease specifying grid search
# parameter spaces and provides a tidier API.

# %%
# HTML representation of estimators
# ---------------------------------
# TODO

# %%
# DecisionTreeRegressor with MAE
# ------------------------------
# TODO

# %%
# ClassicalMDS
# ------------
# TODO
