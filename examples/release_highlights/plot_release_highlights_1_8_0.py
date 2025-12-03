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
# Array API support
# -----------------
# The progressive adoption of the Python array API standard in SciPy and
# scikit-learn allows the user to pass input arrays from conforming
# libraries to scikit-learn estimators and functions and let them
# use those libraries and possibly non-CPU devices such as GPUs to perform
# the computation instead of attempting to convert all inputs to NumPy.
#
# In scikit-learn 1.8, several estimators and functions have been updated to
# support array API compatible inputs, for example PyTorch tensors and CuPy
# arrays.
#
# Array API support was added to the following estimators:
# :class:`preprocessing.StandardScaler`, :class:`preprocessing.PolynomialFeatures`,
# :class:`linear_model.RidgeCV`, :class:`mixture.GaussianMixture` and
# :class:`calibration.CalibratedClassifierCV`.
#
# Array API support was also added to several metrics in :mod:`sklearn.metrics`
# module, see :ref:`array_api_supported` for more details.
#
# Please refer the :ref:`array API support<array_api>` page for instructions
# to use scikit-learn with array API compatible libraries such as PyTorch or CuPy.
# Note that array API support is still experimental and must be
# explicitly be enabled both in SciPy and scikit-learn to work properly.
#
# TODO do we want to write a snippet?
# - which estimators would we feature?
# - we don't have PyTorch in doc build for now ...
# - we don't have GPU in the doc build but we could show a snippet with numpy
#   and commented out code to switch to PyToch on GPU
# - alternative: show only highlighted code without executing it?
# - alternative: add link to Colab notebook?

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
#
# The long term goal of free-threaded Python is to more efficiently leverage
# multi-core CPUs by using thread workers instead of subprocess workers for parallel
# computation when passing `n_jobs>1` in functions or estimators.
# Efficiency gains are expected by removing the need for inter-process communication.
# Note however that process-based parallelism is still the default joblib backend at
# the time of writing. You can call `joblib.parallel_config(backend="threading")` to
# change the default backend to "threading". Be aware that properly testing that
# everything is running smoothly when doing so is still an ongoing effort and that
# there are open issues to fix before considering making this the default.

# %%
# Temperature scaling in `CalibratedClassifierCV`
# -----------------------------------------------
# Probability calibration of classifiers with temperature scaling is available in
# :class:`calibration.CalibratedClassifierCV` by setting `method = "temperature"`.
# This method is particularly well suited for multiclass problems because it provides
# (better-) calibrated probabilities with just one free parameter, in contrast
# to using a "One-vs-Rest" scheme that adds more parameters for each class.

from sklearn.calibration import CalibratedClassifierCV
from sklearn.datasets import make_classification
from sklearn.frozen import FrozenEstimator
from sklearn.model_selection import train_test_split
from sklearn.svm import LinearSVC

X, y = make_classification(
    n_samples=1000,
    n_features=10,
    n_informative=10,
    n_redundant=0,
    n_classes=5,
    n_clusters_per_class=1,
    class_sep=2.0,
    random_state=42,
)
X_train, X_calib, y_train, y_calib = train_test_split(X, y, random_state=42)
clf = LinearSVC(random_state=42)
clf.fit(X_train, y_train)
ts = CalibratedClassifierCV(FrozenEstimator(clf), method="temperature", ensemble=False)
ts.fit(X_calib, y_calib)
beta = ts.calibrated_classifiers_[0].calibrators[0].beta_
print(f"Optimal temperature = {1 / beta:.3}")

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
