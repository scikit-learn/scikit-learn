# ruff: noqa
"""
=======================================
Release Highlights for scikit-learn 1.6
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.6! Many bug fixes
and improvements were added, as well as some key new features. Below we
detail the highlights of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <release_notes_1_6>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# FrozenEstimator: Freezing an estimator
# --------------------------------------
#
# This meta-estimator allows you to take an estimator and freeze its fit method, meaning
# that calling `fit` does not perform any operations; also, `fit_predict` and
# `fit_transform` call `predict` and `transform` respectively without calling `fit`. The
# original estimator's other methods and properties are left unchanged. An interesting
# use case for this is to use a pre-fitted model as a transformer step in a pipeline
# or to pass a pre-fitted model to some of the meta-estimators. Here's a short example:

import time
from sklearn.datasets import make_classification
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import FixedThresholdClassifier

X, y = make_classification(n_samples=1000, random_state=0)

start = time.time()
classifier = SGDClassifier().fit(X, y)
print(f"Fitting the classifier took {(time.time() - start) * 1_000:.2f} milliseconds")

start = time.time()
threshold_classifier = FixedThresholdClassifier(
    estimator=FrozenEstimator(classifier), threshold=0.9
).fit(X, y)
print(
    f"Fitting the threshold classifier took {(time.time() - start) * 1_000:.2f} "
    "milliseconds"
)

# %%
# Fitting the threshold classifier skipped fitting the inner `SGDClassifier`. For more
# details refer to the example :ref:`sphx_glr_auto_examples_frozen_plot_frozen_examples.py`.

# %%
# Transforming data other than X in a Pipeline
# --------------------------------------------
#
# The :class:`~pipeline.Pipeline` now supports transforming passed data other than `X`
# if necessary. This can be done by setting the new `transform_input` parameter. This
# is particularly useful when passing a validation set through the pipeline.
#
# As an example, imagine `EstimatorWithValidationSet` is an estimator which accepts
# a validation set. We can now have a pipeline which will transform the validation set
# and pass it to the estimator::
#
#     sklearn.set_config(enable_metadata_routing=True)
#     est_gs = GridSearchCV(
#         Pipeline(
#             (
#                 StandardScaler(),
#                 EstimatorWithValidationSet(...).set_fit_request(X_val=True, y_val=True),
#             ),
#             # telling pipeline to transform these inputs up to the step which is
#             # requesting them.
#             transform_input=["X_val"],
#         ),
#         param_grid={"estimatorwithvalidationset__param_to_optimize": list(range(5))},
#         cv=5,
#     ).fit(X, y, X_val=X_val, y_val=y_val)
#
# In the above code, the key parts are the call to `set_fit_request` to specify that
# `X_val` and `y_val` are required by the `EstimatorWithValidationSet.fit` method, and
# the `transform_input` parameter to tell the pipeline to transform `X_val` before
# passing it to `EstimatorWithValidationSet.fit`.
#
# Note that at this time scikit-learn estimators have not yet been extended to accept
# user specified validation sets. This feature is released early to collect feedback
# from third-party libraries who might benefit from it.

# %%
# Multiclass support for `LogisticRegression(solver="newton-cholesky")`
# ---------------------------------------------------------------------
#
# The `"newton-cholesky"` solver (originally introduced in scikit-learn version
# 1.2) was previously limited to binary
# :class:`~linear_model.LogisticRegression` and some other generalized linear
# regression estimators (namely :class:`~linear_model.PoissonRegressor`,
# :class:`~linear_model.GammaRegressor` and
# :class:`~linear_model.TweedieRegressor`).
#
# This new release includes support for multiclass (multinomial)
# :class:`~linear_model.LogisticRegression`.
#
# This solver is particularly useful when the number of features is small to
# medium. It has been empirically shown to converge more reliably and faster
# than other solvers on some medium sized datasets with one-hot encoded
# categorical features as can be seen in the `benchmark results of the
# pull-request
# <https://github.com/scikit-learn/scikit-learn/pull/28840#issuecomment-2065368727>`_.

# %%
# Missing value support for Extra Trees
# -------------------------------------
#
# The classes :class:`ensemble.ExtraTreesClassifier` and
# :class:`ensemble.ExtraTreesRegressor` now support missing values. More details in the
# :ref:`User Guide <tree_missing_value_support>`.
import numpy as np
from sklearn.ensemble import ExtraTreesClassifier

X = np.array([0, 1, 6, np.nan]).reshape(-1, 1)
y = [0, 0, 1, 1]

forest = ExtraTreesClassifier(random_state=0).fit(X, y)
forest.predict(X)

# %%
# Download any dataset from the web
# ---------------------------------
#
# The function :func:`datasets.fetch_file` allows downloading a file from any given URL.
# This convenience function provides built-in local disk caching, sha256 digest
# integrity check and an automated retry mechanism on network error.
#
# The goal is to provide the same convenience and reliability as dataset fetchers while
# giving the flexibility to work with data from arbitrary online sources and file
# formats.
#
# The dowloaded file can then be loaded with generic or domain specific functions such
# as `pandas.read_csv`, `pandas.read_parquet`, etc.

# %%
# Array API support
# -----------------
#
# Many more estimators and functions have been updated to support array API compatible
# inputs since version 1.5, in particular the meta-estimators for hyperparameter tuning
# from the :mod:`sklearn.model_selection` module and the metrics from the
# :mod:`sklearn.metrics` module.
#
# Please refer to the :ref:`array API support<array_api>` page for instructions to use
# scikit-learn with array API compatible libraries such as PyTorch or CuPy.

# %%
# Almost complete Metadata Routing support
# ----------------------------------------
#
# Support for routing metadata has been added to all remaining estimators and
# functions except AdaBoost. See :ref:`Metadata Routing User Guide <metadata_routing>`
# for more details.

# %%
# Free-threaded CPython 3.13 support
# ----------------------------------
#
# scikit-learn has preliminary support for free-threaded CPython, in particular
# free-threaded wheels are available for all of our supported platforms.
#
# Free-threaded (also known as nogil) CPython 3.13 is an experimental version of
# CPython 3.13 which aims at enabling efficient multi-threaded use cases by
# removing the Global Interpreter Lock (GIL).
#
# For more details about free-threaded CPython see `py-free-threading doc <https://py-free-threading.github.io>`_,
# in particular `how to install a free-threaded CPython <https://py-free-threading.github.io/installing_cpython/>`_
# and `Ecosystem compatibility tracking <https://py-free-threading.github.io/tracking/>`_.
#
# Feel free to try free-threaded CPython on your use case and report any issues!

# %%
# Improvements to the developer API for third party libraries
# -----------------------------------------------------------
#
# We have been working on improving the developer API for third party libraries.
# This is still a work in progress, but a fair amount of work has been done in this
# release. This release includes:
#
# - :func:`sklearn.utils.validation.validate_data` is introduced and replaces the
#   previously private `BaseEstimator._validate_data` method. This function extends
#   :func:`~sklearn.utils.validation.check_array` and adds support for remembering
#   input feature counts and names.
# - Estimator tags are now revamped and a part of the public API via
#   :class:`sklearn.utils.Tags`. Estimators should now override the
#   :meth:`BaseEstimator.__sklearn_tags__` method instead of implementing a `_more_tags`
#   method. If you'd like to support multiple scikit-learn versions, you can implement
#   both methods in your class.
# - As a consequence of developing a public tag API, we've removed the `_xfail_checks`
#   tag and tests which are expected to fail are directly passed to
#   :func:`~sklearn.utils.estimator_checks.check_estimator` and
#   :func:`~sklearn.utils.estimator_checks.parametrize_with_checks`. See their
#   corresponding API docs for more details.
# - Many tests in the common test suite are updated and raise more helpful error
#   messages. We've also added some new tests, which should help you more easily fix
#   potential issues with your estimators.
#
# An updated version of our :ref:`develop` is also available, which we recommend you
# check out.
