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
# This meta-estimator allows to take an estimator and freeze its fit methods, meaning
# that calling `fit` does not perform any operations; also, `fit_predict` and
# `fit_transform` call `predict` and `transform` respectively without calling `fit`. The
# original estimator's other methods and properties are left unchanged. An interesting
# use case for this is to use a pre-fitted model as a transformer step in a pipeline,
# or to pass a pre-fitted model to some of the meta-estimators. Here's a short example:

from sklearn.datasets import make_classification
from sklearn.frozen import FrozenEstimator
from sklearn.linear_model import SGDClassifier
from sklearn.model_selection import FixedThresholdClassifier

X, y = make_classification(n_samples=1000, random_state=0)
classifier = SGDClassifier().fit(X, y)

threshold_classifier = FixedThresholdClassifier(
    estimator=FrozenEstimator(classifier), threshold=0.9
)

# %%
# For more details refer to :ref:`plot_frozen_estimator_example`.

# %%
# Missing value support for Extra Trees
# -------------------------------------
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
# The function :func:`datasets.fetch_file` allows to download any file from a given url.

# %%
# Array API support
# -----------------
# Many more estimators and functions have been updated to support Array API compatible
# inputs since version 1.5, in particular the meta-estimators for hyperparameter tuning
# from the :mod:`sklearn.model_selection` module and the metrics from the
# :mod:`sklearn.metrics` module.

# %%
# Almost complete Metadata Routing support
# ----------------------------------------
# Support for routing metadata has been added to nearly all remaining estimators and
# functions. See :ref:`Metadata Routing User Guide <metadata_routing>` for more details.

# %%
# Scikit-learn supports Python 3.13 Free-threaded
# -----------------------------------------------

# %%
# Improvements to the developer API for third party libraries
# -----------------------------------------------------------
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
#   messages. We've also added some new tests, which should help you easier fix
#   potential issues with your estimators.
#
# An updated version of our :ref:`develop` is also available which we recommend you
# to check out.
