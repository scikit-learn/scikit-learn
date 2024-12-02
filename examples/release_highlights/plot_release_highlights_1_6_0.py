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
# that calling `fit`, `fit_predict` or `fit_transform` will have no effect. Its other
# methods and properties are left unchanged.
# An interesting use case for this is to use a pre-fitted model as a transformer step in
# a pipeline.

# %%
# Transforming metadata in a Pipeline
# -----------------------------------

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
# A developer API for third party libraries
# -----------------------------------------
