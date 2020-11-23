# flake8: noqa
"""
========================================
Release Highlights for scikit-learn 0.24
========================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 0.24! Many bug fixes
and improvements were added, as well as some new key features. We detail
below a few of the major features of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <changes_0_24>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

.. meta::
   :keywords: highlights

"""

##############################################################################
# Searching parameter space with successive halving
# -------------------------------------------------
# Successive halving, a state of the art method, is now available to
# explore the space of the parameters and identify their best combination.
# It is still experimental, in order to import it `enable_halving_search_cv`
# should be imported from `experimental`.
# :class:`~sklearn.model_selection.HalvingGridSearchCV` and
# :class:`~sklearn.model_selection.HalvingRandomSearchCV` can be
# used as drop-in replacement for
# :class:`~sklearn.model_selection.GridSearchCV` and
# :class:`~sklearn.model_selection.RandomizedSearchCV`.
# Successive halving is an iterative selection process. The first iteration is
# run with a small amount of resources, by resource meaning typically the
# number of training samples, but also arbitrary integer parameters such
# as `n_estimators` in a random forest. Only some of the parameter candidates
# are selected for the next iteration, with an increasing size of the
# allocated resources.
# As illustrated in the figure below,
# only a subset of candidates will last until the end of the iteration process.
# Read more in the :ref:`User Guide <successive_halving_user_guide>`.
# 
# .. figure:: ../model_selection/images/sphx_glr_plot_successive_halving_iterations_001.png
#   :target: ../model_selection/plot_successive_halving_iterations.html
#   :align: center

import numpy as np
from scipy.stats import randint
from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.model_selection import HalvingRandomSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification

rng = np.random.RandomState(0)

X, y = make_classification(n_samples=700, random_state=rng)

clf = RandomForestClassifier(n_estimators=20, random_state=rng)

param_dist = {"max_depth": [3, None],
              "max_features": randint(1, 11),
              "min_samples_split": randint(2, 11),
              "bootstrap": [True, False],
              "criterion": ["gini", "entropy"]}

rsh = HalvingRandomSearchCV(estimator=clf, param_distributions=param_dist,
                            factor=2, random_state=rng)
rsh.fit(X, y)
rsh.best_params_

##############################################################################
# Native support for categorical features in HistGradientBoosting
# ---------------------------------------------------------------
# :class:`sklearn.ensemble.HistGradientBoostingClassifier` and
# :class:`sklearn.ensemble.HistGradientBoostingRegressor` now have native
# support for categorical features: they can consider splits on non-ordered,
# categorical data. Read more in the :ref:`User Guide
# <categorical_support_gbdt>`.
#
# .. figure:: ../ensemble/images/sphx_glr_plot_gradient_boosting_categorical_001.png
#   :target: ../ensemble/plot_gradient_boosting_categorical.html
#   :align: center

##############################################################################
# Improved performances in HistGradientBoosting methods
# -----------------------------------------------------
# Histogram initialization is now done in parallel in
# :class:`ensemble.HistGradientBoostingRegressor` and
# :class:`ensemble.HistGradientBoostingClassifier` which results in speed
# improvement.
# See more in the `Benchmark page
# <https://scikit-learn.org/scikit-learn-benchmarks/>`_.

##############################################################################
# New SequentialFeatureSelector transformer
# -----------------------------------------
# A new iterative transformer to select features is available:
# :class:`~sklearn.feature_selection.SequentialFeatureSelector`.
# Sequential Feature Selection can add features once at a time (forward
# selection) or remove features from the list of the available features
# (backward selection), based on a cross-validated score maximization.
# See the :ref:`User Guide <sequential_feature_selection>`.

from sklearn.feature_selection import SequentialFeatureSelector
from sklearn.neighbors import KNeighborsClassifier
from sklearn.datasets import load_iris

X, y = load_iris(return_X_y=True, as_frame=True)
feature_names = X.columns
knn = KNeighborsClassifier(n_neighbors=3)
sfs = SequentialFeatureSelector(knn, n_features_to_select=2)
sfs.fit(X, y)
print("Features selected by forward sequential selection: "
      f"{feature_names[sfs.get_support()]}")

##############################################################################
# New PolynomialCountSketch kernel approximation function
# -------------------------------------------------------
# The new :class:`~sklearn.kernel_approximation.PolynomialCountSketch`
# approximates a polynomial expansion of a feature space when used with linear
# models, but uses much less memory than
# :class:`~sklearn.preprocessing.PolynomialFeatures`, and is often faster than
# using a polynomial kernel in SVM.

from sklearn.datasets import fetch_covtype
from sklearn.pipeline import make_pipeline
from sklearn.model_selection import train_test_split
from sklearn.kernel_approximation import PolynomialCountSketch
from sklearn.svm import LinearSVC

X, y = fetch_covtype(return_X_y=True)
pipe = make_pipeline(PolynomialCountSketch(degree=4, n_components=1000),
                     LinearSVC())
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=5000,
                                                    test_size=10000,
                                                    random_state=42)
pipe.fit(X_train, y_train)

##############################################################################
# Retrieving datasets from literature as pandas dataframes
# --------------------------------------------------------
# All the relevant functions loading datasets have a new parameter `as_frame`
# allowing to load the dataset as a pandas DataFrame including columns with
# appropriate dtypes (numeric, string, or categorical) and names.
# The target is a pandas DataFrame or Series depending on the number of
# `target_columns`.

from sklearn.datasets import fetch_kddcup99

df_kddcup99 = fetch_kddcup99(as_frame=True)
df_kddcup99.target.head()

##############################################################################
# Individual Conditional Expectation
# ----------------------------------
# A new kind of partial dependence plot is available: the Individual
# Conditional Expectation (ICE) plot. It visualizes the dependence of the
# prediction on a feature for each sample separately with one line per sample.
# See th :ref:`User Guide <individual_conditional>`

from sklearn.linear_model import BayesianRidge
from sklearn.datasets import fetch_california_housing
from sklearn.inspection import plot_partial_dependence

X, y = fetch_california_housing(return_X_y=True, as_frame=True)
print('Computing partial dependence plots...')
features = ['MedInc', 'AveOccup', 'HouseAge', 'AveRooms']
est = BayesianRidge()
est.fit(X, y)
display = plot_partial_dependence(
       est, X, features, kind="individual", subsample=50,
       n_jobs=3, grid_resolution=20, random_state=0
)
display.figure_.suptitle(
    'Partial dependence of house value on non-location features\n'
    'for the California housing dataset, with BayesianRidge'
)
display.figure_.subplots_adjust(hspace=0.3)

##############################################################################
# New metrics available
# ---------------------
# A number of new metric functions are now available, as for example
# :func:`metrics.top_k_accuracy_score` and :func:`metrics.det_curve`.
# For a complete list see the `changelog
# <../../whats_new/v0.24.html#sklearn-metrics>`_.

##############################################################################
# DecisionTreeRegressor now supports the new 'poisson' splitting criterion 
# ------------------------------------------------------------------------
# The integration of Poisson regression estimation continues from version 0.23
# :class:`~sklearn.tree.DecisionTreeRegressor` now supports a new `'poisson'`
# splitting criterion. Setting `criterion="poisson"` might be a good choice
# if your target is a count or a frequency.

from sklearn.tree import DecisionTreeRegressor

n_samples, n_features = 1000, 20
rng = np.random.RandomState(0)
X = rng.randn(n_samples, n_features)
# positive integer target correlated with X[:, 5] with many zeros:
y = rng.poisson(lam=np.exp(X[:, 5]) / 2)
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=rng)
regressor = DecisionTreeRegressor(criterion='poisson', random_state=0)
regressor.fit(X_train, y_train)

##############################################################################
# New documentation improvements
# ------------------------------
#
# `New examples and documentation pages
# <../../search.html?q=versionadded0.24>`_
# have been addeed in a continuous effort
# to improve the understanding of data science practices.
# Among others a new section about :ref:`common pitfalls and recommended
# practices <common_pitfalls>` is now included, and an `example
# <../model_selection/plot_grid_search_stats.html>`_ illustrating
# how statistically compare the performance of models
# evaluated using :class:`~sklearn.model_selection.GridSearchCV`.
