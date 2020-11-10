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
"""

##############################################################################
# Searching parameter space with successive halving
# -------------------------------------------------
# Successive halving, a state of the art method, is now available to
# explore the space of the parameters and identify their best combination.
# :class:`HalvingGridSearchCV` and :class:`HalvingRandomSearchCV` can be
# used as drop-in replacement for :class:`GridSearchCV` and
# :class:`RandomizedSearchCV`.
# Successive halving is an iterative selection process. The first iteration is
# run on a subset of the input sample. Only some of the parameter candidates
# are selected for the next iteration, allowing to increase the size of the
# input subset. As illustrated in the figure below,
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

rsh = HalvingRandomSearchCV(
    estimator=clf,
    param_distributions=param_dist,
    factor=2,
    random_state=rng)
rsh.fit(X, y)
rsh.best_params_

##############################################################################
# New SequentialFeatureSelector transformer
# -----------------------------------------

##############################################################################
# New PolynomialCountSketch kernel approximation function
# -------------------------------------------------------

##############################################################################
# Individual Conditional Expectation
# ----------------------------------

##############################################################################
# New metrics available
# ---------------------

##############################################################################
# DecisionTreeRegressor now supports the new 'poisson' splitting criterion 
# ------------------------------------------------------------------------

##############################################################################
# Retrieving datasets from literature as pandas dataframes
# --------------------------------------------------------

##############################################################################
# HistGradientBoostingClassifier improved performances
# ----------------------------------------------------
