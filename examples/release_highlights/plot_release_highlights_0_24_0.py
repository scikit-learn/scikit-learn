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
# Successive halving, a new, state of the art method, is now available to
# explore the space of the parameters and identify their best combination.
# The parameter space is roughly sampled at the beginning using a small
# amount of resources.
# Only some of the candidates are selected for the next iteration, allowing to
# better sample around the local optimization and to allocate more resources.
# Only a subset of candidates will last until the end of the iteration process.
# Read more in the :ref:`User Guide <successive_halving_user_guide>`.
# 
# .. figure:: ../model_selection/images/sphx_glr_plot_successive_halving_iterations_001.png
#   :target: ../model_selection/plot_successive_halving_iterations.html
#   :align: center

from sklearn.experimental import enable_halving_search_cv  # noqa
from sklearn.model_selection import HalvingGridSearchCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification

param_grid = {'max_depth': [3, 5, 10],
              'min_samples_split': [2, 5, 10]}
base_estimator = RandomForestClassifier(random_state=0)
X, y = make_classification(n_samples=1000, random_state=0)
sh = HalvingGridSearchCV(base_estimator, param_grid, cv=5,
                         factor=2, resource='n_estimators',
                         max_resources=30).fit(X, y)
sh.best_params_

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
