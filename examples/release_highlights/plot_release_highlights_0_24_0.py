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
# :class:`~sklearn.model_selection.HalvingGridSearchCV` and
# :class:`~sklearn.model_selection.HalvingRandomSearchCV` can be
# used as drop-in replacement for
# :class:`~sklearn.model_selection.GridSearchCV` and
# :class:`~sklearn.model_selection.RandomizedSearchCV`.
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
# A new iterative transformer to select features is available:
# :class:`~sklearn.feature_selection.SequentialFeatureSelector`.
# Sequential Feature Selection can add features once at a time (forward
# selection) or remove features from feature list (backward selection),
# based on a cross-validated score maximization.
# See the :ref:`User Guide <sequential_feature_selection>`.

from sklearn.datasets import fetch_openml
from sklearn.compose import make_column_transformer
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import Ridge
from sklearn.feature_selection import SequentialFeatureSelector

survey = fetch_openml(data_id=534, as_frame=True)
X_data = survey.data[survey.feature_names]
y = survey.target.values.ravel()

categorical_columns = ['RACE', 'OCCUPATION', 'SECTOR',
                       'MARR', 'UNION', 'SEX', 'SOUTH']
numerical_columns = ['EDUCATION', 'EXPERIENCE', 'AGE']

preprocessor = make_column_transformer(
    (OneHotEncoder(drop='if_binary'), categorical_columns),
    remainder='passthrough'
)
preprocessor.fit(X_data)
X = preprocessor.transform(X_data)

feature_names = preprocessor.get_feature_names()
ridge = Ridge(alpha=1e-10).fit(X, y)

sfs_forward = SequentialFeatureSelector(ridge, n_features_to_select=2,
                                        direction='forward').fit(X, y)
sfs_backward = SequentialFeatureSelector(ridge, n_features_to_select=2,
                                         direction='backward').fit(X, y)

sfs_forward_features  = [i for (i, v) in zip(feature_names, sfs_forward.get_support()) if v]
sfs_backward_features  = [i for (i, v) in zip(feature_names, sfs_backward.get_support()) if v]

print("Features selected by forward sequential selection: "
      f"{sfs_forward_features}")
print("Features selected by backward sequential selection: "
      f"{sfs_backward_features}")

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
