# ruff: noqa
"""
=======================================
Release Highlights for scikit-learn 1.4
=======================================

.. currentmodule:: sklearn

We are pleased to announce the release of scikit-learn 1.4! Many bug fixes
and improvements were added, as well as some new key features. We detail
below a few of the major features of this release. **For an exhaustive list of
all the changes**, please refer to the :ref:`release notes <changes_1_4>`.

To install the latest version (with pip)::

    pip install --upgrade scikit-learn

or with conda::

    conda install -c conda-forge scikit-learn

"""

# %%
# HistGradientBoosting Natively Supports Categorical DTypes in DataFrames
# -----------------------------------------------------------------------
# :class:`ensemble.HistGradientBoostingClassifier` and
# :class:`ensemble.HistGradientBoostingRegressor` now directly supports dataframes with
# categorical features.  Here we have a dataset with a mixture of
# categorical and numerical features:
from sklearn.datasets import fetch_openml

X_adult, y_adult = fetch_openml("adult", version=2, return_X_y=True)

# Remove redundant and non-feature columns
X_adult = X_adult.drop(["education-num", "fnlwgt"], axis="columns")
X_adult.dtypes


# %%
# By setting `categorical_features="from_dtype"`, the gradient boosting classifier
# treats the columns with categorical dtypes as categorical features in the
# algorithm:
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score

X_train, X_test, y_train, y_test = train_test_split(X_adult, y_adult, random_state=0)
hist = HistGradientBoostingClassifier(categorical_features="from_dtype")

hist.fit(X_train, y_train)
y_decision = hist.decision_function(X_test)
print(f"ROC AUC score is {roc_auc_score(y_test, y_decision)}")

# %%
# Metadata Routing Support
# ------------------------
# Many meta-estimators and cross-validation routines now support metadata
# routing, which are listed in the :ref:`user guide
# <_metadata_routing_models>`. For instance, this is how you can do a nested
# cross-validation with sample weights and :class:`~model_selection.GroupKFold`:
import sklearn
from sklearn.metrics import get_scorer
from sklearn.datasets import make_regression
from sklearn.linear_model import Lasso
from sklearn.model_selection import GridSearchCV, cross_validate, GroupKFold
import numpy as np

# For now by default metadata routing is disabled, and need to be explicitly
# enabled.
sklearn.set_config(enable_metadata_routing=True)

n_samples = 100
X, y = make_regression(n_samples=n_samples, n_features=5, noise=0.5)
rng = np.random.RandomState(7)
groups = rng.randint(0, 10, size=n_samples)
sample_weights = rng.rand(n_samples)
estimator = Lasso().set_fit_request(sample_weight=True)
hyperparameter_grid = {"alpha": [0.1, 0.5, 1.0, 2.0]}
scoring_inner_cv = get_scorer("neg_mean_squared_error").set_score_request(
    sample_weight=True
)
inner_cv = GroupKFold(n_splits=5)

grid_search = GridSearchCV(
    estimator=estimator,
    param_grid=hyperparameter_grid,
    cv=inner_cv,
    scoring=scoring_inner_cv,
)

outer_cv = GroupKFold(n_splits=5)
scorers = {
    "mse": get_scorer("neg_mean_squared_error").set_score_request(sample_weight=True)
}
results = cross_validate(
    grid_search,
    X,
    y,
    cv=outer_cv,
    scoring=scorers,
    return_estimator=True,
    params={"sample_weight": sample_weights, "groups": groups},
)
print("cv error on test sets:", results["test_mse"])

# Setting the flag to the default `False` to avoid interference with other
# scripts.
sklearn.set_config(enable_metadata_routing=False)
