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
# ------------------------------------------------------------------
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
