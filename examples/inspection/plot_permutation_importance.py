"""
================================================================
Permutation Importance vs Random Forest Feature Importance (MDI)
================================================================

In this example, we will show on the titanic dataset how the impurity-based feature
importance (MDI, introduced by Breiman in [RF2001]_) of
:class:`~sklearn.ensemble.RandomForestClassifier` can give misleading results by
favoring high-cardinality features and we will give two alternatives to avoid the
issue.

In a nutshell, the impurity-based feature importance of random forests suffers from
being computed on statistics derived from the training dataset: the importances can be
high even for features that are not predictive of the target variable, as long as the
model has the capacity to use them to overfit. The effect is stronger the more unique
values the feature takes.

A first solution is to use :func:`~sklearn.inspection.permutation_importance` on test
data instead. Although this method is slower, it is not restricted to random forests and
does not suffer from the bias of MDI.

Another solution is to use the `unbiased_feature_importances_` attribute of random
forests, which leverages out-of-bag samples to correct the aforementioned bias. This
method was introduced by Li et al. in [UFI2020]_ and uses the samples that were not used
in the construction of each tree of the forest to modify the MDI.

.. rubric:: References

.. [RF2001] :doi:`"Random Forests" <10.1023/A:1010933404324>` L. Breiman, 2001
.. [UFI2020] :doi:`"Unbiased Measurement of Feature Importance in Tree-Based Methods"
   <10.1145/3429445>` Zhengze Zhou, Giles Hooker, 2020

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Data Loading and Feature Engineering
# ------------------------------------
# Let's use pandas to load a copy of the titanic dataset. The following shows
# how to apply separate preprocessing on numerical and categorical features.
#
# We further include two random variables that are not correlated in any way
# with the target variable (``survived``):
#
# - ``random_num`` is a high cardinality numerical variable (as many unique
#   values as records).
# - ``random_cat`` is a low cardinality categorical variable (3 possible
#   values).
import numpy as np

from sklearn.datasets import fetch_openml
from sklearn.model_selection import train_test_split

X, y = fetch_openml("titanic", version=1, as_frame=True, return_X_y=True)
rng = np.random.RandomState(seed=42)
X["random_cat"] = rng.randint(3, size=X.shape[0])
X["random_num"] = rng.randn(X.shape[0])

categorical_columns = ["pclass", "sex", "embarked", "random_cat"]
numerical_columns = ["age", "sibsp", "parch", "fare", "random_num"]

X = X[categorical_columns + numerical_columns]
X_train, X_test, y_train, y_test = train_test_split(X, y, stratify=y, random_state=42)

# %%
# We define a predictive model based on a random forest. Therefore, we will make
# the following preprocessing steps:
#
# - use :class:`~sklearn.preprocessing.OrdinalEncoder` to encode the
#   categorical features;
# - use :class:`~sklearn.impute.SimpleImputer` to fill missing values for
#   numerical features using a mean strategy.
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OrdinalEncoder

categorical_encoder = OrdinalEncoder(
    handle_unknown="use_encoded_value", unknown_value=-1, encoded_missing_value=-1
)
numerical_pipe = SimpleImputer(strategy="mean")

preprocessing = ColumnTransformer(
    [
        ("cat", categorical_encoder, categorical_columns),
        ("num", numerical_pipe, numerical_columns),
    ],
    verbose_feature_names_out=False,
)

rf = Pipeline(
    [
        ("preprocess", preprocessing),
        ("classifier", RandomForestClassifier(random_state=42, oob_score=True)),
    ]
)
rf.fit(X_train, y_train)

# %%
# Accuracy of the Model
# ---------------------
# Before inspecting the feature importances, it is important to check that
# the model predictive performance is high enough. Indeed, there would be little
# interest in inspecting the important features of a non-predictive model.
#
# By default, random forests subsample a part of the dataset to train each tree, a
# procedure known as bagging, leaving aside "out-of-bag" (oob) samples.
# These samples can be leveraged to compute an accuracy score independently of the
# training samples, when setting the parameter `oob_score = True`.
# This score should be close to the test score.

print(f"RF train accuracy: {rf.score(X_train, y_train):.3f}")
print(f"RF test accuracy: {rf.score(X_test, y_test):.3f}")
print(f"RF out-of-bag accuracy: {rf[-1].oob_score_:.3f}")

# %%
# Here, one can observe that the train accuracy is very high (the forest model
# has enough capacity to completely memorize the training set) but it can still
# generalize well enough to the test set thanks to the built-in bagging of
# random forests.
#
# It might be possible to trade some accuracy on the training set for a
# slightly better accuracy on the test set by limiting the capacity of the
# trees (for instance by setting ``min_samples_leaf=5`` or
# ``min_samples_leaf=10``) so as to limit overfitting while not introducing too
# much underfitting.
#
# However, let us keep our high capacity random forest model for now so that we can
# illustrate some pitfalls about feature importance on variables with many
# unique values.

# %%
# Tree's Feature Importance from Mean Decrease in Impurity (MDI)
# --------------------------------------------------------------
# The impurity-based feature importance ranks the numerical features to be the
# most important features. As a result, the non-predictive ``random_num``
# variable is ranked as one of the most important features!
#
# This problem stems from two limitations of impurity-based feature
# importances:
#
# - impurity-based importances are biased towards high cardinality features;
# - impurity-based importances are computed on training set statistics and
#   therefore do not reflect the ability of feature to be useful to make
#   predictions that generalize to the test set (when the model has enough
#   capacity).
#
# The bias towards high cardinality features explains why the `random_num` has
# a really large importance in comparison with `random_cat` while we would
# expect that both random features have a null importance.
#
# The fact that we use training set statistics explains why both the
# `random_num` and `random_cat` features have a non-null importance.
import pandas as pd

feature_names = rf[:-1].get_feature_names_out()

mdi_importances = pd.Series(
    rf[-1].feature_importances_, index=feature_names
).sort_values(ascending=True)

# %%
ax = mdi_importances.plot.barh()
ax.set_title("Random Forest Feature Importances (MDI)")
ax.set_xlabel("Decrease in impurity")
ax.figure.tight_layout()

# %%
# To avoid this issue, we can compute permutation importance instead. But we need to be
# careful as doing so on the train data will give wrong results.
# Indeed we can see that permutation importance on train data inflates the importance of
# every feature, even the random ones. Therefore one must be careful to use test data.
import matplotlib.pyplot as plt

from sklearn.inspection import permutation_importance

result_train = permutation_importance(
    rf, X_train, y_train, n_repeats=10, random_state=42, n_jobs=2
)
result_test = permutation_importance(
    rf, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
)

sorted_importances_idx = result_test.importances_mean.argsort()
importances_train = pd.DataFrame(
    result_train.importances[sorted_importances_idx].T,
    columns=X.columns[sorted_importances_idx],
)
importances_test = pd.DataFrame(
    result_test.importances[sorted_importances_idx].T,
    columns=X.columns[sorted_importances_idx],
)
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
importances_train.plot.box(vert=False, whis=10, ax=ax[0])
ax[0].set_title("Permutation Importances (train set)")
ax[0].axvline(x=0, color="k", linestyle="--")
ax[0].set_xlabel("Decrease in accuracy score")

importances_test.plot.box(vert=False, whis=10, ax=ax[1])
ax[1].set_title("Permutation Importances (test set)")
ax[1].axvline(x=0, color="k", linestyle="--")
ax[1].set_xlabel("Decrease in accuracy score")
fig.tight_layout()
# %%
# To see how this problem relates to overfitting, we can set `min_samples_leaf` at 20
# data points to reduce the overfitting of the model.
rf.set_params(classifier__min_samples_leaf=20).fit(X_train, y_train)

# %%
# Looking at the accuracy score on the training and testing set, we observe that
# the two metrics are very similar now. Therefore, our model is not overfitting
# anymore.
print(f"RF train accuracy: {rf.score(X_train, y_train):.3f}")
print(f"RF test accuracy: {rf.score(X_test, y_test):.3f}")

# %%
# We can see that our model is now much less reliant on uninformative features and
# therefore assigns lower importance to those. But we still have non zero importance
# values for completely random features when using train data only.
mdi_importances = pd.Series(
    rf[-1].feature_importances_, index=feature_names
).sort_values(ascending=True)

result_train = permutation_importance(
    rf, X_train, y_train, n_repeats=10, random_state=42, n_jobs=2
)

sorted_importances_idx = result_train.importances_mean.argsort()
importances_train = pd.DataFrame(
    result_train.importances[sorted_importances_idx].T,
    columns=X.columns[sorted_importances_idx],
)
fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(12, 6))
ax[0] = mdi_importances.plot.barh(ax=ax[0])
ax[0].set_xlabel("Decrease in impurity")
ax[0].set_title("Random Forest Feature Importances (MDI)")
importances_train.plot.box(vert=False, whis=10, ax=ax[1])
ax[1].set_title("Permutation Importances (train set)")
ax[1].axvline(x=0, color="k", linestyle="--")
ax[1].set_xlabel("Decrease in accuracy score")
fig.tight_layout()

# %%
# To completely ignore irrelevant features we should compute the permutation importance
# of ``rf`` on a held out test set.
# However when test samples are not available, or when permutation importance becomes
# too expensive to compute, there exists a modified version of the MDI,
# `unbiased_feature_importances_` available as soon as `oob_score` is set to `True`,
# that uses the out-of-bag samples of the trees to solve the bias problem.
ufi = rf[-1].unbiased_feature_importances_
mdi_importances = pd.Series(ufi, index=feature_names).sort_values(ascending=True)

# %%
ax = mdi_importances.plot.barh()
ax.set_title("Unbiased Feature Importances (UFI)")
ax.axvline(x=0, color="k", linestyle="--")
ax.set_xlabel("Decrease in impurity")
ax.figure.tight_layout()
# %%
# We can see that the random features have an importance of zero and the important
# features are ordered in the same way as with permutation importance. This method
# is much faster than permutation importances but is limited to random forests.
