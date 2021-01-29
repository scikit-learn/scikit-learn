"""
======================================================
Permutation Importance vs Mean Decrease Impurity (MDI)
======================================================

In this example, we will compare the impurity-based feature importance,
available by default in :class:`~sklearn.ensemble.RandomForestClassifier`,
with the permutation importance on the titanic dataset. This latter strategy
can be computed from two different manner: (i) using the out-of-bag samples
from the random-forest or (ii) by using a held-out dataset and the function
:func:`~sklearn.inspection.permutation_importance`. We will show that the
impurity-based feature importance can inflate the importance of numerical
features.

Furthermore, the impurity-based feature importance of random forests suffers
from being computed on statistics derived from the training dataset: the
importances can be high even for features that are not predictive of the target
variable, as long as the model has the capacity to use them to overfit.

This example shows how to use Permutation Importances as an alternative that
can mitigate those limitations.

.. topic:: References:

   [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32,
       2001. https://doi.org/10.1023/A:1010933404324
"""
print(__doc__)
import matplotlib.pyplot as plt
import sklearn
sklearn.set_config(display="diagram")

# %%
# Data Loading and Feature Engineering
# ------------------------------------
# We will use :func:`~sklearn.datasets.fetch_openml` to fetch the titanic
# dataset from OpenML and load it into a pandas dataframe.
#
# We further include two random variables that are not correlated in any way
# with the target variable (``survived``):
#
# - ``random_num`` is a high cardinality numerical variable (as many unique
#   values as records).
# - ``random_cat`` is a low cardinality categorical variable (3 possible
#   values).
from sklearn.datasets import fetch_openml
import numpy as np

X, y = fetch_openml("titanic", version=1, as_frame=True, return_X_y=True)
rng = np.random.RandomState(seed=42)
X['random_cat'] = rng.randint(3, size=X.shape[0])
X['random_num'] = rng.randn(X.shape[0])

# %%
from sklearn.model_selection import train_test_split
categorical_columns = ['pclass', 'sex', 'embarked', 'random_cat']
numerical_columns = ['age', 'sibsp', 'parch', 'fare', 'random_num']

X = X[categorical_columns + numerical_columns]

X_train, X_test, y_train, y_test = train_test_split(
    X, y, stratify=y, random_state=42)

# %%
# The following shows how to apply separate preprocessing on numerical and
# categorical features.
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder

categorical_encoder = OneHotEncoder(handle_unknown='ignore')
numerical_pipe = Pipeline([
    ('imputer', SimpleImputer(strategy='mean'))
])

preprocessing = ColumnTransformer(
    [('cat', categorical_encoder, categorical_columns),
     ('num', numerical_pipe, numerical_columns)])

rf = Pipeline([
    ('preprocess', preprocessing),
    ('classifier', RandomForestClassifier(random_state=42))
])
rf.fit(X_train, y_train)

# %%
# Accuracy of the Model
# ---------------------
# Prior to inspecting the feature importances, it is important to check that
# the model predictive performance is high enough. Indeed there would be little
# interest of inspecting the important features of a non-predictive model.
#
# Here one can observe that the train accuracy is very high (the forest model
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
# However let's keep our high capacity random forest model for now so as to
# illustrate some pitfalls with feature importance on variables with many
# unique values.
print("RF train accuracy: %0.3f" % rf.score(X_train, y_train))
print("RF test accuracy: %0.3f" % rf.score(X_test, y_test))


# %%
# Tree's Feature Importance from Mean Decrease in Impurity (MDI)
# --------------------------------------------------------------
# The impurity-based feature importance ranks the numerical features to be the
# most important features. As a result, the non-predictive ``random_num``
# variable is ranked the most important!
#
# This problem stems from two limitations of impurity-based feature
# importances:
#
# - impurity-based importances are biased towards high cardinality features;
# - impurity-based importances are computed on training set statistics and
#   therefore do not reflect the ability of feature to be useful to make
#   predictions that generalize to the test set (when the model has enough
#   capacity).
import pandas as pd

ohe = (rf.named_steps['preprocess']
         .named_transformers_['cat'])
feature_names = ohe.get_feature_names(input_features=categorical_columns)
feature_names = np.r_[feature_names, numerical_columns]

tree_feature_importances = pd.DataFrame(
    rf.named_steps['classifier'].importances_.importances.T,
    columns=feature_names)
# sort (reorder columns) the DataFrame for the plotting
tree_feature_importances = tree_feature_importances.reindex(
    tree_feature_importances.mean().sort_values().index,
    axis="columns")

ax = tree_feature_importances.plot.box(vert=False)
ax.set_title("Random Forest Feature Importances (MDI)")
_ = ax.set_xlabel("Impurity decrease")

# %%
# Alternative to MDI using Feature Permutation Importance
# -------------------------------------------------------
# The limitations of MDI pointed out in the previous section can be bypassed
# using an alternative strategy to estimate the feature importances. This
# strategy relies on monitoring the decrease (or not) of a given performance
# metric by randomly permutting the value of a given feature. In short, a
# predictive feature will negatively impact the score when it is randomly
# permuted while a non-predictive feature will not change the score.
#
# This feature permutation importance estimate can be computed in two different
# way: (i) by using the out-of-bag (OOB) samples in the ensemble to perform the
# permutation and the scoring or (ii) by manually splitting and handling a
# train and test set where the latter will be used with permutations.
#
# Feature Permutation Importance on Out-Of-Bag (OOB) samples
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Random-forest exposes a parameter `feature_importances` that allows to switch
# from the MDI to the permutation importance on the OOB samples. The parameter
# need to be set to `"permutation_oob"`.
rf = Pipeline(steps=[
    ("preprocess", preprocessing),
    ("classifier", RandomForestClassifier(
        feature_importances="permutation_oob", random_state=42))
]).fit(X_train, y_train)

# %%
# Once the forest has been train, the permutation importances have been
# estimated internally on the OOB samples. Thus, the fitted attribute
# `importances_` is now displaying the score decrease among all
# trees of the forest for each feature. Thus, we can plot this feature
# importances and compared it with the MDI estimates.
tree_feature_importances = pd.DataFrame(
    rf.named_steps['classifier'].importances_.importances.T,
    columns=feature_names)
# sort (reorder columns) the DataFrame for the plotting
tree_feature_importances = tree_feature_importances.reindex(
    tree_feature_importances.mean().sort_values().index,
    axis="columns")

ax = tree_feature_importances.plot.box(vert=False)
ax.set_title("Random Forest Feature Importances (OOB Permutation)")
_ = ax.set_xlabel("Accuracy decrease")

# %%
# With this strategy, the low cardinality categorical feature, ``sex`` is the
# most important. It gives both random features low importance, confirming that
# it avoids the limitations of MDI feature importances.
#
# Feature Permutation Importance on train-test sets
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In the previous section, we show how one can leverage the OOB samples to
# compute the permutation importance. However, this is also possible to use
# the same strategy but manipulating a train and a test sets.
#
# We illustrate such strategy by using the function
# :func:`~sklearn.inspection.permutation_importance`. Note that this way of
# computing the feature importance is model agnostic while the previous methods
# rely on the forest models.
from sklearn.inspection import permutation_importance

result = permutation_importance(rf, X_test, y_test, n_repeats=10,
                                random_state=42, n_jobs=2)
tree_feature_importances = pd.DataFrame(
    result.importances.T, columns=X_test.columns)
# sort (reorder columns) the DataFrame for the plotting
tree_feature_importances = tree_feature_importances.reindex(
    tree_feature_importances.mean().sort_values().index,
    axis="columns")

ax = tree_feature_importances.plot.box(vert=False)
ax.set_title("Permutation Importances (test set)")
_ = ax.set_xlabel("Accuracy decrease")

# %%
# As with the permutation importance using the OOB samples, the low cardinality
# categorical feature ``sex`` is the most important feature. Also note that
# both random features have very low importances (close to 0) as expected.
#
# It is also possible to compute the permutation importances on the training
# set. This reveals that ``random_num`` gets a significantly higher importance
# ranking than when computed on the test set. The difference between those two
# plots is a confirmation that the RF model has enough capacity to use that
# random numerical feature to overfit. You can further confirm this by
# re-running this example with constrained RF with `min_samples_leaf=10`.
result = permutation_importance(rf, X_train, y_train, n_repeats=10,
                                random_state=42, n_jobs=2)
tree_feature_importances = pd.DataFrame(
    result.importances.T, columns=X_test.columns)
# sort (reorder columns) the DataFrame for the plotting
tree_feature_importances = tree_feature_importances.reindex(
    tree_feature_importances.mean().sort_values().index,
    axis="columns")

ax = tree_feature_importances.plot.box(vert=False)
ax.set_title("Permutation Importances (train set)")
ax.set_xlabel("Accuracy decrease")
plt.show()

# %%
# Final words
# -----------
# As presented, the feature permutation importances can be computed either
# on the OOB samples or on separated datasets.
#
# While they are similar, it should be noted that the variations of the
# importances is estimated differently: the variance of the decrease of the
# score is estimated across the number of trees (i.e.``n_estimators``
# parameter) in the forest while it is estimated via the number of repeated
# permutation (i.e. `n_repeats`) in the other strategy.
#
# Therefore, using the permutation on the OOB samples could be interesting
# when a limited amount of data is at hand. Also, it might provide a faster way
# to evaluate the importances when setting the equivalence
# `n_repeats=n_estimators`.
#
# However, as shown in the previous plots, the permutation importances on the
# OOB will give a score on the random forest input features only. It means that
# this strategy does not allow to get information from original features,
# upstream from the random-forest. Computing the permutation importances on
# held-out train-test sets allows to apply or not a sequence of pre-processing
# and thus to know estimate the feature importances from any step of a
# machine-learning pipeline.
