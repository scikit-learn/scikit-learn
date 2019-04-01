"""
==========================================================
Permutation Importance vs Random Forest Feature Importance
==========================================================

In this example, we will compare the
:class:`sklearn.ensemble.RandomForestClassifier` feature importance with the
permutation importance on the titanic dataset. We will show that the random
forest feature importance can inflate the importance of numerical features.

Furthermore, the built-in feature importance of random forests suffers from
being computed on statistics derived from the training set data: the
importances can be high even for features that are not predictive of the target
variable, as long as the model has the capacity to use them to overfit.

This example shows how to use Permutation Importances as an alternative that
can mitigate those limitations.

.. topic:: References:

   .. [1] Strobl, C., Boulesteix, AL., Zeileis, A. et al. BMC Bioinformatics
        (2007) 8: 25. https://doi.org/10.1186/1471-2105-8-25
"""
##############################################################################
# Training a Random Forest
# ------------------------
# A combination of numerical features and categorical features were used to
# train our random forest model:
print(__doc__)
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.model_inspection import permutation_importance
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder


##############################################################################
# Data Loading and Feature Engineering
# ------------------------------------
# Let's use pandas to load a copy of the titanic dataset. The following shows
# how to apply separate preprocessing on numerical and categorical features.
#
# We further include two random variables that are not correlated in any way
# with the target variable (``survived``):
# - ``random_num`` is a high cardinality numerical variable (as many unique
#   values as records);
# - ``random_cat`` is a low cardinality categorical variable (3 possible
#   values).
titanic_url = ('https://raw.githubusercontent.com/amueller/'
               'scipy-2017-sklearn/091d371/notebooks/datasets/titanic3.csv')
titanic = pd.read_csv(titanic_url)
titanic['random_cat'] = np.random.randint(3, size=titanic.shape[0])
titanic['random_num'] = np.random.randn(titanic.shape[0])

categorical_columns = ['pclass', 'sex', 'embarked', 'random_cat']
numerical_columns = ['age', 'sibsp', 'parch', 'fare', 'random_num']

data = titanic[categorical_columns + numerical_columns]
labels = titanic['survived']

X_train, X_test, y_train, y_test = train_test_split(
    data, labels, stratify=labels, random_state=42)

categorical_pipe = Pipeline([
    ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
    ('onehot', OneHotEncoder(handle_unknown='ignore'))
])
numerical_pipe = Pipeline([
    ('imputer', SimpleImputer(strategy='mean'))
])

preprocessing = ColumnTransformer(
    [('cat', categorical_pipe, categorical_columns),
     ('num', numerical_pipe, numerical_columns)])

rf = Pipeline([
    ('preprocess', preprocessing),
    ('classifier', RandomForestClassifier(n_estimators=100, min_samples_leaf=1,
                                          random_state=42))
])
rf.fit(X_train, y_train)

##############################################################################
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
# trees (for instance by setting ``max_samples_leaf=5`` or
# ``max_samples_leaf=10``) so as to limit overfitting while not introducing too
# much underfitting.
#
# However let's keep our high capacity random forest model for now so as to
# illustrate some pitfalls with feature importance on variables with many
# unique values.
print("RF train accuracy: %0.3f" % rf.score(X_train, y_train))
print("RF test accuracy: %0.3f" % rf.score(X_test, y_test))


##############################################################################
# Tree Based Feature Importance
# -----------------------------
# The tree based feature importance ranks the numerical features, to be the
# most important features. But furthermore, the non-predictive ``random_num``
# variable is ranked the most important!
#
# This problem stems from two limitations of RF feature importances:
# - RF importances are biased towards high cardinality features;
# - RF importances are computed on training set statistics and therefore do not
#   reflect the ability of feature to be useful to make predictions that
#   generalize to the test set (when the model as enough capacity).
ohe = (rf.named_steps['preprocess']
         .named_transformers_['cat']
         .named_steps['onehot'])
feature_names = []
for col, cats in zip(categorical_columns, ohe.categories_):
    for cat in cats:
        feature_names.append("{}_{}".format(col, cat))
feature_names = np.array(feature_names + numerical_columns)

tree_feature_importances = (
    rf.named_steps['classifier'].feature_importances_)
sorted_idx = tree_feature_importances.argsort()

y_ticks = np.arange(0, len(feature_names))
_, ax = plt.subplots(figsize=(10, 8))
ax.barh(y_ticks, tree_feature_importances[sorted_idx])
ax.set_yticklabels(feature_names[sorted_idx])
ax.set_yticks(y_ticks)
ax.set_title("Random Forest Feature Importances")


##############################################################################
# As an alternative, the permutation importances of ``rf`` are computed on a
# held out test set. This shows that the low cardinality categorical feature,
# ``sex`` is the most important feature.
#
# Also note that both random features have very low importances (close to 0) as
# expected.
permute_importance = permutation_importance(rf, X_test, y_test, n_rounds=30)
permute_importance_mean = np.mean(permute_importance, axis=-1)
sorted_idx = permute_importance_mean.argsort()

# sphinx_gallery_thumbnail_number = 2
_, ax = plt.subplots()
ax.boxplot(permute_importance[sorted_idx].T,
           vert=False, labels=X_test.columns[sorted_idx])
ax.set_title("Permutation Importances (test set)")

##############################################################################
# It is also possible to compute the permutation importances on the training
# set. This reveals that ``random_num`` gets a significantly higher importance
# ranking than when computed on the test set. The difference between those two
# plots is a confirmation that the RF model has enough capacity to use that
# random numerical feature to overfit. You can further confirm this by
# re-running this example with constrained RF with min_samples_leaf=10.
permute_importance = permutation_importance(rf, X_train, y_train, n_rounds=30)
permute_importance_mean = np.mean(permute_importance, axis=-1)
sorted_idx = permute_importance_mean.argsort()

# sphinx_gallery_thumbnail_number = 3
_, ax = plt.subplots()
ax.boxplot(permute_importance[sorted_idx].T,
           vert=False, labels=X_train.columns[sorted_idx])
ax.set_title("Permutation Importances (train set)")
plt.show()
