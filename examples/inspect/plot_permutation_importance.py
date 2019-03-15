"""
==========================================================
Permutation Importance vs Random Forest Feature Importance
==========================================================

The permutation importance of a feature is how much a model's score decreases
when the feature's values are permuted. [1]_ In this context, the more
important features will cause the model's score to decrease the most. If a
feature's importance is negative, the model improved when the feature is
permutated. This suggests that the model would benefit from removing the
feature.

.. topic:: References:

   .. [1] Breiman, L. Machine Learning (2001) 45: 5.
        https://doi.org/10.1023/A:1010933404324
   .. [2] Strobl, C., Boulesteix, AL., Zeileis, A. et al. BMC Bioinformatics
        (2007) 8: 25. https://doi.org/10.1186/1471-2105-8-25
"""
print(__doc__)
##############################################################################
# In this example, the :class:`sklearn.ensemble.RandomForestRegressor`'s feature
# importance is compared with the permutation importance using the titantic
# dataset. A combination of numerical features and categorical features were
# used to train the initial model, ``rf_init``:
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.inspect import permutation_importance
from sklearn.compose import ColumnTransformer
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, OneHotEncoder

titanic_url = ('https://raw.githubusercontent.com/amueller/'
               'scipy-2017-sklearn/091d371/notebooks/datasets/titanic3.csv')
titanic = pd.read_csv(titanic_url)

categorical_columns = ['pclass', 'sex', 'embarked']
numerical_columns = ['age', 'sibsp', 'parch', 'fare']

data = titanic[categorical_columns + numerical_columns]
labels = titanic['survived']

X_train, X_test, y_train, y_test = train_test_split(
    data, labels, stratify=labels, random_state=42)

categorical_pipe = Pipeline([
    ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
    ('onehot', OneHotEncoder(handle_unknown='ignore'))
])
numerical_pipe = Pipeline([
    ('scaler', StandardScaler()),
    ('imputer', SimpleImputer(strategy='mean'))
])

preprocessing = ColumnTransformer(
    [('cat', categorical_pipe, categorical_columns),
     ('num', numerical_pipe, numerical_columns)])

rf_init = Pipeline([
    ('preprocess', preprocessing),
    ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
])
rf_init.fit(X_train, y_train)

##############################################################################
# The tree based feature importance ranks the numerical features, age and fare,
# to be the most important features:
ohe = (rf_init.named_steps['preprocess']
              .named_transformers_['cat']
              .named_steps['onehot'])
feature_names = []
for col, cats in zip(categorical_columns, ohe.categories_):
    for cat in cats:
        feature_names.append("{}_{}".format(col, cat))
feature_names = np.array(feature_names + numerical_columns)

tree_feature_importances = (
    rf_init.named_steps['classifier'].feature_importances_)
sorted_idx = tree_feature_importances.argsort()

y_ticks = np.arange(0, len(feature_names))
_, ax = plt.subplots(figsize=(10, 8))
ax.barh(y_ticks, tree_feature_importances[sorted_idx])
ax.set_yticklabels(feature_names[sorted_idx])
ax.set_yticks(y_ticks)
ax.set_title("Random Forest Feature importance")

##############################################################################
# The corresponding test score for, ``rf_init``, is:
print("rf_init test score is", rf_init.score(X_test, y_test))

##############################################################################
# Next, a model, ``rf_new``, is trained with the numerical features removed
# and its test score is calcuated.
preprocessing_cat = ColumnTransformer(
    [('cat', categorical_pipe, categorical_columns)])
rf_new = Pipeline([
    ('preprocess', preprocessing_cat),
    ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
])
print("rf_new test score is",
      rf_new.fit(X_train, y_train).score(X_test, y_test))

##############################################################################
# When the "important" numerical features are removed, the score increases! A
# decision tree's feature importance is known to be inflated for continuous
# variables. [2]_ As an alternatie, the permutation permutation of
# ``rf_init`` is computed, which shows that the categorical feature, ``sex``
# is the most important feature:
permute_importance = permutation_importance(rf_init, X_train, y_train,
                                            n_bootstrap=10)
permute_importance_mean = np.mean(permute_importance, axis=-1)
sorted_idx = permute_importance_mean.argsort()

# sphinx_gallery_thumbnail_number = 2
_, ax = plt.subplots()
ax.boxplot(permute_importance[sorted_idx].T,
           vert=False, labels=X_train.columns[sorted_idx])
ax.set_title("Permutation Importance")

plt.show()
