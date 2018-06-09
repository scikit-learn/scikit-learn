"""
===================================
Column Transformer with Mixed Types
===================================

This example demonstrates how to use
:class:`sklearn.compose.ColumnTransformer` on a dataset consisting of
hereogeneous data types. A subset of numeric and categorical features
are selected and :class:`sklearn.compose.ColumnTransformer` is used to
apply two non-sequential preprocessing pipelines, one for each data
type:

* Numeric data: missing value imputation and standard scaling
* Categorical data: missing value imputation and one-hot encoding.
"""

# Author: Pedro Morales <part.morales@gmail.com>
#
# License: BSD 3 clause

from __future__ import print_function

import pandas as pd

from sklearn.compose import make_column_transformer
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, CategoricalEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV


# Read data from Titanic dataset.
titanic_url = ('https://raw.githubusercontent.com/amueller/'
               'scipy-2017-sklearn/master/notebooks/datasets/titanic3.csv')
data = pd.read_csv(titanic_url)

# We will train our classifier with the following features:
# Numeric Features:
# - age: float.
# - fare: float.
# Categorical Features:
# - embarked: categories encoded as strings {'C', 'S', 'Q'}.
# - sex: categories encoded as strings {'female', 'male'}.
# - plcass: categories encoded as ints {1, 2, 3}.
num_feats = ['age', 'fare']
cat_feats = ['embarked', 'sex', 'pclass']

# Provisionally, use pd.fillna() to impute missing values for categorical
# features; SimpleImputer will eventually support strategy="constant".
data.loc[:, cat_feats] = data.loc[:, cat_feats].fillna(value='missing')

# We create the preprocessing pipelines for both numeric and categorical data.
num_pl = make_pipeline(SimpleImputer(), StandardScaler())
cat_pl = CategoricalEncoder('onehot-dense', handle_unknown='ignore')

preprocessing_pl = make_column_transformer(
    (num_feats, num_pl),
    (cat_feats, cat_pl),
    remainder='drop'
)

# Append classifier to preprocessing pipeline.
# Now we have a full prediction pipeline.
clf_pl = make_pipeline(preprocessing_pl, LogisticRegression())

X = data.drop('survived', axis=1)
y = data.survived.values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
                                                    shuffle=True)

# Fit classifier
clf_pl.fit(X_train, y_train)
print("model score: %f" % clf_pl.score(X_test, y_test))


###############################################################################
# Using our pipeline in a grid search
###############################################################################
# We can also perform grid search on the different preprocessing steps
# defined in our ``ColumnTransformer``, together with the classifier's
# hyperparameters as part of a ``Pipeline``.
# ``ColumnTransformer`` integrates well with the rest of scikit-learn,
# in particular with ``GridSearchCV``.
# We will search for both the imputer strategy of the numeric preprocessing
# as for the regularization parameter of the logistic regression.


param_grid = {
    'columntransformer__pipeline__simpleimputer__strategy': ['mean', 'median'],
    'logisticregression__C': [0.1, 1.0, 1.0],
}

grid_search = GridSearchCV(clf_pl, param_grid, cv=10, iid=False)
grid_search.fit(X_train, y_train)

# We can finally fit the model using the best hyperparameters found
# and assess model performance on holdout set
print(("best logistic regression from grid search: %f"
       % grid_search.best_estimator_.score(X_test, y_test)))
