"""
===================================
Column Transformer with Mixed Types
===================================

This example illustrates how to apply different preprocessing and
feature extraction pipelines to different subsets of features,
using :class:`sklearn.compose.ColumnTransformer`.
This is particularly handy for the case of datasets that contain
heterogeneous data types, since we may want to scale the
numeric features and one-hot encode the categorical ones.

In this example, the numeric data is standard-scaled after
mean-imputation, while the categorical data is one-hot
encoded after imputing missing values with a new category
(``'missing'``).

Finally, the preprocessing pipeline is integrated in a
full prediction pipeline using :class:`sklearn.pipeline.Pipeline`,
together with a simple classification model.
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
               'scipy-2017-sklearn/091d371/notebooks/datasets/titanic3.csv')
data = pd.read_csv(titanic_url)

# We will train our classifier with the following features:
# Numeric Features:
# - age: float.
# - fare: float.
# Categorical Features:
# - embarked: categories encoded as strings {'C', 'S', 'Q'}.
# - sex: categories encoded as strings {'female', 'male'}.
# - pclass: ordinal integers {1, 2, 3}.
numeric_features = ['age', 'fare']
categorical_features = ['embarked', 'sex', 'pclass']

# Provisionally, use pd.fillna() to impute missing values for categorical
# features; SimpleImputer will eventually support strategy="constant".
data[categorical_features] = data[categorical_features].fillna(value='missing')

# We create the preprocessing pipelines for both numeric and categorical data.
numeric_transformer = make_pipeline(SimpleImputer(), StandardScaler())
categorical_transformer = CategoricalEncoder('onehot-dense',
                                             handle_unknown='ignore')

preprocessing_pl = make_column_transformer(
    (numeric_features, numeric_transformer),
    (categorical_features, categorical_transformer),
    remainder='drop'
)

# Append classifier to preprocessing pipeline.
# Now we have a full prediction pipeline.
clf = make_pipeline(preprocessing_pl, LogisticRegression())

X = data.drop('survived', axis=1)
y = data.survived.values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
                                                    shuffle=True)

clf.fit(X_train, y_train)
print("model score: %f" % clf.score(X_test, y_test))


###############################################################################
# Using the prediction pipeline in a grid search
###############################################################################
# Grid search can also be performed on the different preprocessing steps
# defined in the ``ColumnTransformer`` object, together with the classifier's
# hyperparameters as part of the ``Pipeline``.
# We will search for both the imputer strategy of the numeric preprocessing
# and the regularization parameter of the logistic regression using
# :class:`sklearn.model_selection.GridSearchCV`.


param_grid = {
    'columntransformer__pipeline__simpleimputer__strategy': ['mean', 'median'],
    'logisticregression__C': [0.1, 1.0, 1.0],
}

grid_search = GridSearchCV(clf, param_grid, cv=10, iid=False)
grid_search.fit(X_train, y_train)

print(("best logistic regression from grid search: %f"
       % grid_search.score(X_test, y_test)))
