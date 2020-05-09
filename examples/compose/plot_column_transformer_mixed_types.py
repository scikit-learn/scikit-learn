"""
===================================
Column Transformer with Mixed Types
===================================

This example illustrates how to apply different preprocessing and feature
extraction pipelines to different subsets of features, using
:class:`sklearn.compose.ColumnTransformer`. This is particularly handy for the
case of datasets that contain heterogeneous data types, since we may want to
scale the numeric features and one-hot encode the categorical ones.

In this example, the numeric data is standard-scaled after mean-imputation,
while the categorical data is one-hot encoded after imputing missing values
with a new category (``'missing'``).

In addition, we show two different ways to dispatch the columns to the
particular pre-processor: by column names and by column data types.

Finally, the preprocessing pipeline is integrated in a full prediction pipeline
using :class:`sklearn.pipeline.Pipeline`, together with a simple classification
model.
"""

# Author: Pedro Morales <part.morales@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from sklearn.compose import ColumnTransformer
from sklearn.datasets import fetch_openml
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV

np.random.seed(0)

# Load data from https://www.openml.org/d/40945
X, y = fetch_openml("titanic", version=1, as_frame=True, return_X_y=True)

# Alternatively X and y can be obtained directly from the frame attribute:
# X = titanic.frame.drop('survived', axis=1)
# y = titanic.frame['survived']

###############################################################################
# Use ``ColumnTransformer`` by selecting column by names
###############################################################################
# We will train our classifier with the following features:
#
# Numeric Features:
#
# * ``age``: float;
# * ``fare``: float.
#
# Categorical Features:
#
# * ``embarked``: categories encoded as strings ``{'C', 'S', 'Q'}``;
# * ``sex``: categories encoded as strings ``{'female', 'male'}``;
# * ``pclass``: ordinal integers ``{1, 2, 3}``.
#
# We create the preprocessing pipelines for both numeric and categorical data.

numeric_features = ['age', 'fare']
numeric_transformer = Pipeline(steps=[
    ('imputer', SimpleImputer(strategy='median')),
    ('scaler', StandardScaler())])

categorical_features = ['embarked', 'sex', 'pclass']
categorical_transformer = Pipeline(steps=[
    ('imputer', SimpleImputer(strategy='constant', fill_value='missing')),
    ('onehot', OneHotEncoder(handle_unknown='ignore'))])

preprocessor = ColumnTransformer(
    transformers=[
        ('num', numeric_transformer, numeric_features),
        ('cat', categorical_transformer, categorical_features)])

# Append classifier to preprocessing pipeline.
# Now we have a full prediction pipeline.
clf = Pipeline(steps=[('preprocessor', preprocessor),
                      ('classifier', LogisticRegression())])

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

clf.fit(X_train, y_train)
print("model score: %.3f" % clf.score(X_test, y_test))

##############################################################################
# HTML representation of ``Pipeline``
###############################################################################
# When the ``Pipeline`` is printed out in a jupyter notebook an HTML
# representation of the estimator is displayed as follows:
from sklearn import set_config
set_config(display='diagram')
clf

###############################################################################
# Use ``ColumnTransformer`` by selecting column by data types
###############################################################################
# When dealing with a cleaned dataset, the preprocessing can be automatic by
# using the data types of the column to decide whether to treat a column as a
# numerical or categorical feature.
# :func:`sklearn.compose.make_column_selector` gives this possibility.
# First, let's only select a subset of columns to simplify our
# example.

subset_feature = ['embarked', 'sex', 'pclass', 'age', 'fare']
X = X[subset_feature]

###############################################################################
# Then, we introspect the information regarding each column data type.

X.info()

###############################################################################
# We can observe that the `embarked` and `sex` columns were tagged as
# `category` columns when loading the data with ``fetch_openml``. Therefore, we
# can use this information to dispatch the categorical columns to the
# ``categorical_transformer`` and the remaining columns to the
# ``numerical_transformer``.

###############################################################################
# .. note:: In practice, you will have to handle yourself the column data type.
#    If you want some columns to be considered as `category`, you will have to
#    convert them into categorical columns. If you are using pandas, you can
#    refer to their documentation regarding `Categorical data
#    <https://pandas.pydata.org/pandas-docs/stable/user_guide/categorical.html>`_.

from sklearn.compose import make_column_selector as selector

preprocessor = ColumnTransformer(transformers=[
    ('num', numeric_transformer, selector(dtype_exclude="category")),
    ('cat', categorical_transformer, selector(dtype_include="category"))
])

# Reproduce the identical fit/score process
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

clf.fit(X_train, y_train)
print("model score: %.3f" % clf.score(X_test, y_test))

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
    'preprocessor__num__imputer__strategy': ['mean', 'median'],
    'classifier__C': [0.1, 1.0, 10, 100],
}

grid_search = GridSearchCV(clf, param_grid, cv=10)
grid_search.fit(X_train, y_train)

print(("best logistic regression from grid search: %.3f"
       % grid_search.score(X_test, y_test)))
