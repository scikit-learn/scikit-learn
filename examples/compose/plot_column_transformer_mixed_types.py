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

import pandas as pd
import numpy as np

from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split, GridSearchCV

np.random.seed(0)

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
pipeline = Pipeline(steps=[('preprocessor', preprocessor),
                    ('classifier', LogisticRegression())])

X = data.drop('survived', axis=1)
y = data['survived']

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

pipeline.fit(X_train, y_train)
print("model score: %.3f" % pipeline.score(X_test, y_test))


###############################################################################
# Inspecting the coefficients values of the classifier
###############################################################################
# The coefficients of the final classification step of the pipeline gives an
# idea how each feature impacts the likelihood of survival assuming that the
# usual linear model assumptions hold (uncorrelated features, linear
# separability, homoschedastic errors...) which we do not verify in this
# example.
#
# To get error bars we perform cross-validation and compute the mean and
# standard deviation for each coefficient accross CV splits. Because we use a
# standard scaler on the numerical features, the coefficient weights gives us
# an idea on how much the log odds of surviving are impacted by a change in
# this dimension contrasted to the mean. Note that the categorical features
# here are overspecified which makes it slightly harder to interpret because of
# the information redundancy.
#
# We can see that the linear model coefficients are in agreement with the
# historical reports: people in higher classes and therefore in the upper decks
# were the first to reach the lifeboats, and often, priority was given to women
# and children.
#
# Note that conditionned on the "pclass_x" one-hot features, the "fare"
# numerical feature does not seem to be significantly predictive. If we drop
# the "pclass" feature, then higher "fare" values would appear significantly
# correlated with a higher likelihood of survival as the "fare" and "pclass"
# features have a strong statistical dependency.

import matplotlib.pyplot as plt
from sklearn.model_selection import cross_validate
from sklearn.model_selection import StratifiedShuffleSplit

cv = StratifiedShuffleSplit(n_splits=20, test_size=0.25, random_state=42)
cv_results = cross_validate(pipeline, X_train, y_train, cv=cv,
                            return_estimator=True)
cv_coefs = np.concatenate([cv_pipeline.named_steps["classifier"].coef_
                           for cv_pipeline in cv_results["estimator"]])
fig, ax = plt.subplots()
ax.barh(pipeline.named_steps["classifier"].feature_names_in_,
        cv_coefs.mean(axis=0), xerr=cv_coefs.std(axis=0))
plt.tight_layout()
plt.show()


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

grid_search = GridSearchCV(pipeline, param_grid, cv=10)
grid_search.fit(X_train, y_train)

print(("best logistic regression from grid search: %.3f"
       % grid_search.score(X_test, y_test)))
