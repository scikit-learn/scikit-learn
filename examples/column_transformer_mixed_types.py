"""
===================================
Column Transformer with Mixed Types
===================================

This example demonstrates how to use
:class:`sklearn.compose.ColumnTransformer` on a dataset containing mixed
data types. We take a subset of numerical and categorical features
and use :class:`sklearn.compose.ColumnTransformer` to apply two non-sequential
preprocessing pipelines, one for each data type:

* Numerical data: missing value imputation, followed by standard scaling
* Categorical data: missing value imputation (with a custom transformer),
and then one-hot-encode.
"""

# Author: Pedro Morales <part.morales@gmail.com>
#
# License: BSD 3 clause

from __future__ import print_function

import numpy as np
import pandas as pd

from sklearn.base import TransformerMixin
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, CategoricalEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score


class CustomCategoryImputer(TransformerMixin):
    """Impute missing values for categorical features with a desired string.

    Parameters
    ----------
    imputer_value: string, default='missing'
        The placeholder for the missing value.

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array(['a', 'b', np.nan, 'c'], dtype=object)
    >>> imputer = CustomCategoryImputer(imputer_value='missing')
    >>> imputer.transform(X)
    array([['a'],
           ['b'],
           ['missing'],
           ['c']], dtype='<U7')
    """

    def __init__(self, imputer_value='missing'):
        self.imputer_value = imputer_value

    def fit(self, X, y=None):
        return self

    def transform(self, X):
        """Transform input filling missing values with imputer_value.

        Parameters
        ----------
        X: numpy array
            Categorical feature to transform
        """
        X_ = X.copy()

        # convert to DataFrame to handle mixed types
        if isinstance(X_, np.ndarray):
            X_ = pd.DataFrame(X_)

        mask = pd.notnull(X_)
        X_ = X_.where(mask, self.imputer_value)

        return X_.values.astype(str)


# Read data from Titanic dataset
titanic_url = ('https://raw.githubusercontent.com/amueller/'
               'scipy-2017-sklearn/master/notebooks/datasets/titanic3.csv')
data = pd.read_csv(titanic_url)

# We will train our classifier with the following features
num_feats = ['age', 'fare']
cat_feats = ['embarked', 'sex']

# We create the preprocessing pipelines for both numerical and categorical data
num_pl = make_pipeline(SimpleImputer(strategy='median'), StandardScaler())
cat_pl = make_pipeline(
    CustomCategoryImputer(),
    CategoricalEncoder('onehot-dense')
)

preprocessing_pl = ColumnTransformer(
    [
        ('num', num_pl, num_feats),
        ('cat', cat_pl, cat_feats)
    ],
    remainder='drop'
)

# Append classifier to preprocessing process. Use RF with default params
clf_pl = make_pipeline(preprocessing_pl, RandomForestClassifier())

X = data.drop('survived', axis=1)
y = data.survived.values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
                                                    shuffle=True)

# Fit classifier
clf_pl.fit(X_train, y_train)

# Assess model performance
y_pred = clf_pl.predict_proba(X_test)[:, 1]
score = roc_auc_score(y_test, y_pred)
print(score)
