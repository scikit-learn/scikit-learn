# Author: Pedro Morales <part.morales@gmail.com>
#
# License: BSD 3 clause

from __future__ import print_function

import numpy as np
import pandas as pd

from sklearn.compose import make_column_transformer, select_types
from sklearn.pipeline import make_pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler, CategoricalEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split


# Read data from Titanic dataset.
titanic_url = ('https://raw.githubusercontent.com/amueller/'
               'scipy-2017-sklearn/091d371/notebooks/datasets/titanic3.csv')

data = pd.read_csv(titanic_url)

# Provisionally, use pd.fillna() to impute missing values for categorical
# features; SimpleImputer will eventually support strategy="constant".
data[['embarked', 'sex', 'pclass']] = data[
    ['embarked', 'sex', 'pclass']].fillna(value='missing')


# We create the preprocessing pipelines for both numeric and categorical data.
numeric_transformer = make_pipeline(SimpleImputer(), StandardScaler())
categorical_transformer = CategoricalEncoder('onehot-dense',
                                             handle_unknown='ignore')

# We construct a column transformer, using dtype selectors
preprocessing_pl = make_column_transformer(
    (select_types([float]), numeric_transformer),
    (select_types([int, 'object']), categorical_transformer),
    remainder='drop'
)

clf = make_pipeline(preprocessing_pl, LogisticRegression())

np.random.seed(0)

# We will train our classifier with the following features:
# Numeric Features:
# - age: float.
# - fare: float.
# Categorical Features:
# - embarked: categories encoded as strings {'C', 'S', 'Q'}.
# - sex: categories encoded as strings {'female', 'male'}.
# - pclass: ordinal integers {1, 2, 3}.
X = data[['age', 'fare', 'embarked', 'sex', 'pclass']]
y = data.survived.values

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2,
                                                    shuffle=True)

clf.fit(X_train, y_train)
print("model score: %f" % clf.score(X_test, y_test))
