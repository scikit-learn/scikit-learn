"""
===================================
Column Transformer with Mixed Types
===================================

This example demonstrates how to use
:class:`sklearn.compose.ColumnTransformer` on a dataset containing categorical 
and continuous data.  We use the Credit Approval Dataset (UCI ML Repository -
https://archive.ics.uci.edu/ml/machine-learning-databases/credit-screening).
We take a subset of numerical and categorical features and use 
:class:`sklearn.compose.ColumnTransformer` to standard-scale the former and 
one-hot-encode the latter.
The choice of features is not particularly helpful, but serves to illustrate
the technique.
"""

# Author: Pedro Morales <part.morales@gmail.com>
# 
# License: BSD 3 clause

import io
import numpy as np
from sklearn.externals.six.moves.urllib.request import urlopen

def get_clean_data():
    """
    Retrieve data from UCI credit screening database and store it in a numpy array
    """
    url = ('https://archive.ics.uci.edu/ml/machine-learning-databases/'
           'credit-screening/crx.data')

    # Read and decode the data from the url
    f = io.StringIO(urlopen(url).read().decode('utf-8'))

    # Drop rows with missing values (in this dataset NAs are '?')
    f = [row for row in f if not '?' in row]

    # Load data as numpy array
    data = np.loadtxt(f, dtype='object', delimiter=',')

    # Parse numerical data 
    for col in [1, 2, 7, 10, 13, 14]:
        data[:, col] = data[:, col].astype(float)

    return data


data = get_clean_data()

from sklearn.preprocessing import LabelBinarizer
from sklearn.model_selection import train_test_split

# Split train/test sets. The last column of the array contains the class labels
X_train, X_test, y_train, y_test = train_test_split(
    data[:, :-1],
    LabelBinarizer().fit_transform(data[:, -1]),
    test_size=0.2, 
    shuffle=True,
    random_state=0
)

from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import CategoricalEncoder, StandardScaler

# For this example, we will train our classifier with a subset of columns: 
# (numerical), and [0, 5] (categorical). 
# We standard-scale the numerical features and one-hot-encode the categorical ones.
ct = ColumnTransformer(
    [
        ('sc', StandardScaler(), [1, 2]),
        ('ohe', CategoricalEncoder('onehot-dense'), [0, 5])
    ],
    remainder='drop'
)

from sklearn.pipeline import make_pipeline
from sklearn.ensemble import RandomForestClassifier

# We can put our column transformer into a pipeline object with the classifier
pl = make_pipeline(ct, RandomForestClassifier())


pl.fit(X_train, y_train)
y_pred = pl.predict_proba(X_test)[:, 1]

from sklearn.metrics import roc_auc_score

score = roc_auc_score(y_test, y_pred)
print(score)
