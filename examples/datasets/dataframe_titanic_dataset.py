#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
=========================================================
The Titanic Dataset
=========================================================
This dataset is based on the Titanic Passenger List edited by Michael
A. Findlay, originally published in Eaton & Haas (1994) Titanic:
Triumph and Tragedy, Patrick Stephens Ltd, and expanded with the help
of the internet community.  The original HTML files were obtained by
Philip Hind (1999).

See `here
<http://biostat.mc.vanderbilt.edu/wiki/pub/Main/DataSets/titanic3info.txt>`_
for more information about this dataset.

We proceed to find the best parameters (C, gamma) to train a SVM with a RBF
based kernel via grid search over a pipeline estimator which demonstrates the
use of the DataFrameMapper object for feature selection, processing and
normalization.

"""

import numpy as np
import pandas as pd
from sklearn.datasets import titanic
from sklearn.preprocessing import DataFrameMapper
from sklearn.preprocessing import StandardScaler, LabelBinarizer, MinMaxScaler
from sklearn.svm import SVC
from sklearn import pipeline
from sklearn.metrics import accuracy_score, Scorer
from sklearn.cross_validation import ShuffleSplit
from sklearn.grid_search import GridSearchCV


def preprocess(df):
    """
    Simple preprocessor for Titanic Dataset.

    Parameters
    ----------
    df: Titanic Dataframe

    Notes
    ------
    Fills in missing fare values with median value for the class
    Fills in missing embarked values with most common one (S=1)
    Fills in missing ages to static value of 0
    Deriving new title feature (e.g. 'Mrs.') from full name feature
    Deriving new cabin letter feature (e.g. 'T') from cabin feature
    """
    df.fare[pd.isnull(df.fare)] = 0
    df.fare[df.fare == 0] = df.pclass[df.fare == 0].apply(
        lambda x: df.fare[df.pclass == x].median()
    )
    df.embarked = df.embarked.apply(
        lambda x: dict(C=0, S=1, Q=2).get(x, 1)
    )
    df.age[pd.isnull(df.age)] = 0
    df["title"] = df["name"].apply(
        lambda x: x.split(",")[1].split(".")[0].strip()
    )
    df["cabin_letter"] = df["cabin"].apply(
        lambda x: None if pd.isnull(x) else x[0]
    )
    return df

if __name__ == "__main__":

    print(__doc__)

    df = preprocess(titanic.load())

    mapper = DataFrameMapper([
        ('age', StandardScaler()),
        ('pclass', LabelBinarizer()),
        ('title', LabelBinarizer()),
        ('sex', LabelBinarizer()),
        ('parch', LabelBinarizer()),
        ('sibsp', LabelBinarizer()),
        ('cabin', LabelBinarizer()),
        # ('cabin_letter', LabelBinarizer()),
        ('embarked', LabelBinarizer()),
        ('ticket', LabelBinarizer()),
        ('fare', StandardScaler()),
    ], sub=None)

    titanic_estimator = pipeline.Pipeline([
        ('featurize', mapper),
        ('scale', MinMaxScaler(feature_range=(0, 1))),
        ('cl', SVC(kernel="rbf", tol=1e-6))
    ])

    C_range = 10.0 ** np.arange(0, 4)
    gamma_range = 10.0 ** np.arange(-4, 0)
    params = dict(cl__class_weight=["auto"], cl__C=C_range,
                  cl__gamma=gamma_range, featurize__sub=[None])
    cross_validator_iter = ShuffleSplit(len(df), n_iter=10, test_size=0.10,
                                        random_state=0, indices=False)
    accuracy_scorer = Scorer(accuracy_score)

    grid_search = GridSearchCV(titanic_estimator, param_grid=params, verbose=1,
                               cv=cross_validator_iter, n_jobs=-1,
                               scoring=accuracy_scorer, refit=True)
    grid_search.fit(df, df.survived)

    print(grid_search.best_params_)
    print(grid_search.best_score_)
