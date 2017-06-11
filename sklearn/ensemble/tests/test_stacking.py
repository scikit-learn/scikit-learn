"""
Testing for the stacking ensemble module (sklearn.ensemble.stacking).
"""

# Author: Caio Oliveira
# License BSD 3 clause


import numpy as np
from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_raises)
from sklearn.ensemble import (StackMetaEstimator, make_stack_layer)
from sklearn.linear_model import (RidgeClassifier, LinearRegression,
                                  LogisticRegression)
from sklearn.ensemble import RandomForestClassifier, BaggingClassifier
from sklearn.svm import SVC, SVR
from sklearn import datasets
from sklearn.model_selection import (ParameterGrid, StratifiedKFold)

iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target

RANDOM_SEED = 8939

META_ESTIMATOR_PARAMS = {'cv': [2, StratifiedKFold()],
                         'method': ['auto', 'predict', 'predict_proba'],
                         'n_jobs': [1, 2]}
META_ESTIMATOR_FIT_PARAMS = [{}, {"sample_weight": np.ones(y.shape)}]


def _check_estimator(estimator, **fit_params):
    # checks that we can fit_transform to the data
    Xt = estimator.fit_transform(X, y, **fit_params)

    # checks that we get a column vector
    assert_equal(Xt.ndim, 2)

    # checks that `fit` is not available
    assert_raises(NotImplementedError, estimator.fit, X, y)

    # checks that we can transform the data after it's fitted
    Xt2 = estimator.transform(X)

    # checks that transformed data is always a column vector
    assert_equal(Xt2.ndim, 2)


def test_regression():
    # tests regression with various parameter settings

    regressors = [LinearRegression(), SVR()]

    for reg in regressors:
        for params in ParameterGrid(META_ESTIMATOR_PARAMS):
            if params['method'] is 'predict_proba':
                # no need to test this, as it's related to classification
                continue
            blended_reg = StackMetaEstimator(reg, **params)
            for fit_params in META_ESTIMATOR_FIT_PARAMS:
                _check_estimator(blended_reg, **fit_params)


def test_classification():
    # tests classification with various parameter settings

    classifiers_with_proba = [RandomForestClassifier(random_state=RANDOM_SEED),
                              BaggingClassifier(RidgeClassifier(),
                                                random_state=RANDOM_SEED)]
    classifiers_without_proba = [RidgeClassifier(random_state=RANDOM_SEED),
                                 SVC(random_state=RANDOM_SEED)]

    for clf in classifiers_with_proba:
        for params in ParameterGrid(META_ESTIMATOR_PARAMS):
            blended_clf = StackMetaEstimator(clf, **params)
            for fit_params in META_ESTIMATOR_FIT_PARAMS:
                _check_estimator(blended_clf, **fit_params)

    # test method='auto' for classifiers without 'predict_proba'
    for clf in classifiers_without_proba:
        clf = StackMetaEstimator(clf, method='auto')
        _check_estimator(blended_clf, **fit_params)


def test_stacking_shortcuts():
    # Just test if the API generates the expected pipelines
    clf1 = RidgeClassifier(random_state=1)
    clf2 = LogisticRegression(random_state=1)
    clf3 = RandomForestClassifier(random_state=1)

    layer = make_stack_layer(clf1, clf2, clf3)

    assert_array_equal([x[1].base_estimator for x in layer.transformer_list],
                       [clf1, clf2, clf3])


def test_layer_regression():
    pass


def test_restack():
    base_estimators = [LinearRegression(), LinearRegression()]
    blended_layer = make_stack_layer(*base_estimators, restack=True,
                                     method='predict')
    X = np.random.rand(100, 3)
    y = np.random.rand(100)
    Xt = blended_layer.fit_transform(X, y)
    assert_array_equal(X, Xt[:, 2:])
