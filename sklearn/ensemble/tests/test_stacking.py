"""
Testing for the stacking ensemble module (sklearn.ensemble.stacking).
"""

# Author: Caio Oliveira
# License BSD 3 clause


from copy import deepcopy
import numpy as np
from sklearn.utils.testing import (assert_equal, assert_array_equal)
from sklearn.ensemble import (StackingTransformer, StackLayer, make_stack_layer)
from sklearn.linear_model import (RidgeClassifier, LinearRegression,
                                  LogisticRegression)
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC, LinearSVR
from sklearn import datasets
from sklearn.model_selection import (ParameterGrid, StratifiedKFold)

iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target

RANDOM_SEED = 8939

META_ESTIMATOR_PARAMS = {'cv': [2, StratifiedKFold()],
                         'n_jobs': [1, 2]}
META_ESTIMATOR_FIT_PARAMS = [{}, {"sample_weight": np.ones(y.shape)}]


def _check_estimator(estimator, **fit_params):
    # checks that we can fit_transform to the data
    Xt = estimator.fit_transform(X, y, **fit_params)

    # checks that we get a column vector
    assert_equal(Xt.ndim, 2)

    # checks that `fit` is available
    estimator.fit(X, y, **fit_params)

    # checks that we can transform the data after it's fitted
    Xt2 = estimator.transform(X)

    # checks that transformed data is always a column vector
    assert_equal(Xt2.ndim, 2)


def test_regression():
    # tests regression with various parameter settings

    meta_params = {'method': ['auto', 'predict']}
    meta_params.update(META_ESTIMATOR_PARAMS)

    regressors = [LinearRegression(), LinearSVR()]

    for reg in regressors:
        for params in ParameterGrid(meta_params):
            blended_reg = StackingTransformer(reg, **params)
            for fit_params in META_ESTIMATOR_FIT_PARAMS:
                _check_estimator(blended_reg, **fit_params)


def test_classification():
    # tests classification with various parameter settings

    testcases = [{'clf': RandomForestClassifier(random_state=RANDOM_SEED),
                  'extra_params': {'method': ['auto', 'predict',
                                              'predict_proba']}},
                 {'clf': LinearSVC(random_state=RANDOM_SEED),
                  'extra_params': {'method': ['auto', 'predict',
                                              'decision_function']}},
                 {'clf': RidgeClassifier(random_state=RANDOM_SEED),
                  'extra_params': {'method': ['auto', 'predict']}}]

    for testcase in testcases:
        clf = testcase['clf']

        meta_params = deepcopy(testcase['extra_params'])
        meta_params.update(META_ESTIMATOR_PARAMS)

        for params in ParameterGrid(meta_params):
            blended_clf = StackingTransformer(clf, **params)
            for fit_params in META_ESTIMATOR_FIT_PARAMS:
                _check_estimator(blended_clf, **fit_params)


def test_multi_output_classification():
    clf_base = RandomForestClassifier(random_state=RANDOM_SEED)
    clf = StackMetaEstimator(clf_base, method='predict_proba')
    X, y = datasets.make_multilabel_classification()
    Xt = clf.fit_transform(X[:-10], y[:-10])
    print(Xt)
    print(Xt.ndim)


STACK_LAYER_PARAMS = {'restack': [True, False],
                      'cv': [3, StratifiedKFold()],
                      'method': ['auto', 'predict', 'predict_proba']}


def _check_restack(X, Xorig):
    # checks that original data is appended to the rest of the features
    assert_array_equal(Xorig, X[:, -Xorig.shape[1]:])


def _check_layer(l):
    # check that we can fit_transform the data
    Xt = l.fit_transform(X, y)
    if l.restack:
        _check_restack(Xt, X)

    # check that we can transform the data
    Xt = l.transform(X)
    if l.restack:
        _check_restack(Xt, X)

    # check that `fit` is accessible
    l.fit(X, y)


def test_layer_regression():
    base_regs = [('lr', LinearRegression()),
                 ('svr', LinearSVR())]

    for params in ParameterGrid(STACK_LAYER_PARAMS):
        if params['method'] is 'predict_proba':
            continue
        # assert constructor
        reg_layer = StackLayer(base_estimators=base_regs, **params)
        _check_layer(reg_layer)


def test_layer_classification():
    base_clfs = [('rf1', RandomForestClassifier(random_state=RANDOM_SEED,
                                                criterion='gini')),
                 ('rf2', RandomForestClassifier(random_state=RANDOM_SEED,
                                                criterion='entropy'))]

    for params in ParameterGrid(STACK_LAYER_PARAMS):
        # assert constructor
        clf_layer = StackLayer(base_estimators=base_clfs, **params)
        _check_layer(clf_layer)


def test_layer_restack():
    base_estimators = [LinearRegression(), LinearRegression()]
    blended_layer = make_stack_layer(*base_estimators, restack=True,
                                     method='predict')
    X = np.random.rand(100, 3)
    y = np.random.rand(100)
    Xt = blended_layer.fit_transform(X, y)
    assert_array_equal(X, Xt[:, 2:])


def test_stacking_shortcuts():
    # Just test if the API generates the expected pipelines
    clf1 = RidgeClassifier(random_state=1)
    clf2 = LogisticRegression(random_state=1)
    clf3 = RandomForestClassifier(random_state=1)

    layer = make_stack_layer(clf1, clf2, clf3)

    assert_array_equal([x[1].base_estimator for x in layer.transformer_list],
                       [clf1, clf2, clf3])
