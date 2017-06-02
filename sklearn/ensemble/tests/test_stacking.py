"""
Testing for the stacking ensemble module (sklearn.ensemble.stacking).
"""

# Author: Caio Oliveira
# License BSD 3 clause


import numpy as np
from sklearn.utils.testing import (assert_equal, assert_array_equal,
                                   assert_raises)
from sklearn.ensemble import (BlendedEstimator, make_stack_layer,
                              stack_estimators)
from sklearn.linear_model import (LogisticRegression, RidgeClassifier,
                                  LinearRegression)
from sklearn.ensemble import RandomForestClassifier, BaggingClassifier
from sklearn.svm import SVC
from sklearn import datasets
from sklearn.model_selection import (ParameterGrid, StratifiedKFold)

iris = datasets.load_iris()
X, y = iris.data[:, 1:3], iris.target


def test_transformer_init():
    # Tests if attributes are set correctly

    base_clf = RidgeClassifier(random_state=1)
    blended = BlendedEstimator(base_clf)

    assert_equal(blended.base_estimator, base_clf)
    assert_equal(blended.cv, 3)


def test_stacking_api():
    # Just test if the API generates the expected pipelines
    clf1 = RidgeClassifier(random_state=1)
    clf2 = LogisticRegression(random_state=1)
    clf3 = RandomForestClassifier(random_state=1)
    clf4 = BaggingClassifier(random_state=1)
    clf5 = SVC(random_state=1)

    layer = make_stack_layer(clf1, clf2, clf3)

    assert_array_equal([x[1].base_estimator for x in layer.transformer_list],
                       [clf1, clf2, clf3])

    clf_matrix = [[clf1, clf2], [clf3, clf4]]
    final_clf = stack_estimators(clf_matrix, clf5)

    matrix_from_pipeline = [[x[1].base_estimator
                             for x in y[1].transformer_list]
                            for y in final_clf.steps[:-1]]

    assert_array_equal(clf_matrix, matrix_from_pipeline)

    assert_equal(final_clf.steps[-1][1], clf5)


def test_classification():
    # tests classification with various parameter settings

    sample_weight = np.ones(y.shape)

    # grid with some classifiers that don't have `predict_proba`
    grid1 = ParameterGrid({'base_estimator':
                           [RidgeClassifier(random_state=1),
                            LogisticRegression(random_state=1),
                            LinearRegression(),
                            RandomForestClassifier(random_state=1),
                            SVC(random_state=1)],
                           'cv': [2, StratifiedKFold()],
                           'method': ['auto', 'predict'],
                           'n_jobs': [1, 2],
                           'fit_params': [{},
                                          {'sample_weight': sample_weight}]})

    # grid with classifiers that have both predict and predict_proba
    grid2 = ParameterGrid({'base_estimator':
                           [RandomForestClassifier(random_state=1),
                            BaggingClassifier(RidgeClassifier(),
                                              random_state=1)],
                           'cv': [2, StratifiedKFold()],
                           'method': ['auto', 'predict', 'predict_proba'],
                           'n_jobs': [1, 2],
                           'fit_params': [{},
                                          {'sample_weight': sample_weight}]})

    params_list = list(grid1)+list(grid2)

    for params in params_list:
        fit_params = params.pop('fit_params')
        clf = BlendedEstimator(**params)
        Xt = clf.fit_transform(X, y, **fit_params)

        # checks that we get a column vector
        assert_equal(Xt.ndim, 2)

        assert_raises(NotImplementedError, clf.fit, X, y, **fit_params)


def test_restacking():
    base_estimators = [LinearRegression(), LinearRegression()]
    blended_layer = make_stack_layer(*base_estimators, restacking=True,
                                     method='predict')
    X = np.random.rand(100, 3)
    y = np.random.rand(100)
    Xt = blended_layer.fit_transform(X, y)
    assert_array_equal(X, Xt[:, 2:])
