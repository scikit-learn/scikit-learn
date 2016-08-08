# License: BSD 3 clause

import numpy as np

from sklearn import datasets
from sklearn.ensemble import GradientBoostingClassifierCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import GradientBoostingRegressorCV
from sklearn.ensemble import GradientBoostingRegressor

from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split

from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_greater
from sklearn.utils.testing import assert_less
from sklearn.utils.testing import assert_array_equal

rng = np.random.RandomState(0)

iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
X_clf = iris.data[perm]
y_clf = iris.target[perm]

boston = datasets.load_boston()
perm = rng.permutation(boston.target.size)
X_reg = boston.data[perm]
y_reg = boston.target[perm]


def test_no_param():
    X_train, X_test, y_train, y_test = train_test_split(X_clf, y_clf,
                                                        random_state=42)
    gbccv = GradientBoostingClassifierCV(max_iterations=1000,
                                         score_precision=4,
                                         random_state=42)
    gbccv.fit(X_train, y_train)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbccv.best_params_['n_estimators'], 13)
    assert_greater(gbccv.best_estimator_.score(X_test, y_test), 0.9)

    X_train, X_test, y_train, y_test = train_test_split(X_reg, y_reg,
                                                        random_state=42)
    gbrcv = GradientBoostingRegressorCV(max_iterations=1000,
                                        score_precision=4,
                                        random_state=42)
    gbrcv.fit(X_train, y_train)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbrcv.best_params_['n_estimators'], 87)
    assert_greater(gbrcv.best_estimator_.score(X_test, y_test), 0.85)


def test_single_param():
    X_train, X_test, y_train, y_test = train_test_split(X_clf, y_clf,
                                                        random_state=42)
    gbccv = GradientBoostingClassifierCV(max_iterations=1000,
                                         n_stop_rounds=10,
                                         learning_rate=0.1,
                                         max_depth=3,
                                         score_precision=4,
                                         random_state=42)
    gbccv.fit(X_train, y_train)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbccv.best_params_['n_estimators'], 13)
    assert_greater(gbccv.best_estimator_.score(X_test, y_test), 0.9)

    X_train, X_test, y_train, y_test = train_test_split(X_reg, y_reg,
                                                        random_state=42)
    gbrcv = GradientBoostingRegressorCV(max_iterations=1000,
                                        n_stop_rounds=10,
                                        learning_rate=0.1,
                                        max_depth=3,
                                        score_precision=4,
                                        random_state=42)
    gbrcv.fit(X_train, y_train)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbrcv.best_params_['n_estimators'], 87)
    assert_greater(gbrcv.best_estimator_.score(X_test, y_test), 0.80)


def test_param_grid():
    X_train, X_test, y_train, y_test = train_test_split(X_clf, y_clf,
                                                        random_state=42)
    gbccv = GradientBoostingClassifierCV(max_iterations=1000,
                                         learning_rate=[0.1, 0.3],
                                         max_depth=[3, 4],
                                         score_precision=4,
                                         random_state=42)
    gbccv.fit(X_train, y_train)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbccv.best_params_['n_estimators'], 10)
    assert_greater(gbccv.best_estimator_.score(X_test, y_test), 0.9)

    X_train, X_test, y_train, y_test = train_test_split(X_reg, y_reg,
                                                        random_state=42)
    gbrcv = GradientBoostingRegressorCV(max_iterations=1000,
                                        learning_rate=[0.1, 0.3],
                                        max_depth=[3, 4],
                                        score_precision=4,
                                        random_state=42)
    gbrcv.fit(X_train, y_train)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbrcv.best_params_['n_estimators'], 56)
    assert_greater(gbrcv.best_estimator_.score(X_test, y_test), 0.8)


def test_predict():
    X_train, X_val, y_train, y_val = train_test_split(X_clf, y_clf,
                                                      random_state=42)
    params = {'learning_rate': [0.1, 0.3], 'max_depth': [3, 4]}
    gbccv = GradientBoostingClassifierCV(score_precision=4, random_state=42,
                                         **params)
    gbccv.fit(X_train, y_train)

    gs_gbc = GridSearchCV(GradientBoostingClassifier(random_state=42), params)
    gs_gbc.fit(X_train, y_train)

    assert_array_equal(gbccv.predict(X_val), gs_gbc.predict(X_val))

    X_train, X_val, y_train, y_val = train_test_split(X_reg, y_reg,
                                                      random_state=42)
    params = {'learning_rate': [0.1, 0.3], 'max_depth': [3, 4]}
    gbrcv = GradientBoostingRegressorCV(score_precision=4, random_state=42,
                                        **params)
    gbrcv.fit(X_train, y_train)

    gs_gbr = GridSearchCV(GradientBoostingRegressor(random_state=42),
                          param_grid=params)
    gs_gbr.fit(X_train, y_train)

    assert_equal(gs_gbr.best_estimator_.n_estimators, 100)
    # FIXME: The below test depends on numpy/scipy versions
    # assert_equal(gbrcv.best_estimator_.n_estimators, 56)

    pred_1 = gs_gbr.best_estimator_.predict(X_val)
    pred_2 = gbrcv.best_estimator_.predict(X_val)

    # Make sure the MSE between the two best models is less than 0.25
    mse = ((pred_1 - pred_2) ** 2).mean()
    assert_less(0.25, mse)
