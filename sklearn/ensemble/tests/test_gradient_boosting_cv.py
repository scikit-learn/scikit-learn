import numpy as np

from sklearn import datasets
from sklearn.ensemble import GradientBoostingClassifierCV
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.grid_search import GridSearchCV
from sklearn.metrics.cluster import homogeneity_score


rng = np.random.RandomState(0)
iris = datasets.load_iris()
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


def test_no_param():
    gbcv = GradientBoostingClassifierCV(max_iterations=1000)
    gbcv.fit(iris.data, iris.target)
    gbcv.best_params_['n_estimators'] < 1000


def test_single_param():
    gbcv = GradientBoostingClassifierCV(max_iterations=1000, learning_rate=0.1,
                                        max_depth=3)
    gbcv.fit(iris.data, iris.target)
    gbcv.best_params_['n_estimators'] < 1000


def test_param_grid():
    gbcv = GradientBoostingClassifierCV(max_iterations=1000,
                                        learning_rate=[0.1, 0.3],
                                        max_depth=[3, 4])
    gbcv.fit(iris.data, iris.target)
    gbcv.best_params_['n_estimators'] < 1000


def test_predict():
    X_train = iris.data[:100, :]
    y_train = iris.target[:100]

    X_validation = iris.data[100:, :]
    y_validation = iris.target[100:]

    params = {'learning_rate': [0.1, 0.3], 'max_depth': [3, 4]}
    gbcv = GradientBoostingClassifierCV(**params)
    gbcv.fit(X_train, y_train)

    gscv = GridSearchCV(GradientBoostingClassifier(), params)
    gscv.fit(X_train, y_train)

    l1 = gbcv.predict(X_test)
    l2 = gscv.predict(X_test)
    assert homogeneity_score(l1, l2) == 1.0
