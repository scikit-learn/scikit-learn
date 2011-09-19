import numpy as np

from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.utils import check_random_state

from sklearn.tree import DecisionTreeClassifier, DecisionTreeRegressor
from sklearn.ensemble.gradient_boosting import GradientBoostingClassifier
from sklearn.ensemble.gradient_boosting import GradientBoostingRegressor


def repeat(f):
    def wrapper(*args, **kargs):
        scores = []
        for i in range(10):
            scores.append(f(*args, random_state=i, **kargs))
        return np.mean(scores)
    return wrapper


def bench_random_gaussian(random_state=None):
    rs = check_random_state(random_state)
    shape = (12000, 10)
    X = rs.normal(size=shape).reshape(shape)
    y = ((X ** 2.0).sum(axis=1) > 9.34).astype(np.float64)
    y[y == 0.0] = -1.0

    X_train, X_test = X[:2000], X[2000:]
    y_train, y_test = y[:2000], y[2000:]

    stump = DecisionTreeClassifier(max_depth=1)
    stump.fit(X_train, y_train)
    error_rate = (1.0 - stump.score(X_test, y_test)) * 100.0

    print "Stump: %.2f" % error_rate

    tree = DecisionTreeClassifier(max_depth=20, min_split=5)
    tree.fit(X_train, y_train)
    error_rate = (1.0 - tree.score(X_test, y_test)) * 100.0

    print "Tree: %.2f" % error_rate

    tree = DecisionTreeRegressor(max_depth=20, min_split=5)
    tree.fit(X_train, y_train)
    y_pred = tree.predict(X_test)
    y_pred[y_pred > 0.0] = 1.0
    y_pred[y_pred <= 0.0] = -1.0
    error_rate = (1.0 - np.mean(y_pred == y_test)) * 100.0

    print "Tree: %.2f" % error_rate

    gbrt = GradientBoostingClassifier(loss='deviance', n_iter=100,
                                      min_split=5, max_depth=1,
                                      learn_rate=0.1)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test)) * 100.0
    print "GBRT: %.2f" % error_rate


regression_params = {'n_iter': 100, 'max_depth': 4,
                     'min_split': 1, 'learn_rate': 0.1,
                     'loss': 'ls'}


@repeat
def bench_boston(random_state=None):
    boston = datasets.load_boston()
    X, y = shuffle(boston.data, boston.target, random_state=random_state)
    offset = int(X.shape[0] * 0.9)
    X_train = X[:offset]
    y_train = y[:offset]
    X_test = X[offset:]
    y_test = y[offset:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    return mse


@repeat
def bench_friedman1(random_state=None):
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=random_state, noise=1.0)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    return mse


@repeat
def bench_friedman2(random_state=None):
    X, y = datasets.make_friedman2(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    return mse


@repeat
def bench_friedman3(random_state=None):
    X, y = datasets.make_friedman3(n_samples=1200, random_state=random_state)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]
    clf = GradientBoostingRegressor(**regression_params)
    clf.fit(X_train, y_train)
    mse = np.mean((clf.predict(X_test) - y_test) ** 2.0)
    return mse


if __name__ == "__main__":
    print "Example 10.2"
    bench_random_gaussian()

    #print "Boston", bench_boston()
    #print "Friedman#1", bench_friedman1()
    #print "Friedman#2", bench_friedman2()
    #print "Friedman#3", bench_friedman3()
