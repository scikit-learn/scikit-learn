import numpy as np

from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble.gradient_boosting import GradientBoostingClassifier


def bench_random_gaussian():
    rs = np.random.RandomState(13)
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

    gbrt = GradientBoostingClassifier(n_iter=100, max_depth=1)
    gbrt.fit(X_train, y_train)
    error_rate = (1.0 - gbrt.score(X_test, y_test)) * 100.0
    print "GBRT: %.2f" % error_rate

if __name__ == "__main__":
    bench_random_gaussian()
